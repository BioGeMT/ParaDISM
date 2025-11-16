#!/bin/bash

set -euo pipefail

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
SEED_START=1
SEED_END=1000
SIM_OUTPUT_BASE="sim_output"
REFERENCE="ref.fa"
THREADS=16
NUM_READS=100000
ERROR_RATE=0.01             # default sequencing error rate (per base)
MINIMAP2_PROFILE="short"     # sr preset
SNP_PROB=0.001
INDEL_PROB=0.0001
FREEBAYES_PLOIDY=1
FREEBAYES_MIN_ALT=5

ALIGNERS=("bowtie2" "bwa-mem2" "minimap2")

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
format_error_suffix() {
    local rate="$1"
    local formatted
    formatted=$(printf "%.3f" "$rate")
    echo "${formatted/./_}"
}

run_mapper() {
    local aligner="$1"
    local r1="$2"
    local r2="$3"
    local output_dir="$4"
    local prefix="$5"

    local cmd=(python mapper.py --read1 "$r1" --read2 "$r2" --reference "$REFERENCE" --aligner "$aligner" --threads "$THREADS" --output-dir "$output_dir" --prefix "$prefix")
    if [[ "$aligner" == "minimap2" ]]; then
        cmd+=(--minimap2-profile "$MINIMAP2_PROFILE")
    fi

    "${cmd[@]}"
}

run_base_alignment() {
    local aligner="$1"
    local r1="$2"
    local r2="$3"
    local base_dir="$4"
    local bam_path="$5"

    mkdir -p "$base_dir"
    local tmp_sam="$base_dir/base.tmp.sam"

    case "$aligner" in
        bowtie2)
            local index_prefix="$base_dir/ref_index"
            bowtie2-build "$REFERENCE" "$index_prefix"
            bowtie2 -p "$THREADS" -x "$index_prefix" -1 "$r1" -2 "$r2" -S "$tmp_sam"
            ;;
        bwa-mem2)
            local index_prefix="$base_dir/ref_index"
            bwa-mem2 index -p "$index_prefix" "$REFERENCE"
            bwa-mem2 mem -t "$THREADS" "$index_prefix" "$r1" "$r2" > "$tmp_sam"
            ;;
        minimap2)
            local preset
            case "$MINIMAP2_PROFILE" in
                short) preset="sr" ;;
                pacbio-hifi) preset="map-hifi" ;;
                pacbio-clr) preset="map-pb" ;;
                ont-q20) preset="lr:hq" ;;
                ont-standard) preset="map-ont" ;;
                *) preset="sr" ;;
            esac
            local index_mmi="$base_dir/ref_index.mmi"
            minimap2 -x "$preset" -d "$index_mmi" "$REFERENCE"
            minimap2 -ax "$preset" --MD -t "$THREADS" "$index_mmi" "$r1" "$r2" > "$tmp_sam"
            ;;
        *)
            echo "Unsupported aligner: $aligner" >&2
            exit 1
            ;;
    esac

    samtools view -b "$tmp_sam" | samtools sort -o "$bam_path"
    samtools index "$bam_path"
    rm -f "$tmp_sam"
    rm -rf "$base_dir"
}

run_variant_calling() {
    local aligner="$1"
    local seed="$2"
    local prefix="$3"
    local paradism_output="$4"
    local base_bam="$5"
    local aligner_dir="$6"
    local seed_dir="$7"

    local bam_dir="${paradism_output}/${prefix}_bam"
    if [[ ! -d "$bam_dir" ]]; then
        echo "Per-gene BAM directory not found: $bam_dir" >&2
        return 1
    fi

    local vc_dir="$aligner_dir/variant_calling_analysis"
    local per_gene_dir="$vc_dir/per_gene_vcfs"
    mkdir -p "$per_gene_dir"

    local combined_vcf="$vc_dir/${prefix}_paradism_combined.vcf"
    local direct_vcf="$vc_dir/${prefix}_${aligner}_base.vcf"

    shopt -s nullglob
    local bam_files=("$bam_dir"/*.sorted.bam)
    if (( ${#bam_files[@]} == 0 )); then
        shopt -u nullglob
        echo "No per-gene BAMs in $bam_dir" >&2
        return 1
    fi

    echo "  Calling variants for ParaDISM per-gene BAMs..."
    for bam_file in "${bam_files[@]}"; do
        local filename="${bam_file##*/}"
        local gene="${filename%.sorted.bam}"
        gene="${gene#${prefix}_}"
        local tmp_vcf="$per_gene_dir/${gene}.tmp.vcf"
        local final_vcf="$per_gene_dir/${gene}.vcf"

        freebayes \
            --bam "$bam_file" \
            --fasta-reference "$REFERENCE" \
            --ploidy "$FREEBAYES_PLOIDY" \
            --min-alternate-count "$FREEBAYES_MIN_ALT" \
            --region "$gene" \
            -i \
            > "$tmp_vcf" 2> "$per_gene_dir/${gene}.log" || {
                echo "    Warning: FreeBayes failed for $gene (seed $seed, $aligner); creating empty VCF"
                printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$final_vcf"
                rm -f "$tmp_vcf"
                continue
            }

        awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {ref=$4; alt=$5; if(length(ref)==1 && length(alt)==1) print}' "$tmp_vcf" > "$final_vcf"
        rm -f "$tmp_vcf"
    done

    echo "  Building combined ParaDISM VCF..."
    local first_vcf
    first_vcf=$(printf "%s\n" "$per_gene_dir"/*.vcf | head -n1)
    if [[ -n "$first_vcf" ]]; then
        {
            grep '^#' "$first_vcf"
            for vcf in "$per_gene_dir"/*.vcf; do
                grep -v '^#' "$vcf" || true
            done | sort -k1,1 -k2,2n
        } > "$combined_vcf"
    else
        printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$combined_vcf"
    fi

    echo "  Calling variants from base alignment..."
    local tmp_direct="$vc_dir/${prefix}_${aligner}_base.tmp.vcf"
    freebayes \
        --bam "$base_bam" \
        --fasta-reference "$REFERENCE" \
        --ploidy "$FREEBAYES_PLOIDY" \
        --min-alternate-count "$FREEBAYES_MIN_ALT" \
        -i \
        > "$tmp_direct" 2> "$vc_dir/${prefix}_${aligner}_base.log" || {
            echo "    Warning: FreeBayes failed for base alignment (seed $seed, $aligner); creating empty VCF"
            printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$direct_vcf"
            rm -f "$tmp_direct"
        }
    if [[ -f "$tmp_direct" ]]; then
        awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {ref=$4; alt=$5; if(length(ref)==1 && length(alt)==1) print}' "$tmp_direct" > "$direct_vcf"
        rm -f "$tmp_direct"
    fi
    shopt -u nullglob

    echo "  Running variant calling analysis..."
    python variant_calling_analysis.py \
        --aligner "$aligner" \
        --ground-truth "$seed_dir/all_genes_mutations.vcf" \
        --combined-vcf "$combined_vcf" \
        --direct-vcf "$direct_vcf" \
        --positions-tsv "$seed_dir/all_genes_mutations.tsv" \
        --analysis-dir "$vc_dir" \
        --output-prefix "${prefix}" \
        > "$vc_dir/${prefix}_analysis.log" 2>&1
}

# ------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------
ERROR_SUFFIX=$(format_error_suffix "$ERROR_RATE")

for seed in $(seq "$SEED_START" "$SEED_END"); do
    echo "=============================="
    echo "Processing seed: $seed"
    echo "=============================="

    seed_dir="$SIM_OUTPUT_BASE/seed_${seed}"
    mkdir -p "$seed_dir"

    echo "Simulating reads..."
    python read_simulator.py \
        --reference "$REFERENCE" \
        --output-dir "$seed_dir" \
        --seed "$seed" \
        --num-reads "$NUM_READS" \
        --error-rate "$ERROR_RATE" \
        --snp-prob "$SNP_PROB" \
        --indel-prob "$INDEL_PROB"

    r1="$seed_dir/simulated_r1_err_${ERROR_SUFFIX}.fq"
    r2="$seed_dir/simulated_r2_err_${ERROR_SUFFIX}.fq"

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "Missing simulated FASTQs for seed $seed" >&2
        exit 1
    fi

    for aligner in "${ALIGNERS[@]}"; do
        echo "--- Running aligner: $aligner (seed $seed) ---"
        aligner_dir="$seed_dir/$aligner"
        paradism_output="$aligner_dir/paradism_output"
        prefix="seed_${seed}_${aligner}"
        paradism_prefix="paradism_${prefix}"
        mkdir -p "$paradism_output"

        run_mapper "$aligner" "$r1" "$r2" "$paradism_output" "$paradism_prefix"

        base_bam="$aligner_dir/${aligner}_base.sorted.bam"
        base_tmp="$aligner_dir/base_align"
        run_base_alignment "$aligner" "$r1" "$r2" "$base_tmp" "$base_bam"

        mapper_tsv="$paradism_output/${paradism_prefix}_unique_mappings.tsv"
        analysis_dir="$aligner_dir/read_mapping_analysis"
        output_prefix="seed_${seed}_${aligner}"

        python read_mapping_analysis.py \
            --aligner "$aligner" \
            --fastq-r1 "$r1" \
            --mapper-tsv "$mapper_tsv" \
            --direct-sam "$base_bam" \
            --analysis-dir "$analysis_dir" \
            --output-prefix "$output_prefix"

        run_variant_calling "$aligner" "$seed" "$paradism_prefix" "$paradism_output" "$base_bam" "$aligner_dir" "$seed_dir"
    done
done

# ------------------------------------------------------------------
# Aggregate results across all seeds
# ------------------------------------------------------------------
echo "=============================="
echo "Aggregating results across all seeds..."
echo "=============================="

python aggregate_results.py \
    --sim-output-base "$SIM_OUTPUT_BASE" \
    --seed-start "$SEED_START" \
    --seed-end "$SEED_END" \
    --aligners "${ALIGNERS[@]}" \
    --output-dir "${SIM_OUTPUT_BASE}/aggregated_results"

echo "Aggregation complete!"
