#!/usr/bin/env bash
# Call RAW variants for G60 output only (no filtering, no deletions).
# Produces per-gene raw VCFs and a merged raw VCF for ParaDISM,
# and a raw VCF for the base aligner.
#
# Run from anywhere: bash giab_benchmark/call_variants_raw_g60.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

REFERENCE="$SCRIPT_DIR/ref.fa"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate paradism_env

GENES=($(grep "^>" "$REFERENCE" | sed 's/^>//'))

call_raw_from_gene_bams() {
    local bam_dir=$1
    local out_dir=$2
    local prefix=$3

    echo "Calling RAW variants per-gene: $bam_dir"
    mkdir -p "$out_dir/per_gene"

    local gene_vcfs=()
    for gene in "${GENES[@]}"; do
        local gene_bam="${bam_dir}/${prefix}_${gene}.sorted.bam"
        local gene_vcf="${out_dir}/per_gene/${gene}.raw.vcf"
        local gene_vcfgz="${gene_vcf}.gz"
        if [[ -f "$gene_bam" ]]; then
            echo "  FreeBayes: $gene"
            freebayes --bam "$gene_bam" --fasta-reference "$REFERENCE" \
                --ploidy 2 --min-alternate-count 5 > "$gene_vcf"
            bgzip -f -c "$gene_vcf" > "$gene_vcfgz"
            bcftools index -f "$gene_vcfgz"
            gene_vcfs+=("$gene_vcfgz")
        fi
    done

    if [[ ${#gene_vcfs[@]} -eq 0 ]]; then
        echo "  No gene BAMs found; skipping merge"
        return
    fi

    echo "  Merging ${#gene_vcfs[@]} per-gene VCFs..."
    bcftools concat -a -O v -o "${out_dir}/variants_raw.vcf" "${gene_vcfs[@]}"
    bgzip -f -c "${out_dir}/variants_raw.vcf" > "${out_dir}/variants_raw.vcf.gz"
    bcftools index -f "${out_dir}/variants_raw.vcf.gz"
}

call_raw_from_sam() {
    local sam_file=$1
    local out_dir=$2

    echo "Calling RAW variants from original SAM: $sam_file"
    mkdir -p "$out_dir"

    if [[ ! -f "$sam_file" ]]; then
        echo "  SAM file not found: $sam_file"
        return
    fi

    local sorted_bam="${out_dir}/mapped_reads.sorted.bam"
    if [[ ! -f "$sorted_bam" ]]; then
        echo "  Converting SAM to sorted BAM..."
        samtools view -bS "$sam_file" | samtools sort -@ 8 -o "$sorted_bam"
        samtools index "$sorted_bam"
    fi

    echo "  FreeBayes (raw)..."
    freebayes --bam "$sorted_bam" --fasta-reference "$REFERENCE" \
        --ploidy 2 --min-alternate-count 5 > "${out_dir}/variants_raw.vcf"

    bgzip -f -c "${out_dir}/variants_raw.vcf" > "${out_dir}/variants_raw.vcf.gz"
    bcftools index -f "${out_dir}/variants_raw.vcf.gz"
}

echo "=========================================="
echo "RAW Variant Calling (G60 only)"
echo "=========================================="
echo ""

G60_DIR="$SCRIPT_DIR/giab_hg002_output_bowtie2_G60_min5_qfilters"
G60_PREFIX="giab_hg002_output_bowtie2_G60_min5_qfilters"

if [[ -d "$G60_DIR/final_outputs" ]]; then
    echo "=== G60 ParaDISM (raw per-gene) ==="
    call_raw_from_gene_bams "$G60_DIR/final_outputs/${G60_PREFIX}_bam" \
                            "$G60_DIR/final_outputs/variant_calling_raw" \
                            "$G60_PREFIX"
fi

if [[ -f "$G60_DIR/iteration_1/mapped_reads.sam" ]]; then
    echo "=== G60 Base Aligner (raw) ==="
    call_raw_from_sam "$G60_DIR/iteration_1/mapped_reads.sam" \
                      "$G60_DIR/iteration_1/variant_calling_basealigner_raw"
fi

echo ""
echo "Done!"
