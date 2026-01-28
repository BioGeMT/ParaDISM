#!/usr/bin/env bash
# Call variants with BALANCED filters for G30 and G60 outputs
# Filters: QUAL>=500, DP>=30, AF>=0.90, RO<=5, MQMR<=10
# Run from anywhere: bash giab_benchmark/call_variants_balanced_filters.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

REFERENCE="ref.fa"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate paradism_env

GENES=($(grep "^>" "$REFERENCE" | sed 's/^>//'))

call_variants() {
    local bam_dir=$1
    local variant_dir=$2
    local prefix=$3

    echo "Calling variants: $bam_dir"
    mkdir -p "${variant_dir}/per_gene_vcfs"
    local per_gene_vcfs=()

    for gene in "${GENES[@]}"; do
        local gene_bam="${bam_dir}/${prefix}_${gene}.sorted.bam"
        [[ ! -f "$gene_bam" ]] && continue

        local gene_vcf="${variant_dir}/per_gene_vcfs/${gene}.vcf"
        local gene_vcf_tmp="${variant_dir}/per_gene_vcfs/${gene}.tmp.vcf"

        freebayes --bam "$gene_bam" --fasta-reference "$REFERENCE" \
            --ploidy 2 --min-alternate-count 5 --region "$gene" -i 2>/dev/null | \
            bcftools view -v snps > "$gene_vcf_tmp" 2>/dev/null

        if [[ -s "$gene_vcf_tmp" ]]; then
            bcftools filter -e "QUAL<500 || INFO/DP<30 || INFO/AF<0.90 || INFO/RO>5 || INFO/MQMR>10" \
                "$gene_vcf_tmp" > "$gene_vcf" 2>/dev/null || cp "$gene_vcf_tmp" "$gene_vcf"
            rm -f "$gene_vcf_tmp"
        else
            rm -f "$gene_vcf_tmp"
        fi

        if [[ -f "$gene_vcf" ]] && [[ -s "$gene_vcf" ]]; then
            variant_count=$(grep -v "^#" "$gene_vcf" | wc -l)
            if [[ $variant_count -gt 0 ]]; then
                per_gene_vcfs+=("$gene_vcf")
            else
                rm -f "$gene_vcf"
            fi
        else
            rm -f "$gene_vcf"
        fi
    done

    if [[ ${#per_gene_vcfs[@]} -eq 0 ]]; then
        echo "##fileformat=VCFv4.2" > "${variant_dir}/variants.vcf"
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> "${variant_dir}/variants.vcf"
        return
    fi

    if [[ ${#per_gene_vcfs[@]} -eq 1 ]]; then
        cp "${per_gene_vcfs[0]}" "${variant_dir}/variants.vcf"
    else
        for vcf in "${per_gene_vcfs[@]}"; do
            bgzip -c "$vcf" > "${vcf}.gz"
            bcftools index "${vcf}.gz"
        done
        bcftools merge -m all --force-samples "${per_gene_vcfs[@]/%/.gz}" | bcftools sort > "${variant_dir}/variants.vcf"
    fi

    bgzip -c "${variant_dir}/variants.vcf" > "${variant_dir}/variants.vcf.gz"
    bcftools index "${variant_dir}/variants.vcf.gz"

    echo "  Variants: $(grep -v '^#' "${variant_dir}/variants.vcf" | wc -l)"
}

echo "=========================================="
echo "Variant Calling with Balanced Filters"
echo "=========================================="
echo ""

G60_DIR="giab_benchmark/giab_hg002_output_bowtie2_G60_min5_qfilters"
G60_PREFIX="giab_hg002_output_bowtie2_G60_min5_qfilters"
if [[ -d "$G60_DIR/final_outputs" ]]; then
    echo "Processing G60 ParaDISM (final)..."
    call_variants "$G60_DIR/final_outputs/${G60_PREFIX}_bam" \
                  "$G60_DIR/final_outputs/variant_calling_balanced" \
                  "$G60_PREFIX"
fi
if [[ -d "$G60_DIR/iteration_1" ]]; then
    echo "Processing G60 Base Aligner (iteration_1)..."
    call_variants "$G60_DIR/iteration_1/${G60_PREFIX}_bam" \
                  "$G60_DIR/iteration_1/variant_calling_balanced" \
                  "$G60_PREFIX"
fi

G40_DIR="giab_benchmark/giab_hg002_output_bowtie2_G40_min5_qfilters"
G40_PREFIX="giab_hg002_output_bowtie2_G40_min5_qfilters"
if [[ -d "$G40_DIR/final_outputs" ]]; then
    echo "Processing G40 ParaDISM (final)..."
    call_variants "$G40_DIR/final_outputs/${G40_PREFIX}_bam" \
                  "$G40_DIR/final_outputs/variant_calling_balanced" \
                  "$G40_PREFIX"
fi
if [[ -d "$G40_DIR/iteration_1" ]]; then
    echo "Processing G40 Base Aligner (iteration_1)..."
    call_variants "$G40_DIR/iteration_1/${G40_PREFIX}_bam" \
                  "$G40_DIR/iteration_1/variant_calling_balanced" \
                  "$G40_PREFIX"
fi

echo ""
echo "Done!"
