#!/usr/bin/env bash
# Filter RAW VCFs for G60 output only (no deletions).
# Filters: SNPs only, GT=1/1 (no DP/MQM filters).
#
# Run from anywhere: bash giab_benchmark/filter_variants_g60.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate paradism_env

filter_vcf() {
    local raw_vcf=$1
    local out_dir=$2

    if [[ ! -f "$raw_vcf" ]]; then
        echo "  Raw VCF not found: $raw_vcf"
        return
    fi

    mkdir -p "$out_dir"
    echo "  Filtering: $raw_vcf"

    local out_vcf="${out_dir}/variants.vcf"
    bcftools view -v snps "$raw_vcf" | \
        bcftools filter -i "GT=\"1/1\"" \
        > "$out_vcf"

    bgzip -f -c "$out_vcf" > "${out_vcf}.gz"
    bcftools index -f "${out_vcf}.gz"
}

echo "=========================================="
echo "Filter Variants (G60 only)"
echo "=========================================="
echo ""

G60_DIR="$SCRIPT_DIR/giab_hg002_output_bowtie2_G60_min5_qfilters"

echo "=== G60 ParaDISM (final) ==="
filter_vcf "$G60_DIR/final_outputs/variant_calling_raw/variants_raw.vcf" \
           "$G60_DIR/final_outputs/variant_calling_balanced"

echo "=== G60 Base Aligner (iteration_1) ==="
filter_vcf "$G60_DIR/iteration_1/variant_calling_basealigner_raw/variants_raw.vcf" \
           "$G60_DIR/iteration_1/variant_calling_basealigner_balanced"

echo ""
echo "Done!"
