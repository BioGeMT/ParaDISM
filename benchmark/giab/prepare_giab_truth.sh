#!/bin/bash
# Download and prepare GIAB HG002 ground truth VCF for PKD1 region
# Run from anywhere: bash benchmark/prepare_giab_truth.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
cd "$PROJECT_ROOT"

OUTPUT_DIR="giab_hg002_vcf"
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Preparing GIAB HG002 Ground Truth VCF"
echo "=========================================="
echo ""

BENCHMARK_VCF="$OUTPUT_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
if [[ ! -f "$BENCHMARK_VCF" ]]; then
    echo "Downloading GIAB HG002 benchmark VCF..."
    wget -P "$OUTPUT_DIR" \
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    wget -P "$OUTPUT_DIR" \
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
else
    echo "Benchmark VCF already exists"
fi

PKD1_START=2138711
PKD1_END=2185899

echo ""
echo "Extracting PKD1 region variants (chr16:${PKD1_START}-${PKD1_END})..."

OUTPUT_VCF="$OUTPUT_DIR/HG002_PKD1_genes_SNPs_illumina_only.vcf.gz"

bcftools view -r "chr16:${PKD1_START}-${PKD1_END}" "$BENCHMARK_VCF" | \
    bcftools view -v snps | \
    bgzip -c > "$OUTPUT_VCF"

bcftools index "$OUTPUT_VCF"

VARIANT_COUNT=$(bcftools view -H "$OUTPUT_VCF" | wc -l)

echo ""
echo "=========================================="
echo "Ground Truth VCF Ready"
echo "=========================================="
echo "Output: $OUTPUT_VCF"
echo "Variants: $VARIANT_COUNT SNPs"
