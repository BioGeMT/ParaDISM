#!/bin/bash
# Download and prepare GIAB HG002 ground truth VCF for PKD1 region
# Run from anywhere: bash giab_benchmark/prepare_giab_truth.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/giab_hg002_vcf"
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

REGIONS="chr16:2084108-2132498,chr16:16309741-16334790,chr16:16355624-16378107,chr16:14910951-14936308,chr16:18329800-18349076,chr16:18369921-18398540,chr16:15120642-15151164"

echo ""
echo "Extracting PKD1 + pseudogene variants (${REGIONS})..."

OUTPUT_VCF="$OUTPUT_DIR/HG002_PKD1_genes_SNPs_exact.vcf.gz"

bcftools view -r "$REGIONS" "$BENCHMARK_VCF" | \
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
