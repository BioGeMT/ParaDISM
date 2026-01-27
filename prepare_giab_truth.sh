#!/bin/bash
# Download and prepare GIAB HG002 ground truth VCF for PKD1 region
# Creates filtered VCF with Illumina-supported SNPs only

set -euo pipefail

OUTPUT_DIR="giab_hg002_vcf"
mkdir -p "$OUTPUT_DIR"

echo "=========================================="
echo "Preparing GIAB HG002 Ground Truth VCF"
echo "=========================================="
echo ""

# Download GIAB benchmark VCF if not present
BENCHMARK_VCF="$OUTPUT_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
if [[ ! -f "$BENCHMARK_VCF" ]]; then
    echo "Downloading GIAB HG002 benchmark VCF..."
    wget -P "$OUTPUT_DIR" \
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    wget -P "$OUTPUT_DIR" \
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
    echo "  Done"
else
    echo "Benchmark VCF already exists"
fi

# PKD1 gene region coordinates (chr16)
# PKD1: 2138711-2185899 (minus strand)
PKD1_START=2138711
PKD1_END=2185899

echo ""
echo "Extracting PKD1 region variants (chr16:${PKD1_START}-${PKD1_END})..."

# Extract PKD1 region SNPs
OUTPUT_VCF="$OUTPUT_DIR/HG002_PKD1_genes_SNPs_illumina_only.vcf.gz"

bcftools view -r "chr16:${PKD1_START}-${PKD1_END}" "$BENCHMARK_VCF" | \
    bcftools view -v snps | \
    bgzip -c > "$OUTPUT_VCF"

bcftools index "$OUTPUT_VCF"

# Count variants
VARIANT_COUNT=$(bcftools view -H "$OUTPUT_VCF" | wc -l)

echo ""
echo "=========================================="
echo "Ground Truth VCF Ready"
echo "=========================================="
echo ""
echo "Output: $OUTPUT_VCF"
echo "Variants: $VARIANT_COUNT SNPs"
echo ""
echo "Use this for comparison with ParaDISM variant calls"
