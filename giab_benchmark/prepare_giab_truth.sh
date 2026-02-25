#!/bin/bash
# Download and prepare GIAB HG002 ground truth VCF for PKD1 region
# Run from anywhere: bash giab_benchmark/prepare_giab_truth.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/giab_hg002_vcf"
mkdir -p "$OUTPUT_DIR"

BCFTOOLS="$(command -v bcftools || true)"
BGZIP="$(command -v bgzip || true)"
ATOMIZE_SCRIPT="$SCRIPT_DIR/atomize_equal_length_substitutions.py"
if [[ -z "$BCFTOOLS" || -z "$BGZIP" ]]; then
    if command -v conda >/dev/null 2>&1; then
        BCFTOOLS="$(conda run -n paradism_env which bcftools 2>/dev/null || true)"
        BGZIP="$(conda run -n paradism_env which bgzip 2>/dev/null || true)"
    fi
fi
if [[ -z "$BCFTOOLS" || -z "$BGZIP" ]]; then
    echo "ERROR: bcftools/bgzip not found on PATH (and conda env 'paradism_env' not available)." >&2
    exit 1
fi
if [[ ! -f "$ATOMIZE_SCRIPT" ]]; then
    echo "ERROR: atomize script not found: $ATOMIZE_SCRIPT" >&2
    exit 1
fi

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

# Regions match ParaDISM reference contigs (giab_benchmark/ref.fa) and
# the GENE_COORDS mappings in make_truth_gene_coords_vcf.py.
REGIONS="chr16:2088707-2135898,chr16:16310133-16334190,chr16:16356223-16377507,chr16:14911550-14935708,chr16:18334399-18352476,chr16:18374520-18402014,chr16:15125138-15154873"

echo ""
echo "Extracting PKD1 + pseudogene variants (${REGIONS})..."

OUTPUT_VCF="$OUTPUT_DIR/HG002_PKD1_genes_SNPs_exact.vcf.gz"

$BCFTOOLS view -r "$REGIONS" "$BENCHMARK_VCF" | \
    $BCFTOOLS norm -m -any -O v | \
    python3 "$ATOMIZE_SCRIPT" | \
    $BCFTOOLS view -v snps -m2 -M2 -i 'REF~"^[ACGT]$" && ALT~"^[ACGT]$"' | \
    $BGZIP -c > "$OUTPUT_VCF"

$BCFTOOLS index "$OUTPUT_VCF"

VARIANT_COUNT=$($BCFTOOLS view -H "$OUTPUT_VCF" | wc -l)

BENCHMARK_BED="$OUTPUT_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
BENCHMARKABLE_VCF="$OUTPUT_DIR/HG002_PKD1_genes_SNPs_exact_benchmarkable.vcf.gz"
$BCFTOOLS view -R "$BENCHMARK_BED" -Oz -o "$BENCHMARKABLE_VCF" "$OUTPUT_VCF"
$BCFTOOLS index "$BENCHMARKABLE_VCF"
BENCHMARKABLE_COUNT=$($BCFTOOLS view -H "$BENCHMARKABLE_VCF" | wc -l)

echo ""
echo "=========================================="
echo "Ground Truth VCF Ready"
echo "=========================================="
echo "Output: $OUTPUT_VCF"
echo "Variants: $VARIANT_COUNT SNPs"

echo ""
echo "Converting truth VCF to ParaDISM gene-contig coordinates..."
echo "Benchmarkable (BED-filtered): $BENCHMARKABLE_VCF"
echo "Benchmarkable variants: $BENCHMARKABLE_COUNT SNPs"
python3 "$SCRIPT_DIR/make_truth_gene_coords_vcf.py" --truth-vcf "$BENCHMARKABLE_VCF"
