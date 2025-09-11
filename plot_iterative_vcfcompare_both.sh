#!/bin/bash

# Run VCFCompare for both FreeBayes and bcftools iteratively (SNV-only plots)
# Usage: ./plot_iterative_vcfcompare_both.sh <iterations> <FREEBAYES_OUT> <BCFTOOLS_OUT> <truth_vcf> [vcfcompare_dir]

set -euo pipefail

ITERATIONS="${1:?iterations required}"
FREEBAYES_DIR="${2:?FREEBAYES output dir required}"
BCFTOOLS_DIR="${3:?BCFTOOLS output dir required}"
TRUTH_VCF="${4:?truth VCF required}"
VCFCOMPARE_DIR="${5:-./VCFCompare}"

echo "Running VCFCompare for FreeBayes..."
./plot_iterative_vcfcompare.sh freebayes "$ITERATIONS" "$FREEBAYES_DIR" "$TRUTH_VCF" "$VCFCOMPARE_DIR"

echo "Running VCFCompare for bcftools..."
./plot_iterative_vcfcompare.sh bcftools "$ITERATIONS" "$BCFTOOLS_DIR" "$TRUTH_VCF" "$VCFCOMPARE_DIR"

echo "Done. Outputs:"
echo "  - plots_iterative_vcfcompare_freebayes/"
echo "  - plots_iterative_vcfcompare_bcftools/"

