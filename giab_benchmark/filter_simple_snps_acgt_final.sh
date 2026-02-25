#!/usr/bin/env bash
# Create final filtered VCF from a simple SNP A/C/G/T VCF.
#
# Default filter:
#   INFO/AO>=3 && QUAL>=10 && INFO/SAF>=1 && INFO/SAR>=1
#
# Usage:
#   bash giab_benchmark/filter_simple_snps_acgt_final.sh \
#     --input-vcf <in.vcf.gz> \
#     --output-vcf <out.vcf.gz> \
#     --benchmark-bed <mask.bed>
#
# Notes:
# - If --benchmark-bed is provided, records are first restricted to BED regions.
# - Output is bgzipped and indexed (.csi).

set -euo pipefail

INPUT_VCF=""
OUTPUT_VCF=""
BENCHMARK_BED=""
FILTER_EXPR='INFO/AO>=3 && QUAL>=10 && INFO/SAF>=1 && INFO/SAR>=1'

usage() {
  cat <<'EOF'
Usage:
  bash giab_benchmark/filter_simple_snps_acgt_final.sh \
    --input-vcf <in.vcf.gz> \
    --output-vcf <out.vcf.gz> \
    --benchmark-bed <mask.bed>

Options:
  --input-vcf FILE      Input VCF.GZ (simple SNP A/C/G/T set)
  --output-vcf FILE     Output VCF.GZ
  --benchmark-bed FILE  BED mask to restrict to confident regions (required)
  -h, --help            Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-vcf)
      INPUT_VCF="$2"
      shift 2
      ;;
    --output-vcf)
      OUTPUT_VCF="$2"
      shift 2
      ;;
    --benchmark-bed)
      BENCHMARK_BED="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$INPUT_VCF" || -z "$OUTPUT_VCF" || -z "$BENCHMARK_BED" ]]; then
  echo "Error: --input-vcf, --output-vcf, and --benchmark-bed are required." >&2
  usage >&2
  exit 1
fi

if [[ ! -f "$INPUT_VCF" ]]; then
  echo "Error: input VCF not found: $INPUT_VCF" >&2
  exit 1
fi

if [[ ! -f "$BENCHMARK_BED" ]]; then
  echo "Error: benchmark BED not found: $BENCHMARK_BED" >&2
  exit 1
fi

BCFTOOLS="$(command -v bcftools || true)"
if [[ -z "$BCFTOOLS" && -x "$HOME/miniconda3/envs/paradism_env/bin/bcftools" ]]; then
  BCFTOOLS="$HOME/miniconda3/envs/paradism_env/bin/bcftools"
fi
if [[ -z "$BCFTOOLS" ]]; then
  echo "Error: bcftools not found in PATH and not found in paradism_env." >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_VCF")"

echo "Input:  $INPUT_VCF"
echo "Output: $OUTPUT_VCF"
echo "Filter: $FILTER_EXPR"
echo "BED:    $BENCHMARK_BED"

"$BCFTOOLS" view -R "$BENCHMARK_BED" -i "$FILTER_EXPR" -Oz -o "$OUTPUT_VCF" "$INPUT_VCF"

"$BCFTOOLS" index -f "$OUTPUT_VCF"

count="$("$BCFTOOLS" view -H "$OUTPUT_VCF" | wc -l)"
echo "Done. Variants: $count"
