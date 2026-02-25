#!/usr/bin/env bash
# Run GIAB post-processing from variant calling to final vcf_out plots/metrics.
#
# This script assumes ParaDISM alignment/iterations are already finished in --run-dir.
#
# It performs:
#   1) call_variants_raw_g60.sh
#   2) sort/index simple SNP VCFs
#   3) final benchmark filtering
#   4) plot_balanced_comparison.py
#   5) plot_coverage_split_confusion.py

set -euo pipefail

CALLER_CWD="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RUN_DIR=""
OUT_DIR=""
THREADS=8
SKIP_CALL=0
BENCHMARK_BED="$SCRIPT_DIR/giab_hg002_vcf/HG002_PKD1_genes_benchmarkable_gene_coords.bed"

CALL_SCRIPT="$SCRIPT_DIR/call_variants_raw_g60.sh"
FILTER_SCRIPT="$SCRIPT_DIR/filter_simple_snps_acgt_final.sh"
PLOT_SCRIPT="$SCRIPT_DIR/plot_balanced_comparison.py"
PLOT_COV_SCRIPT="$SCRIPT_DIR/plot_coverage_split_confusion.py"

usage() {
    cat <<'EOF'
Usage:
  bash giab_benchmark/run_variant_calling_to_vcf_out.sh \
    --run-dir <paradism_run_dir> \
    --out-dir <vcf_out_dir> \
    [options]

Required:
  --run-dir DIR          ParaDISM run directory (already completed)
  --out-dir DIR          Output directory for plots/metrics (vcf_out*)

Optional:
  --threads N            Threads for variant calling sort steps (default: 8)
  --benchmark-bed FILE   Benchmark BED in gene coordinates
                         (default: giab_benchmark/giab_hg002_vcf/HG002_PKD1_genes_benchmarkable_gene_coords.bed)
  --skip-call            Skip call_variants_raw_g60.sh and only do sort/filter/plots
  -h, --help             Show help
EOF
}

to_abs_path() {
    local path="$1"
    if [[ "$path" == /* ]]; then
        printf "%s\n" "$path"
    else
        printf "%s\n" "$CALLER_CWD/$path"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-dir)
            RUN_DIR="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --benchmark-bed)
            BENCHMARK_BED="$2"
            shift 2
            ;;
        --skip-call)
            SKIP_CALL=1
            shift
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

if [[ -z "$RUN_DIR" || -z "$OUT_DIR" ]]; then
    echo "Error: --run-dir and --out-dir are required." >&2
    usage >&2
    exit 1
fi

RUN_DIR="$(to_abs_path "$RUN_DIR")"
OUT_DIR="$(to_abs_path "$OUT_DIR")"
BENCHMARK_BED="$(to_abs_path "$BENCHMARK_BED")"

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS <= 0 )); then
    echo "Error: --threads must be a positive integer (got: $THREADS)" >&2
    exit 1
fi

if [[ ! -d "$RUN_DIR" ]]; then
    echo "Error: run directory not found: $RUN_DIR" >&2
    exit 1
fi

for required in "$CALL_SCRIPT" "$FILTER_SCRIPT" "$PLOT_SCRIPT" "$PLOT_COV_SCRIPT"; do
    if [[ ! -f "$required" ]]; then
        echo "Error: required script not found: $required" >&2
        exit 1
    fi
done

if [[ ! -f "$BENCHMARK_BED" ]]; then
    echo "Error: benchmark BED not found: $BENCHMARK_BED" >&2
    exit 1
fi

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate paradism_env

echo "=========================================="
echo "GIAB post-processing (variant_calling -> vcf_out)"
echo "=========================================="
echo "Run dir:       $RUN_DIR"
echo "Out dir:       $OUT_DIR"
echo "Threads:       $THREADS"
echo "Benchmark BED: $BENCHMARK_BED"
echo ""

if [[ $SKIP_CALL -eq 0 ]]; then
    echo "1) Calling raw + simple SNP A/C/G/T variants..."
    bash "$CALL_SCRIPT" \
        --input-dir "$RUN_DIR" \
        --output-dir "$RUN_DIR/variant_calling" \
        --threads "$THREADS"
else
    echo "1) Skipping raw variant calling (--skip-call)"
fi

PARADISM_SIMPLE="$RUN_DIR/variant_calling/paradism_raw/variants_simple_snps_acgt.vcf.gz"
BASE_SIMPLE="$RUN_DIR/variant_calling/basealigner_raw/variants_simple_snps_acgt.vcf.gz"
PARADISM_SORTED="$RUN_DIR/variant_calling/paradism_raw/variants_simple_snps_acgt.sorted.vcf.gz"
BASE_SORTED="$RUN_DIR/variant_calling/basealigner_raw/variants_simple_snps_acgt.sorted.vcf.gz"
PARADISM_FINAL="$RUN_DIR/variant_calling/paradism_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
BASE_FINAL="$RUN_DIR/variant_calling/basealigner_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"

if [[ ! -f "$PARADISM_SIMPLE" || ! -f "$BASE_SIMPLE" ]]; then
    echo "Error: simple SNP VCFs not found after variant calling." >&2
    echo "  Missing one of:" >&2
    echo "    $PARADISM_SIMPLE" >&2
    echo "    $BASE_SIMPLE" >&2
    exit 1
fi

echo ""
echo "2) Sorting/indexing simple SNP VCFs..."
bcftools sort -Oz -o "$PARADISM_SORTED" "$PARADISM_SIMPLE"
bcftools index -f "$PARADISM_SORTED"
bcftools sort -Oz -o "$BASE_SORTED" "$BASE_SIMPLE"
bcftools index -f "$BASE_SORTED"

echo ""
echo "3) Applying final benchmark filters..."
bash "$FILTER_SCRIPT" \
    --input-vcf "$PARADISM_SORTED" \
    --output-vcf "$PARADISM_FINAL" \
    --benchmark-bed "$BENCHMARK_BED"

bash "$FILTER_SCRIPT" \
    --input-vcf "$BASE_SORTED" \
    --output-vcf "$BASE_FINAL" \
    --benchmark-bed "$BENCHMARK_BED"

echo ""
echo "4) Generating comparison plots/metrics..."
python "$PLOT_SCRIPT" \
    --dataset-dir "$RUN_DIR" \
    --out-dir "$OUT_DIR"

echo ""
echo "5) Generating coverage-split confusion matrices..."
python "$PLOT_COV_SCRIPT" \
    --dataset-dir "$RUN_DIR" \
    --out-dir "$OUT_DIR"

echo ""
echo "Done."
echo "Final ParaDISM VCF: $PARADISM_FINAL"
echo "Final Base VCF:     $BASE_FINAL"
echo "Plots/metrics dir:  $OUT_DIR"
