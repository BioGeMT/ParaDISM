#!/usr/bin/env bash
# Liftover IGV per-sample VCFs and exons BED from gene-local to chr16 coordinates.
#
# Example:
#   bash hts_analysis/liftover_igv_hts.sh \
#     --run-root hts_analysis/HTS_rerun_bowtie2_1ref_t4_p10 \
#     --positions PKD1_b38_pseudogene_positions.txt

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash hts_analysis/liftover_igv_hts.sh [options]

Options:
  --run-root PATH       HTS run root directory
  --positions PATH      Gene positions file (e.g., PKD1_b38_pseudogene_positions.txt)
  --out-dir PATH        Output directory (default: <run-root>/igv_per_sample_chr16)
  -h, --help            Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

RUN_ROOT=""
POSITIONS=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root)
            RUN_ROOT="$2"
            shift 2
            ;;
        --positions)
            POSITIONS="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$RUN_ROOT" || -z "$POSITIONS" ]]; then
    echo "Error: --run-root and --positions are required."
    usage
    exit 1
fi

if [[ ! -d "$RUN_ROOT" ]]; then
    echo "Error: run root not found: $RUN_ROOT"
    exit 1
fi

if [[ ! -f "$POSITIONS" ]]; then
    echo "Error: positions file not found: $POSITIONS"
    exit 1
fi

IGV_DIR="${RUN_ROOT}/igv_per_sample"
if [[ ! -d "$IGV_DIR" ]]; then
    echo "Error: igv_per_sample directory not found: $IGV_DIR"
    exit 1
fi

if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR="${RUN_ROOT}/igv_per_sample_chr16"
fi

echo "=========================================="
echo "Liftover IGV VCFs to chr16 coordinates"
echo "Run root:   $RUN_ROOT"
echo "Positions:  $POSITIONS"
echo "Output:     $OUT_DIR"
echo "=========================================="

mkdir -p "$OUT_DIR/clinical" "$OUT_DIR/control" "$OUT_DIR/reference_assets"

# Liftover exons BED
EXONS_BED="${IGV_DIR}/reference_assets/PKD1_exons_1ref.fa.bed"
if [[ -f "$EXONS_BED" ]]; then
    echo "Lifting exons BED..."
    python paradism.py liftover --bed "$EXONS_BED" --positions "$POSITIONS" \
        -o "$OUT_DIR/reference_assets/PKD1_exons_chr16.bed"
else
    echo "Warning: exons BED not found at $EXONS_BED, skipping"
fi

# Liftover VCFs
SUMMARY_TSV="${OUT_DIR}/liftover_summary.tsv"
echo -e "group\tsample\tstatus\toutput" > "$SUMMARY_TSV"

TOTAL=0
OK=0
FAIL=0

liftover_group() {
    local group="$1"
    local src_dir="${IGV_DIR}/${group}"

    if [[ ! -d "$src_dir" ]]; then
        echo "Warning: $src_dir not found, skipping"
        return
    fi

    for vcf in "$src_dir"/*.vcf; do
        [[ -f "$vcf" ]] || continue
        local sample
        sample="$(basename "$vcf" .vcf)"
        local out_vcf="${OUT_DIR}/${group}/${sample}.vcf"

        TOTAL=$((TOTAL + 1))

        if python paradism.py liftover --vcf "$vcf" --positions "$POSITIONS" -o "$out_vcf" 2>/dev/null; then
            OK=$((OK + 1))
            printf "%s\t%s\t%s\t%s\n" "$group" "$sample" "ok" "$out_vcf" >> "$SUMMARY_TSV"
        else
            FAIL=$((FAIL + 1))
            printf "%s\t%s\t%s\t%s\n" "$group" "$sample" "FAIL" "" >> "$SUMMARY_TSV"
            echo "  FAILED: $sample"
        fi
    done
}

echo ""
echo "Lifting clinical VCFs..."
liftover_group "clinical"

echo ""
echo "Lifting control VCFs..."
liftover_group "control"

echo ""
echo "=========================================="
echo "Done."
echo "Total samples: $TOTAL"
echo "Success:       $OK"
echo "Failed:        $FAIL"
echo "Summary:       $SUMMARY_TSV"
echo "=========================================="
