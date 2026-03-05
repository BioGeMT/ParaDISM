#!/usr/bin/env bash
# Filter HTS final RAW VCFs using the same final GIAB expression:
#   INFO/AO>=3 && QUAL>=10 && INFO/SAF>=1 && INFO/SAR>=1
#
# Example:
#   bash hts_analysis/filter_variants_hts.sh \
#     --run-root hts_analysis/HTS_rerun_ref1_fixed_all

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash hts_analysis/filter_variants_hts.sh [options]

Options:
  --run-root PATH      HTS run root (default: hts_analysis/HTS_rerun_ref1_fixed_all)
  --aligner NAME       Output aligner tag (default: bowtie2)
  -h, --help           Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

RUN_ROOT="hts_analysis/HTS_rerun_ref1_fixed_all"
ALIGNER="bowtie2"
FILTER_EXPR='INFO/AO>=3 && QUAL>=10 && INFO/SAF>=1 && INFO/SAR>=1'

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root)
            RUN_ROOT="$2"
            shift 2
            ;;
        --aligner)
            ALIGNER="$2"
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

BCFTOOLS="$(command -v bcftools || true)"
if [[ -z "$BCFTOOLS" && -x "$HOME/miniconda3/envs/paradism_env/bin/bcftools" ]]; then
    BCFTOOLS="$HOME/miniconda3/envs/paradism_env/bin/bcftools"
fi
if [[ -z "$BCFTOOLS" ]]; then
    echo "Error: bcftools not found in PATH and not found in paradism_env."
    exit 1
fi

filter_vcf() {
    local sample="$1"
    local raw_vcf="$2"
    local out_dir="$3"

    if [[ ! -f "$raw_vcf" ]]; then
        return 1
    fi

    mkdir -p "$out_dir"
    local out_vcfgz="${out_dir}/${sample}_variants_final.vcf.gz"

    "$BCFTOOLS" view -i "$FILTER_EXPR" -Oz -o "$out_vcfgz" "$raw_vcf"
    "$BCFTOOLS" index -f "$out_vcfgz"
    return 0
}

process_group() {
    local group_label="$1"
    local group_root="$2"

    if [[ ! -d "$group_root" ]]; then
        echo "[${group_label}] Directory not found, skipping: $group_root"
        return
    fi

    local sample_dir
    local sample
    local final_ok

    while IFS= read -r sample_dir; do
        sample="$(basename "$sample_dir")"
        if [[ "$sample" == "mapper_logs" ]]; then
            continue
        fi
        if [[ ! -d "${sample_dir}/final_outputs" ]]; then
            continue
        fi

        final_ok=0

        if filter_vcf \
            "$sample" \
            "${sample_dir}/final_outputs/variant_calling_raw/variants_raw.vcf.gz" \
            "${sample_dir}/final_outputs/variant_calling_balanced"; then
            final_ok=1
        fi

        printf "%s\t%s\t%s\n" "$group_label" "$sample" "$final_ok" >> "$SUMMARY_FILE"
        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
        FINAL_FILTER_SUCCESS=$((FINAL_FILTER_SUCCESS + final_ok))
    done < <(find "$group_root" -mindepth 1 -maxdepth 1 -type d | sort)
}

CLINICAL_ROOT="${RUN_ROOT}/HTS_output/HTS_${ALIGNER}_output"
CONTROL_ROOT="${RUN_ROOT}/HTS_CNTRL/HTS_${ALIGNER}_output"
SUMMARY_FILE="${RUN_ROOT}/variant_filtering_summary_${ALIGNER}.tsv"

TOTAL_SAMPLES=0
FINAL_FILTER_SUCCESS=0

echo -e "group\tsample\tfinal_filter_ok" > "$SUMMARY_FILE"

echo "=========================================="
echo "HTS Variant Filtering (GIAB-style final filter)"
echo "Run root: $RUN_ROOT"
echo "Aligner: $ALIGNER"
echo "Filter expression: $FILTER_EXPR"
echo "Summary: $SUMMARY_FILE"
echo "=========================================="

process_group "clinical" "$CLINICAL_ROOT"
process_group "control" "$CONTROL_ROOT"

echo ""
echo "Done."
echo "Samples processed: $TOTAL_SAMPLES"
echo "Final filter success: $FINAL_FILTER_SUCCESS"
echo "Summary TSV: $SUMMARY_FILE"
