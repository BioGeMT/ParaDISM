#!/usr/bin/env bash
# Build IGV-ready per-sample outputs from HTS pipeline results:
# - one merged BAM per sample (from per-gene final BAMs)
# - one uncompressed final VCF per sample
#
# Example:
#   bash hts_analysis/prepare_igv_per_sample_hts.sh \
#     --run-root hts_analysis/HTS_rerun_bowtie2_t4_p10 \
#     --aligner bowtie2

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash hts_analysis/prepare_igv_per_sample_hts.sh [options]

Options:
  --run-root PATH     HTS run root (default: hts_analysis/HTS_rerun_bowtie2_t4_p10)
  --aligner NAME      Output aligner tag (default: bowtie2)
  --out-dir PATH      Output directory for IGV-ready files
                      (default: <run-root>/igv_per_sample)
  --threads N         Threads for samtools merge/index (default: 8)
  -h, --help          Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

RUN_ROOT="hts_analysis/HTS_rerun_bowtie2_t4_p10"
ALIGNER="bowtie2"
OUT_DIR=""
THREADS=8

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
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
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

if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR="${RUN_ROOT}/igv_per_sample"
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS <= 0 )); then
    echo "Error: --threads must be a positive integer (got: $THREADS)"
    exit 1
fi

SAMTOOLS="$(command -v samtools || true)"
BCFTOOLS="$(command -v bcftools || true)"

if [[ -z "$SAMTOOLS" && -x "$HOME/miniconda3/envs/paradism_env/bin/samtools" ]]; then
    SAMTOOLS="$HOME/miniconda3/envs/paradism_env/bin/samtools"
fi
if [[ -z "$BCFTOOLS" && -x "$HOME/miniconda3/envs/paradism_env/bin/bcftools" ]]; then
    BCFTOOLS="$HOME/miniconda3/envs/paradism_env/bin/bcftools"
fi

if [[ -z "$SAMTOOLS" ]]; then
    echo "Error: samtools not found in PATH and not found in paradism_env."
    exit 1
fi
if [[ -z "$BCFTOOLS" ]]; then
    echo "Error: bcftools not found in PATH and not found in paradism_env."
    exit 1
fi

CLINICAL_ROOT="${RUN_ROOT}/HTS_output/HTS_${ALIGNER}_output"
CONTROL_ROOT="${RUN_ROOT}/HTS_CNTRL/HTS_${ALIGNER}_output"

for root in "$CLINICAL_ROOT" "$CONTROL_ROOT"; do
    if [[ ! -d "$root" ]]; then
        echo "Error: expected directory not found: $root"
        exit 1
    fi
done

mkdir -p "$OUT_DIR/clinical" "$OUT_DIR/control"

SUMMARY_TSV="${OUT_DIR}/igv_export_summary.tsv"
echo -e "group\tsample\tbam_ok\tvcf_ok\tbam_path\tvcf_path" > "$SUMMARY_TSV"

TOTAL=0
BAM_OK=0
VCF_OK=0

process_group() {
    local group_label="$1"
    local group_root="$2"
    local group_out_dir="$OUT_DIR/$group_label"

    local sample_dir sample bam_dir final_vcfgz out_bam out_vcf
    local bam_ok vcf_ok
    local -a gene_bams

    while IFS= read -r sample_dir; do
        sample="$(basename "$sample_dir")"
        if [[ "$sample" == "mapper_logs" ]]; then
            continue
        fi
        if [[ ! -d "${sample_dir}/final_outputs" ]]; then
            continue
        fi

        TOTAL=$((TOTAL + 1))
        bam_ok=0
        vcf_ok=0

        bam_dir="${sample_dir}/final_outputs/${sample}_bam"
        final_vcfgz="${sample_dir}/final_outputs/variant_calling_balanced/${sample}_variants_final.vcf.gz"

        out_bam="${group_out_dir}/${sample}.bam"
        out_vcf="${group_out_dir}/${sample}.vcf"

        mkdir -p "$group_out_dir"

        if [[ -d "$bam_dir" ]]; then
            mapfile -t gene_bams < <(find "$bam_dir" -maxdepth 1 -type f -name "${sample}_*.sorted.bam" | sort)
            if [[ "${#gene_bams[@]}" -gt 0 ]]; then
                if "$SAMTOOLS" merge -f -@ "$THREADS" -o "$out_bam" "${gene_bams[@]}"; then
                    rm -f "${out_bam}.bai"
                    if "$SAMTOOLS" index -@ "$THREADS" "$out_bam"; then
                        bam_ok=1
                        BAM_OK=$((BAM_OK + 1))
                    fi
                fi
            fi
        fi

        if [[ -f "$final_vcfgz" ]]; then
            if "$BCFTOOLS" view -Ov -o "$out_vcf" "$final_vcfgz"; then
                vcf_ok=1
                VCF_OK=$((VCF_OK + 1))
            fi
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$group_label" "$sample" "$bam_ok" "$vcf_ok" "$out_bam" "$out_vcf" >> "$SUMMARY_TSV"
    done < <(find "$group_root" -mindepth 1 -maxdepth 1 -type d | sort)
}

echo "=========================================="
echo "Preparing IGV per-sample outputs"
echo "Run root: $RUN_ROOT"
echo "Aligner: $ALIGNER"
echo "Out dir: $OUT_DIR"
echo "Threads: $THREADS"
echo "=========================================="

process_group "clinical" "$CLINICAL_ROOT"
process_group "control" "$CONTROL_ROOT"

echo ""
echo "Done."
echo "Samples processed: $TOTAL"
echo "Merged BAM success: $BAM_OK"
echo "Final VCF success: $VCF_OK"
echo "Summary: $SUMMARY_TSV"
