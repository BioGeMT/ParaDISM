#!/usr/bin/env bash
# Call RAW variants for HTS outputs.
# Produces per-gene raw VCFs and merged raw VCFs for final ParaDISM outputs.
#
# Example:
#   bash hts_analysis/call_variants_raw_hts.sh \
#     --run-root hts_analysis/HTS_rerun_ref1_fixed_all \
#     --reference /homes/dtzim01/ref.fa

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash hts_analysis/call_variants_raw_hts.sh [options]

Options:
  --run-root PATH           HTS run root (default: hts_analysis/HTS_rerun_ref1_fixed_all)
  --reference PATH          Reference FASTA used by ParaDISM (default: ref.fa in repo root)
  --aligner NAME            Output aligner tag (default: bowtie2)
  -h, --help                Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

RUN_ROOT="hts_analysis/HTS_rerun_ref1_fixed_all"
REFERENCE="${PROJECT_ROOT}/ref.fa"
ALIGNER="bowtie2"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root)
            RUN_ROOT="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
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

if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: Reference file not found: $REFERENCE"
    exit 1
fi

if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate paradism_env 2>/dev/null || true
fi

for tool in freebayes bcftools; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: required tool not found in PATH: $tool"
        exit 1
    fi
done

mapfile -t GENES < <(grep '^>' "$REFERENCE" | sed 's/^>//')
if [[ "${#GENES[@]}" -eq 0 ]]; then
    echo "Error: No FASTA headers found in $REFERENCE"
    exit 1
fi

FREEBAYES_ARGS=(--ploidy 2)

call_raw_from_gene_bams() {
    local sample_dir="$1"
    local sample="$2"

    local bam_dir="${sample_dir}/final_outputs/${sample}_bam"
    local out_dir="${sample_dir}/final_outputs/variant_calling_raw"
    local gene
    local gene_bam
    local gene_vcf_raw
    local gene_vcfgz_raw
    local -a gene_vcfs_raw=()

    if [[ ! -d "$bam_dir" ]]; then
        echo "    [final] Missing BAM directory: $bam_dir"
        return 1
    fi

    mkdir -p "${out_dir}/per_gene"

    for gene in "${GENES[@]}"; do
        gene_bam="${bam_dir}/${sample}_${gene}.sorted.bam"
        gene_vcf_raw="${out_dir}/per_gene/${gene}.raw.vcf"
        gene_vcfgz_raw="${gene_vcf_raw}.gz"

        if [[ -f "$gene_bam" ]]; then
            freebayes --bam "$gene_bam" --fasta-reference "$REFERENCE" \
                "${FREEBAYES_ARGS[@]}" > "$gene_vcf_raw"
            bcftools sort -Oz -o "$gene_vcfgz_raw" "$gene_vcf_raw"
            bcftools index -f "$gene_vcfgz_raw"
            gene_vcfs_raw+=("$gene_vcfgz_raw")
        fi
    done

    if [[ ${#gene_vcfs_raw[@]} -eq 0 ]]; then
        echo "    [final] No gene BAMs found for $sample"
        return 1
    fi

    bcftools concat -a -O v -o "${out_dir}/variants_raw.vcf" "${gene_vcfs_raw[@]}"
    bcftools sort -Oz -o "${out_dir}/variants_raw.vcf.gz" "${out_dir}/variants_raw.vcf"
    bcftools index -f "${out_dir}/variants_raw.vcf.gz"
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

        echo "[${group_label}] $sample"
        final_ok=0

        if call_raw_from_gene_bams "$sample_dir" "$sample"; then
            final_ok=1
        fi

        printf "%s\t%s\t%s\n" "$group_label" "$sample" "$final_ok" >> "$SUMMARY_FILE"
        TOTAL_SAMPLES=$((TOTAL_SAMPLES + 1))
        FINAL_SUCCESS=$((FINAL_SUCCESS + final_ok))
    done < <(find "$group_root" -mindepth 1 -maxdepth 1 -type d | sort)
}

CLINICAL_ROOT="${RUN_ROOT}/HTS_output/HTS_${ALIGNER}_output"
CONTROL_ROOT="${RUN_ROOT}/HTS_CNTRL/HTS_${ALIGNER}_output"
SUMMARY_FILE="${RUN_ROOT}/variant_calling_raw_summary_${ALIGNER}.tsv"

TOTAL_SAMPLES=0
FINAL_SUCCESS=0

echo -e "group\tsample\tfinal_raw_ok" > "$SUMMARY_FILE"

echo "=========================================="
echo "HTS RAW Variant Calling"
echo "Run root: $RUN_ROOT"
echo "Reference: $REFERENCE"
echo "Aligner: $ALIGNER"
echo "FreeBayes args: ${FREEBAYES_ARGS[*]}"
echo "Summary: $SUMMARY_FILE"
echo "=========================================="

process_group "clinical" "$CLINICAL_ROOT"
process_group "control" "$CONTROL_ROOT"

echo ""
echo "Done."
echo "Samples processed: $TOTAL_SAMPLES"
echo "Final raw success: $FINAL_SUCCESS"
echo "Summary TSV: $SUMMARY_FILE"
