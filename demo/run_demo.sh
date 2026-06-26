#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

READ1="${SCRIPT_DIR}/reads_R1.fq"
READ2="${SCRIPT_DIR}/reads_R2.fq"
REFERENCE="${SCRIPT_DIR}/ref.fa"
DEFAULT_OUTPUT_DIR="${SCRIPT_DIR}/output"
OUTPUT_DIR="${OUTPUT_DIR:-${DEFAULT_OUTPUT_DIR}}"
PREFIX="tiny_demo"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

check_command() {
    local command_name="$1"
    command -v "$command_name" >/dev/null 2>&1 || fail "$command_name not found on PATH. Activate the paradism environment first."
}

check_file() {
    local path="$1"
    [[ -f "$path" ]] || fail "Required input not found: $path"
}

check_nonempty_glob() {
    local label="$1"
    shift
    local matches=("$@")
    [[ ${#matches[@]} -gt 0 ]] || fail "No ${label} found"
}

check_command python
check_command mafft
check_command bowtie2
check_command bowtie2-build
check_command samtools

python - <<'PY' || fail "Required Python packages missing. Activate the paradism environment first."
import Bio
import pysam
import rich
PY

check_file "$READ1"
check_file "$READ2"
check_file "$REFERENCE"

case "$OUTPUT_DIR" in
    ""|"/"|"/tmp"|"/var/tmp")
        fail "Refusing unsafe OUTPUT_DIR: ${OUTPUT_DIR}"
        ;;
esac

cd "$PROJECT_ROOT"
if [[ -e "$OUTPUT_DIR" ]]; then
    if [[ "$OUTPUT_DIR" == "$DEFAULT_OUTPUT_DIR" ]]; then
        rm -rf "$OUTPUT_DIR"
    else
        fail "OUTPUT_DIR already exists; remove it first or choose a new path: $OUTPUT_DIR"
    fi
fi

python paradism.py \
    --read1 "$READ1" \
    --read2 "$READ2" \
    --reference "$REFERENCE" \
    --aligner bowtie2 \
    --threads 1 \
    --iterations 1 \
    --output-dir "$OUTPUT_DIR" \
    --prefix "$PREFIX"

[[ -d "$OUTPUT_DIR" ]] || fail "Output directory was not created: $OUTPUT_DIR"
[[ -f "${OUTPUT_DIR}/mapped_reads.sam" ]] || fail "Expected one-iteration SAM not found: ${OUTPUT_DIR}/mapped_reads.sam"

FINAL_FASTQ_DIR="${OUTPUT_DIR}/final_outputs/${PREFIX}_fastq"
FINAL_BAM_DIR="${OUTPUT_DIR}/final_outputs/${PREFIX}_bam"
[[ -d "$FINAL_FASTQ_DIR" ]] || fail "Final FASTQ directory not found: $FINAL_FASTQ_DIR"
[[ -d "$FINAL_BAM_DIR" ]] || fail "Final BAM directory not found: $FINAL_BAM_DIR"

shopt -s nullglob
final_fastqs=("${FINAL_FASTQ_DIR}"/*.fq)
sorted_bams=("${FINAL_BAM_DIR}"/*.sorted.bam)
bam_indexes=("${FINAL_BAM_DIR}"/*.sorted.bam.bai)
shopt -u nullglob

check_nonempty_glob "final FASTQ files in ${FINAL_FASTQ_DIR}" "${final_fastqs[@]}"
check_nonempty_glob "sorted BAM files in ${FINAL_BAM_DIR}" "${sorted_bams[@]}"
check_nonempty_glob "BAM index files in ${FINAL_BAM_DIR}" "${bam_indexes[@]}"

echo "Tiny demo completed successfully."
echo "Output directory: ${OUTPUT_DIR}"
echo "SAM: ${OUTPUT_DIR}/mapped_reads.sam"
echo "Final FASTQs: ${FINAL_FASTQ_DIR}"
echo "Final BAMs: ${FINAL_BAM_DIR}"
