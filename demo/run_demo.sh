#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

READ1="${SCRIPT_DIR}/pkd1_sim_R1.fq"
READ2="${SCRIPT_DIR}/pkd1_sim_R2.fq"
REFERENCE="${SCRIPT_DIR}/ref.fa"
DEFAULT_OUTPUT_DIR="${SCRIPT_DIR}/output"
OUTPUT_DIR="${OUTPUT_DIR:-${DEFAULT_OUTPUT_DIR}}"
PREFIX="pkd1_demo"
if [[ -n "${CONDA_PREFIX:-}" && -x "${CONDA_PREFIX}/bin/python" ]]; then
    PYTHON_BIN="${CONDA_PREFIX}/bin/python"
else
    PYTHON_BIN="$(command -v python || true)"
fi

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

[[ -n "$PYTHON_BIN" ]] || fail "python not found. Activate the paradism environment first."
check_command mafft
check_command bowtie2
check_command bowtie2-build
check_command samtools

"$PYTHON_BIN" - <<'PY' || fail "Required Python packages missing. Activate or update the environment with: conda env update -f environment.yml --prune"
import importlib.util
import sys

missing = [module for module in ("Bio", "pysam", "rich") if importlib.util.find_spec(module) is None]
if missing:
    print(f"Python executable: {sys.executable}", file=sys.stderr)
    print("Missing Python modules: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
PY

check_file "$REFERENCE"

if [[ "${REGENERATE_READS:-0}" == "1" ]]; then
    bash "${SCRIPT_DIR}/generate_demo_reads.sh"
fi

check_file "$READ1"
check_file "$READ2"

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

"$PYTHON_BIN" paradism.py \
    --read1 "$READ1" \
    --read2 "$READ2" \
    --reference "$REFERENCE" \
    --aligner bowtie2 \
    --threads 1 \
    --iterations 1 \
    --output-dir "$OUTPUT_DIR" \
    --prefix "$PREFIX"

[[ -d "$OUTPUT_DIR" ]] || fail "Output directory was not created: $OUTPUT_DIR"
INITIAL_SAM="${OUTPUT_DIR}/iteration_1/mapped_reads.sam"
[[ -f "$INITIAL_SAM" ]] || fail "Expected initial SAM not found: $INITIAL_SAM"

FINAL_FASTQ_DIR="${OUTPUT_DIR}/final_outputs/${PREFIX}_fastq"
FINAL_BAM_DIR="${OUTPUT_DIR}/final_outputs/${PREFIX}_bam"
FINAL_NONE_DIR="${OUTPUT_DIR}/final_outputs/${PREFIX}_none"
UNASSIGNED_READ1="${FINAL_NONE_DIR}/${PREFIX}_NONE_r1.fq"
SUMMARY_PATH="${OUTPUT_DIR}/demo_assignment_summary.tsv"
[[ -d "$FINAL_FASTQ_DIR" ]] || fail "Final FASTQ directory not found: $FINAL_FASTQ_DIR"
[[ -d "$FINAL_BAM_DIR" ]] || fail "Final BAM directory not found: $FINAL_BAM_DIR"
check_file "$UNASSIGNED_READ1"

shopt -s nullglob
final_fastqs=("${FINAL_FASTQ_DIR}"/*.fq)
sorted_bams=("${FINAL_BAM_DIR}"/*.sorted.bam)
bam_indexes=("${FINAL_BAM_DIR}"/*.sorted.bam.bai)
shopt -u nullglob

check_nonempty_glob "final FASTQ files in ${FINAL_FASTQ_DIR}" "${final_fastqs[@]}"
check_nonempty_glob "sorted BAM files in ${FINAL_BAM_DIR}" "${sorted_bams[@]}"
check_nonempty_glob "BAM index files in ${FINAL_BAM_DIR}" "${bam_indexes[@]}"

"$PYTHON_BIN" "${SCRIPT_DIR}/summarize_demo.py" \
    --read1 "$READ1" \
    --assigned-dir "$FINAL_FASTQ_DIR" \
    --unassigned-read1 "$UNASSIGNED_READ1" \
    --prefix "$PREFIX" \
    --summary "$SUMMARY_PATH" \
    || fail "Demo assignments did not match the documented expected result."

echo "PKD1 demo completed successfully."
echo "Output directory: ${OUTPUT_DIR}"
echo "Initial SAM: ${INITIAL_SAM}"
echo "Final FASTQs: ${FINAL_FASTQ_DIR}"
echo "Final BAMs: ${FINAL_BAM_DIR}"
echo "Assignment summary: ${SUMMARY_PATH}"
