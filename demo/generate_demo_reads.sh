#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFERENCE="${SCRIPT_DIR}/ref.fa"
WORK_DIR="${WORK_DIR:-${SCRIPT_DIR}/generated_reads}"
READ1_OUT="${READ1_OUT:-${SCRIPT_DIR}/reads_R1.fq}"
READ2_OUT="${READ2_OUT:-${SCRIPT_DIR}/reads_R2.fq}"
PREFIX="${WORK_DIR}/demo_seed1"
READ_PAIRS=20

command -v dwgsim >/dev/null 2>&1 || {
    echo "ERROR: dwgsim not found on PATH. Activate the paradism environment first." >&2
    exit 1
}

mkdir -p "$WORK_DIR"
dwgsim \
    -z 1 \
    -N "$READ_PAIRS" \
    -1 150 -2 150 \
    -d 350 -s 35 \
    -y 0 \
    -e 0.01 -E 0.01 \
    -r 0.001 -R 0.0001 -X 0.5 \
    "$REFERENCE" "$PREFIX"

gzip -cd "${PREFIX}.bwa.read1.fastq.gz" > "$READ1_OUT"
gzip -cd "${PREFIX}.bwa.read2.fastq.gz" > "$READ2_OUT"

echo "Generated ${READ_PAIRS} paired-end demo read pairs:"
echo "  Reference: ${REFERENCE}"
echo "  ${READ1_OUT}"
echo "  ${READ2_OUT}"
