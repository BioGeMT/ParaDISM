#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFERENCE="${SCRIPT_DIR}/ref.fa"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/generated_reads}"
PREFIX="${OUT_DIR}/demo_seed1"

command -v dwgsim >/dev/null 2>&1 || {
    echo "ERROR: dwgsim not found on PATH. Activate the paradism environment first." >&2
    exit 1
}

mkdir -p "$OUT_DIR"
dwgsim \
    -z 1 \
    -N 20 \
    -1 150 -2 150 \
    -d 350 -s 35 \
    -y 0 \
    -e 0.01 -E 0.01 \
    -r 0.001 -R 0.0001 -X 0.5 \
    "$REFERENCE" "$PREFIX"

gzip -cd "${PREFIX}.bwa.read1.fastq.gz" > "${OUT_DIR}/reads_R1.fq"
gzip -cd "${PREFIX}.bwa.read2.fastq.gz" > "${OUT_DIR}/reads_R2.fq"

echo "Generated demo reads:"
echo "  ${OUT_DIR}/reads_R1.fq"
echo "  ${OUT_DIR}/reads_R2.fq"
