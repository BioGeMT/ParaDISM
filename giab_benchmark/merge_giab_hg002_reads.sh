#!/bin/bash
# Merge GIAB HG002 lane chunks into single R1/R2 files.
# Output: giab_benchmark/giab_hg002_reads/HG002_R1.fq.gz and HG002_R2.fq.gz

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
READS_DIR="$SCRIPT_DIR/giab_hg002_reads"
cd "$READS_DIR"

shopt -s nullglob

R1_FILES=(D1_S1_L00*_R1_*.fastq.gz)
R2_FILES=(D1_S1_L00*_R2_*.fastq.gz)

if [[ ${#R1_FILES[@]} -eq 0 || ${#R2_FILES[@]} -eq 0 ]]; then
    echo "Error: No lane chunk FASTQs found in $READS_DIR"
    echo "Expected files like D1_S1_L001_R1_001.fastq.gz"
    exit 1
fi

echo "Merging ${#R1_FILES[@]} R1 chunks..."
printf '%s\n' "${R1_FILES[@]}" | LC_ALL=C sort | xargs cat > HG002_R1.fq.gz

echo "Merging ${#R2_FILES[@]} R2 chunks..."
printf '%s\n' "${R2_FILES[@]}" | LC_ALL=C sort | xargs cat > HG002_R2.fq.gz

echo "Done: HG002_R1.fq.gz and HG002_R2.fq.gz"
