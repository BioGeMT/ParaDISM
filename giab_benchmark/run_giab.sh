#!/bin/bash
# Run ParaDISM on GIAB HG002 with G60 threshold
# Uses bowtie2, minalt 5, qual filtered
# Run from anywhere: bash giab_benchmark/run_giab.sh

set -euo pipefail

# Get script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$SCRIPT_DIR"

READS_DIR="$SCRIPT_DIR/giab_hg002_reads"
REFERENCE="$SCRIPT_DIR/ref.fa"
PARADISM="$PROJECT_ROOT/paradism.py"
THREADS="${1:-8}"
ITERATIONS="${2:-10}"
MIN_ALT_COUNT=5

# Find merged reads
if [[ -f "${READS_DIR}/HG002_R1.fq.gz" ]]; then
    R1_MERGED="${READS_DIR}/HG002_R1.fq.gz"
    R2_MERGED="${READS_DIR}/HG002_R2.fq.gz"
elif [[ -f "${READS_DIR}/HG002_R1.fq" ]]; then
    R1_MERGED="${READS_DIR}/HG002_R1.fq"
    R2_MERGED="${READS_DIR}/HG002_R2.fq"
else
    echo "Error: Merged reads not found in $READS_DIR"
    echo "Expected: HG002_R1.fq.gz and HG002_R2.fq.gz"
    exit 1
fi

echo "=========================================="
echo "ParaDISM GIAB HG002 - G60"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Read1: $R1_MERGED"
echo "  Read2: $R2_MERGED"
echo "  Reference: $REFERENCE"
echo "  Threads: $THREADS"
echo "  Iterations: $ITERATIONS"
echo "  Min-alternate-count: $MIN_ALT_COUNT"
echo ""

# Activate conda environment
if command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate paradism_env 2>/dev/null || true
fi

# Run G60 (recommended threshold)
OUTPUT_G60="$SCRIPT_DIR/giab_hg002_output_bowtie2_G60_min5_qfilters"
if [[ -d "$OUTPUT_G60/final_outputs" ]]; then
    echo "G60 already complete, skipping..."
else
    echo "=========================================="
    echo "Running G60 threshold..."
    echo "=========================================="
    python "$PARADISM" \
        --read1 "$R1_MERGED" \
        --read2 "$R2_MERGED" \
        --reference "$REFERENCE" \
        --aligner bowtie2 \
        --threads "$THREADS" \
        --iterations "$ITERATIONS" \
        --min-alternate-count "$MIN_ALT_COUNT" \
        --add-quality-filters \
        --qual-threshold 20 \
        --dp-threshold 10 \
        --af-threshold 0.05 \
        --threshold "G,60,60" \
        --output-dir "$OUTPUT_G60"
fi

echo ""
echo "=========================================="
echo "Complete!"
echo "=========================================="
echo ""
echo "Output directories:"
echo "  G60: $OUTPUT_G60"
echo ""
echo "Next: Run variant calling"
echo "  bash giab_benchmark/call_variants_raw_g60.sh"
echo "  bash giab_benchmark/filter_variants_g60.sh"
