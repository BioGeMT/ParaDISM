#!/bin/bash
# Run ParaDISM on GNAQ samples with bwa-mem2 threshold ~G,40,40 equivalent
# Run from anywhere: bash gnaq_analysis/run_gnaq_samples_bwa_160.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$SCRIPT_DIR"

# Configuration
DATA_DIR="$SCRIPT_DIR/GNAQ_reads"
REFERENCE="$SCRIPT_DIR/gnaq-gnaqp_ref.fa"
ALIGNER="bwa-mem2"
THREADS=12
ITERATIONS=10
MIN_ALT_COUNT=5
# For ~100bp reads, use score threshold 160 (ParaDISM default for 100bp)
THRESHOLD=180
OUTPUT_BASE="$SCRIPT_DIR/bwa-mem2_180_output"
LOG_DIR="$OUTPUT_BASE/mapper_logs"

if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE not found!"
    exit 1
fi

# Activate conda environment
if command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate paradism_env 2>/dev/null || true
fi

mkdir -p "$LOG_DIR"

SAMPLES=(
    "SRR5602384"
    "SRR5602389"
    "SRR5602393"
    "SRR5602414"
    "SRR5602419"
)

echo "=========================================="
echo "Processing ${#SAMPLES[@]} GNAQ samples (bwa-mem2, threshold=${THRESHOLD})"
echo "=========================================="

for sample in "${SAMPLES[@]}"; do
    (
        r1_file="$DATA_DIR/${sample}_1.fastq"
        r2_file="$DATA_DIR/${sample}_2.fastq"
        output_dir="$OUTPUT_BASE/$sample"

        if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
            echo "Error: Files not found for $sample"
            exit 1
        fi

        echo "Processing: $sample"
        log_file="$LOG_DIR/${sample}.log"

        mapper_cmd=(
            python "$PROJECT_ROOT/paradism.py"
            --read1 "$r1_file"
            --read2 "$r2_file"
            --reference "$REFERENCE"
            --aligner "$ALIGNER"
            --threads "$THREADS"
            --output-dir "$output_dir"
            --iterations "$ITERATIONS"
            --min-alternate-count "$MIN_ALT_COUNT"
            --threshold "$THRESHOLD"
            --add-quality-filters
        )

        if "${mapper_cmd[@]}" >> "$log_file" 2>&1; then
            echo "Completed: $sample"
        else
            echo "Failed: $sample"
        fi
    ) &
done

wait
echo "All GNAQ samples processed!"
