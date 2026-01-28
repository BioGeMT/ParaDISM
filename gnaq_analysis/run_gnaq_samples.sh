#!/bin/bash
# Run ParaDISM on GNAQ samples
# Run from anywhere: bash gnaq_analysis/run_gnaq_samples.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# Configuration
DATA_DIR="gnaq_analysis/gnaq_reads"
REFERENCE="gnaq_analysis/gnaq-gnaqp_ref.fa"
ALIGNER="minimap2"
THREADS=2
OUTPUT_BASE="gnaq_analysis/GNAQ_minimap2_output"
LOG_DIR="$OUTPUT_BASE/mapper_logs"
ITERATIONS=10
MINIMAP2_PROFILE="short"

if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE not found!"
    exit 1
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
echo "Processing ${#SAMPLES[@]} GNAQ samples"
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
            python paradism.py
            --read1 "$r1_file"
            --read2 "$r2_file"
            --reference "$REFERENCE"
            --aligner "$ALIGNER"
            --threads "$THREADS"
            --output-dir "$output_dir"
            --iterations "$ITERATIONS"
            --minimap2-profile "$MINIMAP2_PROFILE"
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
