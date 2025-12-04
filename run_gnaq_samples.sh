#!/bin/bash

set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Configuration
DATA_DIR="$SCRIPT_DIR/gnaq_reads"
REFERENCE="gnaq-gnaqp_ref.fa"
ALIGNER="minimap2"
THREADS=2
OUTPUT_BASE="GNAQ_minimap2_output"
LOG_DIR="$OUTPUT_BASE/mapper_logs"
ITERATIONS=10
# Only used when ALIGNER=minimap2 (valid presets: short, pacbio-hifi, pacbio-clr, ont-q20, ont-standard)
MINIMAP2_PROFILE="short"

# Check if reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE not found!"
    exit 1
fi

# Create log directory inside the output hierarchy
mkdir -p "$LOG_DIR"
echo "Logs will be saved to: $LOG_DIR/"

# GNAQ samples (SRA accessions)
SAMPLES=(
    "SRR5602384"
    "SRR5602389"
    "SRR5602393"
    "SRR5602414"
    "SRR5602419"
)

echo "=========================================="
echo "Processing ${#SAMPLES[@]} GNAQ samples with $ITERATIONS iterations"
echo "=========================================="

# Process all samples in parallel
for sample in "${SAMPLES[@]}"; do
    (
        r1_file="$DATA_DIR/${sample}_1.fastq"
        r2_file="$DATA_DIR/${sample}_2.fastq"
        output_dir="$OUTPUT_BASE/$sample"

        # Check if files exist
        if [ ! -f "$r1_file" ]; then
            echo "Error: R1 file not found: $r1_file"
            exit 1
        fi
        if [ ! -f "$r2_file" ]; then
            echo "Error: R2 file not found: $r2_file"
            exit 1
        fi

        echo "=========================================="
        echo "Processing sample: $sample"
        echo "R1: $r1_file"
        echo "R2: $r2_file"
        echo "Output: $output_dir"
        echo "=========================================="

        # Set log file path
        log_file="$LOG_DIR/${sample}.log"

        # Capture start time
        start_time=$(date +%s)
        start_date=$(date '+%Y-%m-%d %H:%M:%S')
        echo "Start time: $start_date" | tee -a "$log_file"

        # Run paradism.py with time command, save all output to log
        mapper_cmd=(
            python paradism.py
            --read1 "$r1_file"
            --read2 "$r2_file"
            --reference "$REFERENCE"
            --aligner "$ALIGNER"
            --threads "$THREADS"
            --output-dir "$output_dir"
            --iterations "$ITERATIONS"
        )

        if [[ "$ALIGNER" == "minimap2" ]]; then
            mapper_cmd+=(--minimap2-profile "$MINIMAP2_PROFILE")
        fi

        if /usr/bin/time -v "${mapper_cmd[@]}" >> "$log_file" 2>&1; then
            status="SUCCESS"
        else
            status="FAILED"
        fi

        # Capture end time and calculate duration
        end_time=$(date +%s)
        end_date=$(date '+%Y-%m-%d %H:%M:%S')
        duration=$((end_time - start_time))
        hours=$((duration / 3600))
        minutes=$(((duration % 3600) / 60))
        seconds=$((duration % 60))

        echo "End time: $end_date" | tee -a "$log_file"
        echo "Status: $status" | tee -a "$log_file"
        printf "Wall clock time: %02d:%02d:%02d (%d seconds)\n" $hours $minutes $seconds $duration | tee -a "$log_file"

        echo "Completed: $sample (Status: $status, Time: ${hours}h ${minutes}m ${seconds}s, Log: $log_file)"
    ) &
done

# Wait for all samples to complete
wait

echo "=========================================="
echo "All GNAQ samples processed!"
echo "Logs saved to: $LOG_DIR/"
echo "=========================================="
