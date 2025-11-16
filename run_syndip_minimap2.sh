#!/bin/bash

set -euo pipefail

# Configuration
DATA_DIR="/homes/dtzim01/ParaDISM/syndip_reads"
REFERENCE="ref.fa"
ALIGNER="minimap2"
THREADS=16
OUTPUT_BASE="syndip_minimap2_output"
LOG_DIR="$OUTPUT_BASE/mapper_logs"
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

SAMPLES=(
    "ERR1341793"
    "ERR1341794"
    "ERR1341795"
    "ERR1341796"
)

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

        # Run mapper.py with time command, save all output to log
        mapper_cmd=(
            python mapper.py
            --read1 "$r1_file"
            --read2 "$r2_file"
            --reference "$REFERENCE"
            --aligner "$ALIGNER"
            --threads "$THREADS"
            --output-dir "$output_dir"
        )

        if [[ "$ALIGNER" == "minimap2" ]]; then
            mapper_cmd+=(--minimap2-profile "$MINIMAP2_PROFILE")
        fi

        # Run mapper (may exit with status 1 even if successful)
        /usr/bin/time -v "${mapper_cmd[@]}" >> "$log_file" 2>&1 || true

        # Check if output files were created successfully (this is the real success indicator)
        expected_output="$output_dir/${sample}_unique_mappings.tsv"
        if [ -f "$expected_output" ]; then
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
        echo ""
        
        # Exit with 0 if successful, 1 if failed
        if [ "$status" = "SUCCESS" ]; then
            exit 0
        else
            exit 1
        fi
    ) &
done

# Wait for all samples to complete
wait

echo "=========================================="
echo "All samples processed!"
echo "Logs saved to: $LOG_DIR/"
echo "=========================================="

