#!/bin/bash

set -euo pipefail

# Configuration
DATA_DIR="/mnt/STORAGE-BioGeMT-01/pkd_data"
REFERENCE="ref.fa"
ALIGNER="bwa-mem2"
THREADS=64
OUTPUT_BASE="output"
LOG_DIR="mapper_logs"

# Check if reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE not found!"
    exit 1
fi

# Create log directory
mkdir -p "$LOG_DIR"
echo "Logs will be saved to: $LOG_DIR/"

SAMPLES=(
    "HTS009"
    "HTS131_L1"
    "HTS131_L2"
    "HTS132_L1"
    "HTS132_L2"
    "HTS158"
    "HTS159"
    "HTS160"
    "HTS161"
    "HTS162"
    "HTS163"
    "HTS164"
    "HTS165"
    "HTS166"
    "HTS167"
    "HTS168"
    "HTS169"
    "HTS170"
    "HTS171"
)

# Process each sample
for sample in "${SAMPLES[@]}"; do
    r1_file="$DATA_DIR/${sample}_R1.fq"
    r2_file="$DATA_DIR/${sample}_R2.fq"
    output_dir="$OUTPUT_BASE/$sample"

    # Check if files exist
    if [ ! -f "$r1_file" ]; then
        echo "Error: R1 file not found: $r1_file"
        continue
    fi
    if [ ! -f "$r2_file" ]; then
        echo "Error: R2 file not found: $r2_file"
        continue
    fi

    echo "=========================================="
    echo "Processing sample: $sample"
    echo "R1: $r1_file"
    echo "R2: $r2_file"
    echo "Output: $output_dir"
    echo "=========================================="

    # Set log file path
    log_file="$LOG_DIR/${sample}.log"

    # Run mapper.py with time command, save all output to log
    if /usr/bin/time -v python mapper.py \
        --read1 "$r1_file" \
        --read2 "$r2_file" \
        --reference "$REFERENCE" \
        --aligner "$ALIGNER" \
        --threads "$THREADS" \
        --output-dir "$output_dir" > "$log_file" 2>&1; then
        status="SUCCESS"
    else
        status="FAILED"
    fi

    echo "Status: $status" >> "$log_file"

    echo "Completed: $sample (Status: $status, Log: $log_file)"
    echo ""
done

echo "=========================================="
echo "All samples processed!"
echo "Logs saved to: $LOG_DIR/"
echo "=========================================="

