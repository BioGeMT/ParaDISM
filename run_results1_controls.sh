#!/bin/bash

# Run ParaDISM on a set of controls from /mnt/ngs-data/Results 1
# Uses *_1.fq / *_2.fq naming inside each sample directory.

set -euo pipefail

DATA_BASE="/mnt/ngs-data/Results 1"
REFERENCE="ref.fa"
ALIGNER="bwa-mem2"
THREADS=4
ITERATIONS=10
OUTPUT_BASE="HTS_results1_output"
LOG_DIR="$OUTPUT_BASE/mapper_logs"
MINIMAP2_PROFILE="short"  # only used if ALIGNER=minimap2

# Samples to run (exclude ones already processed elsewhere)
SAMPLES=(
    "HTS003"
    "HTS005"
    "HTS007"
    "HTS010"
    "HTS011"
    "HTS012"
    "HTS013"
    "HTS014"
    "HTS015"
    "HTS016"
    "HTS017"
    "HTS018"
    "HTS019"
    "HTS020"
    "HTS021"
    "HTS022"
)

# Check reference
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE not found."
    exit 1
fi

mkdir -p "$LOG_DIR"
echo "Logs will be saved to: $LOG_DIR/"

SAMPLES_PER_BATCH=16
total_samples=${#SAMPLES[@]}

for ((batch_start=0; batch_start<total_samples; batch_start+=SAMPLES_PER_BATCH)); do
    batch_end=$((batch_start + SAMPLES_PER_BATCH))
    if [[ $batch_end -gt $total_samples ]]; then
        batch_end=$total_samples
    fi

    echo "=========================================="
    echo "Processing samples batch: ${batch_start} to $((batch_end - 1))"
    echo "=========================================="

    for ((i=batch_start; i<batch_end; i++)); do
        sample="${SAMPLES[$i]}"
        (
            sample_dir="$DATA_BASE/$sample"
            if [[ ! -d "$sample_dir" ]]; then
                echo "Error: Sample directory not found: $sample_dir"
                exit 1
            fi

            # Find R1/R2 (first match of *_1.fq / *_2.fq)
            r1_file=$(ls "$sample_dir"/*_1.fq 2>/dev/null | head -n 1 || true)
            r2_file=$(ls "$sample_dir"/*_2.fq 2>/dev/null | head -n 1 || true)

            if [[ -z "$r1_file" || -z "$r2_file" ]]; then
                echo "Error: Missing FASTQs for $sample in $sample_dir"
                exit 1
            fi

            output_dir="$OUTPUT_BASE/$sample"
            mkdir -p "$output_dir"

            echo "=========================================="
            echo "Processing sample: $sample"
            echo "R1: $r1_file"
            echo "R2: $r2_file"
            echo "Output: $output_dir"
            echo "=========================================="

            log_file="$LOG_DIR/${sample}.log"
            start_time=$(date +%s)
            echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$log_file"

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

            end_time=$(date +%s)
            duration=$((end_time - start_time))
            hours=$((duration / 3600))
            minutes=$(((duration % 3600) / 60))
            seconds=$((duration % 60))
            echo "End time: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "$log_file"
            echo "Status: $status" | tee -a "$log_file"
            printf "Wall clock time: %02d:%02d:%02d (%d seconds)\n" $hours $minutes $seconds $duration | tee -a "$log_file"
            echo "Completed: $sample (Status: $status, Log: $log_file)"
        ) &
    done

    wait
    echo "=========================================="
    echo "Completed batch: ${batch_start} to $((batch_end - 1))"
    echo "=========================================="
    echo ""
done

echo "=========================================="
echo "All samples processed!"
echo "Logs saved to: $LOG_DIR/"
echo "=========================================="
