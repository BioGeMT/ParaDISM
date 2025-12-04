#!/bin/bash

set -euo pipefail

# Configuration
DATA_DIR="/mnt/STORAGE-BioGeMT-01/pkd_data"
REFERENCE="ref.fa"
ALIGNER="bwa-mem2"
THREADS=4
OUTPUT_BASE="HTS_bwa_output"
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
    "HTS009"
    "HTS131"
    "HTS132"
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
    "IndexCHKPEPI00000968_HTS002"
)

# Process samples in batches to use 90 CPUs (18 samples × 5 threads = 90 CPUs)
SAMPLES_PER_BATCH=18
total_samples=${#SAMPLES[@]}

for ((batch_start=0; batch_start<total_samples; batch_start+=SAMPLES_PER_BATCH)); do
    batch_end=$((batch_start + SAMPLES_PER_BATCH))
    if [[ $batch_end -gt $total_samples ]]; then
        batch_end=$total_samples
    fi

    echo "=========================================="
    echo "Processing samples batch: ${batch_start} to $((batch_end - 1))"
    echo "=========================================="

    # Process all samples in this batch in parallel
    for ((i=batch_start; i<batch_end; i++)); do
        sample="${SAMPLES[$i]}"
        (
            # Handle different naming patterns
            if [[ "$sample" == "IndexCHKPEPI00000968_HTS002" ]]; then
                r1_file="$DATA_DIR/${sample}_1.fq"
                r2_file="$DATA_DIR/${sample}_2.fq"
            else
                r1_file="$DATA_DIR/${sample}_R1.fq"
                r2_file="$DATA_DIR/${sample}_R2.fq"

                # For HTS131 and HTS132, merge L1/L2 lanes once and reuse
                if [[ "$sample" == "HTS131" || "$sample" == "HTS132" ]]; then
                    base="$sample"
                    r1_merged="$OUTPUT_BASE/${base}_merged_R1.fq"
                    r2_merged="$OUTPUT_BASE/${base}_merged_R2.fq"
                    if [[ ! -f "$r1_merged" || ! -f "$r2_merged" ]]; then
                        cat "$DATA_DIR/${base}_L1_R1.fq" "$DATA_DIR/${base}_L2_R1.fq" > "$r1_merged"
                        cat "$DATA_DIR/${base}_L1_R2.fq" "$DATA_DIR/${base}_L2_R2.fq" > "$r2_merged"
                    fi
                    r1_file="$r1_merged"
                    r2_file="$r2_merged"
                fi
            fi

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
                --iterations 10
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

    # Wait for all samples in this batch to complete
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
