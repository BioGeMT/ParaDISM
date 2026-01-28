#!/bin/bash
# Run ParaDISM on HTS samples
# Run from anywhere: bash hts_analysis/run_hts_samples.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

# Configuration
DATA_DIR="/mnt/STORAGE-BioGeMT-01/pkd_data"
REFERENCE="ref.fa"
ALIGNER="bowtie2"
THREADS=8
OUTPUT_BASE="hts_analysis/HTS_bowtie2_output"
LOG_DIR="$OUTPUT_BASE/mapper_logs"
ITERATIONS=10

if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file $REFERENCE not found!"
    exit 1
fi

if [ ! -d "$DATA_DIR" ]; then
    echo "Error: Data directory $DATA_DIR not accessible!"
    exit 1
fi

mkdir -p "$LOG_DIR"

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
    "IndexCHKPEPI00000968_HTS002"
)

echo "=========================================="
echo "Processing ${#SAMPLES[@]} HTS samples"
echo "=========================================="

for sample in "${SAMPLES[@]}"; do
    # Handle different naming patterns
    if [[ "$sample" == "IndexCHKPEPI00000968_HTS002" ]]; then
        r1_file="$DATA_DIR/${sample}_1.fq"
        r2_file="$DATA_DIR/${sample}_2.fq"
    else
        r1_file="$DATA_DIR/${sample}_R1.fq"
        r2_file="$DATA_DIR/${sample}_R2.fq"
    fi

    output_dir="$OUTPUT_BASE/$sample"

    if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
        echo "Error: Files not found for $sample"
        continue
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
    )

    if "${mapper_cmd[@]}" >> "$log_file" 2>&1; then
        echo "Completed: $sample"
    else
        echo "Failed: $sample"
    fi
done

echo "All HTS samples processed!"
