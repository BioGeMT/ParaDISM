#!/bin/bash

# Get script directory (simulation/) and ParaDISM root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARADISM_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Activate conda environment if available
if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
    source "${HOME}/miniconda3/etc/profile.d/conda.sh"
    conda activate paradism_env 2>/dev/null || {
        echo "Warning: Could not activate paradism_env. Make sure it exists: conda env create -f ${PARADISM_ROOT}/paradism_env.yml"
    }
fi

set -euo pipefail

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
SEED_START=1
SEED_END=1             # Number of seeds (for error bars)
SIM_OUTPUT_BASE="${SCRIPT_DIR}/sim_output_iterations"
REFERENCE="${PARADISM_ROOT}/ref.fa"
THREADS=1                   # Threads per ParaDISM run (10 seeds × 1 aligner × 1 thread = 10 CPUs)
NUM_READS=1000
ERROR_RATE=0.01             # default sequencing error rate (per base)
MINIMAP2_PROFILE="short"     # sr preset
SNP_RATE=0.005              # DWGSIM SNP rate (per base)
INDEL_RATE=0.0005           # DWGSIM indel rate (per base)
INDEL_EXT=0.5               # DWGSIM indel extension probability
READ_LEN=150                # Read length
FRAG_MEAN=350               # Mean fragment length
FRAG_SD=35                  # Fragment length std dev
DWGSIM_DIR="${PARADISM_ROOT}/../dwgsim"  # DWGSIM directory

ALIGNERS=("bwa-mem2")
ITERATIONS=5               # Number of ParaDISM runs (10 = run 10 times, 9 refinement iterations)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
format_error_suffix() {
    local rate="$1"
    local formatted
    formatted=$(printf "%.3f" "$rate")
    echo "${formatted/./_}"
}

find_latest_iteration_dir() {
    local base="$1"
    local latest=""
    local latest_num=0

    shopt -s nullglob
    for path in "${base}"/iteration_*; do
        [[ -d "$path" ]] || continue
        local suffix="${path##*/iteration_}"
        if [[ "$suffix" =~ ^([0-9]+)$ ]]; then
            local num="${BASH_REMATCH[1]}"
            if (( num > latest_num )); then
                latest_num=$num
                latest="$path"
            fi
        fi
    done
    shopt -u nullglob

    echo "$latest"
}

run_mapper() {
    local aligner="$1"
    local r1="$2"
    local r2="$3"
    local output_dir="$4"
    local prefix="$5"

    local cmd=(python "${PARADISM_ROOT}/paradism.py" --read1 "$r1" --read2 "$r2" --reference "$REFERENCE" --aligner "$aligner" --threads "$THREADS" --output-dir "$output_dir" --prefix "$prefix" --iterations "$ITERATIONS")
    if [[ "$aligner" == "minimap2" ]]; then
        cmd+=(--minimap2-profile "$MINIMAP2_PROFILE")
    fi

    # Run mapper (may exit with status 1 even if successful)
    "${cmd[@]}" || true
    
    # Verify output file was created; allow early convergence (final iteration may be < ITERATIONS)
    local final_iter_dir
    final_iter_dir=$(find_latest_iteration_dir "$output_dir")
    if [[ -z "$final_iter_dir" ]]; then
        echo "ERROR: No iteration directories found under $output_dir" >&2
        exit 1
    fi

    local expected_output="${final_iter_dir}/${prefix}_unique_mappings.tsv"
    if [[ ! -f "$expected_output" ]]; then
        echo "ERROR: Expected output file not found after convergence: $expected_output" >&2
        exit 1
    fi
}

# ------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------
ERROR_SUFFIX=$(format_error_suffix "$ERROR_RATE")

# Process seeds in batches of 100 in parallel (100 seeds × 1 aligner × 1 thread = 100 CPUs)
SEEDS_PER_BATCH=100
seed_array=($(seq "$SEED_START" "$SEED_END"))
total_seeds=${#seed_array[@]}

for ((batch_start=0; batch_start<total_seeds; batch_start+=SEEDS_PER_BATCH)); do
    batch_end=$((batch_start + SEEDS_PER_BATCH))
    if [[ $batch_end -gt $total_seeds ]]; then
        batch_end=$total_seeds
    fi
    
    echo "=============================="
    echo "Processing seeds batch: ${seed_array[$batch_start]} to ${seed_array[$((batch_end-1))]}"
    echo "=============================="
    
    # Process all seeds in this batch in parallel
    for ((i=batch_start; i<batch_end; i++)); do
        seed=${seed_array[$i]}
        (
            # Activate conda environment in subshell for background process
            if [ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
                source "${HOME}/miniconda3/etc/profile.d/conda.sh"
                conda activate paradism_env 2>/dev/null
            fi

            echo "=============================="
            echo "Processing seed: $seed"
            echo "=============================="

    seed_dir="$SIM_OUTPUT_BASE/seed_${seed}"
    mkdir -p "$seed_dir"

    echo "Simulating reads for seed $seed using DWGSIM..."
    
    # Set up DWGSIM binary path
    DWGSIM_BIN=""
    if command -v dwgsim >/dev/null 2>&1; then
        DWGSIM_BIN="$(command -v dwgsim)"
    elif [ -x "${DWGSIM_DIR}/dwgsim" ]; then
        DWGSIM_BIN="${DWGSIM_DIR}/dwgsim"
    else
        echo "Error: DWGSIM not found. Please run run_dwgsim_simulation.sh first or install DWGSIM." >&2
        exit 1
    fi
    
    # Run DWGSIM simulation
    prefix_seed="${seed_dir}/simulated_seed${seed}"
    "${DWGSIM_BIN}" \
        -z "${seed}" \
        -N "${NUM_READS}" \
        -1 "${READ_LEN}" -2 "${READ_LEN}" \
        -d "${FRAG_MEAN}" -s "${FRAG_SD}" \
        -y 0 \
        -e "${ERROR_RATE}" -E "${ERROR_RATE}" \
        -r "${SNP_RATE}" -R "${INDEL_RATE}" -X "${INDEL_EXT}" \
        "${REFERENCE}" "${prefix_seed}"
    
    # DWGSIM outputs .bwa.read1.fastq.gz and .bwa.read2.fastq.gz
    # Convert to uncompressed and rename to match expected format
    dwgsim_r1="${prefix_seed}.bwa.read1.fastq.gz"
    dwgsim_r2="${prefix_seed}.bwa.read2.fastq.gz"
    r1="${seed_dir}/simulated_r1_err_${ERROR_SUFFIX}.fq"
    r2="${seed_dir}/simulated_r2_err_${ERROR_SUFFIX}.fq"
    
    if [[ -f "$dwgsim_r1" && -f "$dwgsim_r2" ]]; then
        # Decompress and rename
        gunzip -c "$dwgsim_r1" > "$r1"
        gunzip -c "$dwgsim_r2" > "$r2"
        # Delete compressed files to save disk space
        rm -f "$dwgsim_r1" "$dwgsim_r2"
    else
        echo "Error: DWGSIM output files not found for seed $seed" >&2
        exit 1
    fi
    
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "Missing simulated FASTQs for seed $seed" >&2
        exit 1
    fi

    # Process aligner
    for aligner in "${ALIGNERS[@]}"; do
        echo "--- Running aligner: $aligner (seed $seed) ---"
        aligner_dir="$seed_dir/$aligner"
        paradism_output="$aligner_dir/paradism_output"
        prefix="seed_${seed}_${aligner}"
        paradism_prefix="paradism_${prefix}"
        mkdir -p "$paradism_output"

        run_mapper "$aligner" "$r1" "$r2" "$paradism_output" "$paradism_prefix"

        # Reuse ParaDISM's SAM file from iteration 1 for direct alignment comparison
        # (This is the base alignment before iterative refinement)
        paradism_sam="$paradism_output/iteration_1/mapped_reads.sam"
        base_bam="$aligner_dir/${aligner}_base.sorted.bam"
        
        if [[ ! -f "$paradism_sam" ]]; then
            echo "Error: ParaDISM SAM file not found: $paradism_sam" >&2
            exit 1
        fi
        
        # Convert ParaDISM's SAM to sorted BAM for analysis
        samtools view -b "$paradism_sam" | samtools sort -o "$base_bam"
        samtools index "$base_bam"

        # Analyze each iteration (1-based: iter1, iter2, ..., iterN)
        analysis_dir="$aligner_dir/read_mapping_analysis"
        mkdir -p "$analysis_dir"
        
        # Analyze all iterations (1 to ITERATIONS)
        for iter in $(seq 1 $ITERATIONS); do
            mapper_tsv="$paradism_output/iteration_${iter}/${paradism_prefix}_unique_mappings.tsv"
            if [[ -f "$mapper_tsv" ]]; then
                output_prefix="seed_${seed}_${aligner}_iter${iter}"
                
                python "${SCRIPT_DIR}/read_mapping_analysis.py" \
                    --aligner "$aligner" \
                    --fastq-r1 "$r1" \
                    --mapper-tsv "$mapper_tsv" \
                    --direct-sam "$base_bam" \
                    --analysis-dir "$analysis_dir" \
                    --output-prefix "$output_prefix"
            fi
        done
        
        echo "--- Completed aligner: $aligner (seed $seed) ---"
    done
    
            echo "=============================="
            echo "Completed seed: $seed"
            echo "=============================="
        ) &
    done
    
    # Wait for all seeds in this batch to complete
    wait
    echo "=============================="
    echo "Completed batch: ${seed_array[$batch_start]} to ${seed_array[$((batch_end-1))]}"
    echo "=============================="
done

# ------------------------------------------------------------------
# Create iteration progression plots
# ------------------------------------------------------------------
echo "=============================="
echo "Creating iteration progression plots..."
echo "=============================="

python "${SCRIPT_DIR}/plot_iteration_progression.py" \
    --sim-output-base "$SIM_OUTPUT_BASE" \
    --seed-start "$SEED_START" \
    --seed-end "$SEED_END" \
    --output-dir "${SIM_OUTPUT_BASE}/iteration_plots" \
    --iterations "$ITERATIONS" \
    --aligner "${ALIGNERS[0]}"

echo "=============================="
echo "Simulation complete!"
echo "Results saved to: ${SIM_OUTPUT_BASE}/"
echo "Iteration plots saved to: ${SIM_OUTPUT_BASE}/iteration_plots/"
echo "=============================="
