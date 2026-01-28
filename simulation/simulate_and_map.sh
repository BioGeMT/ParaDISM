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
# Configuration (can override via env, e.g. SEED_END=10 ITERATIONS=10 ALIGNERS="bwa-mem2")
# ------------------------------------------------------------------
SEED_START=${SEED_START:-1}
SEED_END=${SEED_END:-1000}
SIM_OUTPUT_BASE="${SIM_OUTPUT_BASE:-${SCRIPT_DIR}/simulations_outputs}"
REFERENCE="${REFERENCE:-${PARADISM_ROOT}/ref.fa}"
THREADS=${THREADS:-1}                   # Threads per ParaDISM run (30 seeds × 3 aligners × 1 thread = 90 CPUs)
NUM_READS=${NUM_READS:-100000}
ERROR_RATE=${ERROR_RATE:-0.01}             # default sequencing error rate (per base)
MINIMAP2_PROFILE="${MINIMAP2_PROFILE:-short}"     # sr preset
SNP_RATE=${SNP_RATE:-0.005}              # DWGSIM SNP rate (per base)
INDEL_RATE=${INDEL_RATE:-0.0005}           # DWGSIM indel rate (per base)
INDEL_EXT=${INDEL_EXT:-0.5}               # DWGSIM indel extension probability
READ_LEN=${READ_LEN:-150}                # Read length
FRAG_MEAN=${FRAG_MEAN:-350}               # Mean fragment length
FRAG_SD=${FRAG_SD:-35}                  # Fragment length std dev
DWGSIM_DIR="${DWGSIM_DIR:-${SCRIPT_DIR}/dwgsim}"  # DWGSIM directory

ALIGNERS_STR="${ALIGNERS:-bwa-mem2 bowtie2 minimap2}"  # Space- or comma-separated
IFS=', ' read -r -a ALIGNERS <<< "$ALIGNERS_STR"
ITERATIONS=${ITERATIONS:-10}                # Number of ParaDISM runs (2 = run twice, 1 refinement iteration)

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

    # Verify final outputs exist (prefer final_outputs, fall back to last iteration for early convergence)
    local final_outputs_dir="${output_dir}/final_outputs"
    local final_tsv="${final_outputs_dir}/${prefix}_unique_mappings.tsv"

    if [[ ! -f "$final_tsv" ]]; then
        echo "ERROR: Final outputs TSV not found: $final_tsv" >&2
        exit 1
    fi
    echo "  Using final outputs TSV: $final_tsv" >&2
}

# ------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------
ERROR_SUFFIX=$(format_error_suffix "$ERROR_RATE")

# Timing CSV setup
TIMING_CSV="${SIM_OUTPUT_BASE}/timing_data.csv"
TIMING_TMP_DIR="${SIM_OUTPUT_BASE}/timing_tmp"
mkdir -p "$TIMING_TMP_DIR"

# Process seeds in batches of 30 in parallel (30 seeds × 3 aligners × 1 thread = 90 CPUs)
SEEDS_PER_BATCH=30
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
                echo "Error: DWGSIM not found. Please run simulation/run_dwgsim_simulation.sh first or install DWGSIM." >&2
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

            # Process all aligners in parallel for this seed
            for aligner in "${ALIGNERS[@]}"; do
                (
                    echo "--- Running aligner: $aligner (seed $seed) ---"
                    aligner_dir="$seed_dir/$aligner"
                    paradism_output="$aligner_dir/paradism_output"
                    prefix="seed_${seed}_${aligner}"
                    paradism_prefix="paradism_${prefix}"
                    mkdir -p "$paradism_output"

                    # Start timing for ParaDISM mapper only
                    paradism_start_time=$(date +%s)
                    paradism_start_date=$(date '+%Y-%m-%d %H:%M:%S')

                    run_mapper "$aligner" "$r1" "$r2" "$paradism_output" "$paradism_prefix"

                    # End timing for ParaDISM mapper
                    paradism_end_time=$(date +%s)
                    paradism_end_date=$(date '+%Y-%m-%d %H:%M:%S')
                    paradism_duration=$((paradism_end_time - paradism_start_time))
                    hours=$((paradism_duration / 3600))
                    minutes=$(((paradism_duration % 3600) / 60))
                    seconds=$((paradism_duration % 60))
                    formatted_time=$(printf "%02d:%02d:%02d" $hours $minutes $seconds)
                    
                    # Write timing to temporary file (for parallel safety)
                    timing_tmp_file="$TIMING_TMP_DIR/seed_${seed}_${aligner}.timing"
                    echo "$seed,$aligner,$paradism_duration,$formatted_time" > "$timing_tmp_file"

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

                    # Use final iteration's mappings for ParaDISM comparison
                    # Find the actual last iteration that was produced (handles early convergence)
                    final_iter_dir=$(find_latest_iteration_dir "$paradism_output")
                    if [[ -z "$final_iter_dir" ]]; then
                        echo "ERROR: No iteration directories found under $paradism_output" >&2
                        exit 1
                    fi

                    # Run iteration-by-iteration read mapping analysis (for progression plots)
                    analysis_dir="$aligner_dir/read_mapping_analysis"
                    mkdir -p "$analysis_dir"

                    # Sorted list of iteration directories
                    mapfile -t iter_dirs < <(find "$paradism_output" -maxdepth 1 -type d -name 'iteration_*' | sort -t_ -k2,2n)
                    for iter_dir in "${iter_dirs[@]}"; do
                        iter_num="${iter_dir##*/iteration_}"
                        iter_tsv="${iter_dir}/${paradism_prefix}_unique_mappings.tsv"
                        if [[ -f "$iter_tsv" ]]; then
                            python "${SCRIPT_DIR}/read_mapping_analysis.py" \
                                --aligner "$aligner" \
                                --fastq-r1 "$r1" \
                                --mapper-tsv "$iter_tsv" \
                                --direct-sam "$base_bam" \
                                --analysis-dir "$analysis_dir" \
                                --output-prefix "${analysis_dir}/seed_${seed}_${aligner}_iter${iter_num}"
                        fi
                    done

                    # Also run analysis on the final iteration for aggregate summary
                    mapper_tsv="${final_iter_dir}/${paradism_prefix}_unique_mappings.tsv"
                    if [[ ! -f "$mapper_tsv" ]]; then
                        echo "ERROR: Expected output file not found: $mapper_tsv" >&2
                        exit 1
                    fi
                    final_iter_num="${final_iter_dir##*/iteration_}"
                    echo "  Using final iteration ${final_iter_num} for analysis" >&2
                    python "${SCRIPT_DIR}/read_mapping_analysis.py" \
                        --aligner "$aligner" \
                        --fastq-r1 "$r1" \
                        --mapper-tsv "$mapper_tsv" \
                        --direct-sam "$base_bam" \
                        --analysis-dir "$analysis_dir" \
                        --output-prefix "seed_${seed}_${aligner}"

                    echo "--- Completed aligner: $aligner (seed $seed) ---"
                    echo "  ParaDISM time: ${formatted_time} (${paradism_duration} seconds)"
                ) &
            done
            
            # Wait for all aligners to complete for this seed
            wait
            echo "=============================="
            echo "Completed seed: $seed"
            echo "=============================="
        ) &
    done
    
    # Wait for all seeds in this batch to complete before moving to next batch
    wait
    echo "=============================="
    echo "Completed batch: ${seed_array[$batch_start]} to ${seed_array[$((batch_end-1))]}"
    echo "=============================="
done

# ------------------------------------------------------------------
# Aggregate timing data
# ------------------------------------------------------------------
echo "=============================="
echo "Aggregating timing data..."
echo "=============================="

# Create CSV with header and combine all timing files
echo "seed,aligner,wall_clock_seconds,wall_clock_formatted" > "$TIMING_CSV"

# Append all timing files
for timing_file in "$TIMING_TMP_DIR"/*.timing; do
    if [[ -f "$timing_file" ]]; then
        cat "$timing_file" >> "$TIMING_CSV"
    fi
done

# Sort CSV by seed, then aligner
sort -t',' -k1,1n -k2,2 "$TIMING_CSV" > "${TIMING_CSV}.tmp"
mv "${TIMING_CSV}.tmp" "$TIMING_CSV"

# Clean up temporary timing files
rm -rf "$TIMING_TMP_DIR"

echo "Timing data aggregated to: $TIMING_CSV"

# ------------------------------------------------------------------
# Aggregate results across all seeds
# ------------------------------------------------------------------
echo "=============================="
echo "Aggregating results across all seeds..."
echo "=============================="

python "${SCRIPT_DIR}/aggregate_results.py" \
    --sim-output-base "$SIM_OUTPUT_BASE" \
    --seed-start "$SEED_START" \
    --seed-end "$SEED_END" \
    --aligners "${ALIGNERS[@]}" \
    --output-dir "${SIM_OUTPUT_BASE}/aggregated_results"

# Per-seed overall (ParaDISM vs direct aligner)
python "${SCRIPT_DIR}/aggregate_per_seed_overall.py" \
    --sim-output-base "$SIM_OUTPUT_BASE" \
    --seed-start "$SEED_START" \
    --seed-end "$SEED_END" \
    --aligners "${ALIGNERS[@]}" \
    --output "${SIM_OUTPUT_BASE}/aggregated_results/per_seed_overall_metrics.csv"

echo "Aggregation complete!"

# ------------------------------------------------------------------
# Plot iteration progression per aligner (auto-detects iteration counts per seed)
# ------------------------------------------------------------------
echo "=============================="
echo "Creating iteration progression plots..."
echo "=============================="

PLOT_DIR="${SIM_OUTPUT_BASE}/iteration_plots"
mkdir -p "$PLOT_DIR"

for aligner in "${ALIGNERS[@]}"; do
    python "${SCRIPT_DIR}/plot_iteration_progression.py" \
        --sim-output-base "$SIM_OUTPUT_BASE" \
        --seed-start "$SEED_START" \
        --seed-end "$SEED_END" \
        --output-dir "$PLOT_DIR" \
        --aligner "$aligner"
done

echo "Iteration plots saved to: $PLOT_DIR"
