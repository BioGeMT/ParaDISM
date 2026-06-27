#!/usr/bin/env bash
set -euo pipefail

# Get script directory (benchmark/simulation/) and ParaDISM root.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARADISM_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

usage() {
    cat <<'EOF'
Usage:
  bash benchmark/simulation/run_simulation.sh [options]

Runs one simulated-read benchmark for one reference/gene group. The script
generates reads with dwgsim, runs ParaDISM with the requested aligner(s), and
summarizes read-assignment metrics.

Options:
  --group NAME             Group label used in output paths (default: pkd1_panel)
  --reference PATH         Reference FASTA (default: benchmark/references/<group>.fa)
  --out PATH               Output directory (default: results/simulation/<group>)
  --seeds N                Run seeds 1..N (default: 1)
  --seed-start N           First seed (default: 1)
  --seed-end N             Last seed (default: 1)
  --reads-per-seed N       Read pairs per seed (default: 1000)
  --aligners LIST          Comma/space separated aligners (default: bwa-mem2,bowtie2,minimap2)
  --threads N              Threads per ParaDISM run (default: 2)
  --iterations N           ParaDISM iterations (default: 1)
  -h, --help               Show this help

Example:
  bash benchmark/simulation/run_simulation.sh \
    --group hba_pair \
    --reference benchmark/references/hba_pair.fa \
    --seeds 2 \
    --reads-per-seed 1000 \
    --out results/simulation/hba_pair
EOF
}

# ------------------------------------------------------------------
# Configuration (can override via env, e.g. SEED_END=10 ITERATIONS=2 ALIGNERS="bwa-mem2")
# ------------------------------------------------------------------
GROUP="${GROUP:-pkd1_panel}"
SEED_START=${SEED_START:-1}
SEED_END=${SEED_END:-1}
REFERENCE="${REFERENCE:-}"
SIM_OUTPUT_BASE="${SIM_OUTPUT_BASE:-}"
THREADS=${THREADS:-2}
NUM_READS=${NUM_READS:-1000}
ERROR_RATE=${ERROR_RATE:-0.01}             # default sequencing error rate (per base)
QUAL_THRESHOLD=${QUAL_THRESHOLD:-20}
DP_THRESHOLD=${DP_THRESHOLD:-10}
AF_THRESHOLD=${AF_THRESHOLD:-0.05}
MINIMAP2_PROFILE="${MINIMAP2_PROFILE:-short}"     # sr preset
SNP_RATE=${SNP_RATE:-0.005}              # DWGSIM SNP rate (per base)
INDEL_RATE=${INDEL_RATE:-0.0005}           # DWGSIM indel rate (per base)
INDEL_EXT=${INDEL_EXT:-0.5}               # DWGSIM indel extension probability
READ_LEN=${READ_LEN:-150}                # Read length
FRAG_MEAN=${FRAG_MEAN:-350}               # Mean fragment length
FRAG_SD=${FRAG_SD:-35}                  # Fragment length std dev
DWGSIM_DIR="${DWGSIM_DIR:-${SCRIPT_DIR}/dwgsim}"  # DWGSIM directory

ALIGNERS_STR="${ALIGNERS:-bwa-mem2,bowtie2,minimap2}"  # Space- or comma-separated
ITERATIONS=${ITERATIONS:-1}                # Number of ParaDISM runs (1 = no refinement)
SEEDS_PER_BATCH=${SEEDS_PER_BATCH:-1}
reference_was_set=0
output_was_set=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --group)
            GROUP="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            reference_was_set=1
            shift 2
            ;;
        --out)
            SIM_OUTPUT_BASE="$2"
            output_was_set=1
            shift 2
            ;;
        --seeds)
            SEED_START=1
            SEED_END="$2"
            shift 2
            ;;
        --seed-start)
            SEED_START="$2"
            shift 2
            ;;
        --seed-end)
            SEED_END="$2"
            shift 2
            ;;
        --reads-per-seed)
            NUM_READS="$2"
            shift 2
            ;;
        --aligners)
            ALIGNERS_STR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --iterations)
            ITERATIONS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ "$reference_was_set" -eq 0 ]]; then
    REFERENCE="${PARADISM_ROOT}/benchmark/references/${GROUP}.fa"
fi
if [[ "$output_was_set" -eq 0 ]]; then
    SIM_OUTPUT_BASE="${PARADISM_ROOT}/results/simulation/${GROUP}"
fi
IFS=', ' read -r -a ALIGNERS <<< "$ALIGNERS_STR"

if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: reference FASTA not found: $REFERENCE" >&2
    exit 1
fi

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
format_error_suffix() {
    local rate="$1"
    local formatted
    formatted=$(printf "%.3f" "$rate")
    echo "${formatted/./_}"
}

run_mapper() {
    local aligner="$1"
    local r1="$2"
    local r2="$3"
    local output_dir="$4"
    local prefix="$5"

    local cmd=(
        python "${PARADISM_ROOT}/paradism.py"
        --read1 "$r1"
        --read2 "$r2"
        --reference "$REFERENCE"
        --aligner "$aligner"
        --threads "$THREADS"
        --output-dir "$output_dir"
        --prefix "$prefix"
        --iterations "$ITERATIONS"
        --add-quality-filters
        --qual-threshold "$QUAL_THRESHOLD"
        --dp-threshold "$DP_THRESHOLD"
        --af-threshold "$AF_THRESHOLD"
    )
    if [[ "$aligner" == "minimap2" ]]; then
        cmd+=(--minimap2-profile "$MINIMAP2_PROFILE")
    fi

    # Run mapper (may exit with status 1 even if successful)
    "${cmd[@]}" || true

    # Verify final outputs exist
    local final_outputs_dir="${output_dir}/final_outputs"
    local final_fastq_dir="${final_outputs_dir}/${prefix}_fastq"

    if [[ ! -d "$final_outputs_dir" ]]; then
        echo "ERROR: Final outputs directory not found: $final_outputs_dir" >&2
        exit 1
    fi
    if [[ ! -d "$final_fastq_dir" ]]; then
        echo "ERROR: Final ParaDISM FASTQ directory not found: $final_fastq_dir" >&2
        exit 1
    fi
    echo "  Using final outputs FASTQ dir: $final_fastq_dir" >&2
}

write_timing_csv() {
    local output_csv="$1"
    local tmp_csv="${output_csv}.tmp"

    mkdir -p "$(dirname "$output_csv")"
    echo "seed,aligner,wall_clock_seconds,wall_clock_formatted" > "$tmp_csv"

    shopt -s nullglob
    local timing_files=("$TIMING_TMP_DIR"/*.timing)
    if [[ ${#timing_files[@]} -gt 0 ]]; then
        for timing_file in "${timing_files[@]}"; do
            cat "$timing_file" >> "$tmp_csv"
        done
        sort -t',' -k1,1n -k2,2 "$tmp_csv" > "${tmp_csv}.sorted"
        mv "${tmp_csv}.sorted" "$tmp_csv"
    fi
    shopt -u nullglob

    mv "$tmp_csv" "$output_csv"
}

run_postprocessing() {
    local seed_end="$1"
    local out_base="$2"
    local label="$3"

    echo "=============================="
    echo "${label}: Aggregating timing data (through seed ${seed_end})"
    echo "=============================="

    local timing_csv="${out_base}/timing_data.csv"
    write_timing_csv "$timing_csv"
    echo "Timing data written to: $timing_csv"

    echo "=============================="
    echo "${label}: Aggregating results (through seed ${seed_end})"
    echo "=============================="

    python "${SCRIPT_DIR}/aggregate_results.py" \
        --sim-output-base "$SIM_OUTPUT_BASE" \
        --seed-start "$SEED_START" \
        --seed-end "$seed_end" \
        --aligners "${ALIGNERS[@]}" \
        --output-dir "${out_base}/aggregated_results"

    python "${SCRIPT_DIR}/aggregate_per_seed_overall.py" \
        --sim-output-base "$SIM_OUTPUT_BASE" \
        --seed-start "$SEED_START" \
        --seed-end "$seed_end" \
        --aligners "${ALIGNERS[@]}" \
        --output "${out_base}/aggregated_results/per_seed_overall_metrics.csv"

    echo "${label}: Aggregation complete!"
}

# ------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------
ERROR_SUFFIX=$(format_error_suffix "$ERROR_RATE")

if (( SEED_START > SEED_END )); then
    echo "Error: SEED_START (${SEED_START}) must be <= SEED_END (${SEED_END})" >&2
    exit 1
fi

# Timing CSV setup
TIMING_CSV="${SIM_OUTPUT_BASE}/timing_data.csv"
TIMING_TMP_DIR="${SIM_OUTPUT_BASE}/timing_tmp"
mkdir -p "$TIMING_TMP_DIR"

# Process seeds in batches in parallel.
# Effective CPU use is approximately: SEEDS_PER_BATCH × number_of_aligners × THREADS.
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
    seed_pids=()
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
                echo "Error: dwgsim not found. Install the pinned environment with: mamba env create -f environment.yml" >&2
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
            aligner_pids=()
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

                    # Reuse ParaDISM's initial SAM for direct alignment comparison.
                    paradism_sam="$paradism_output/iteration_1/mapped_reads.sam"
                    if [[ ! -f "$paradism_sam" ]]; then
                        echo "Error: ParaDISM SAM file not found: $paradism_sam" >&2
                        exit 1
                    fi
                    base_bam="$aligner_dir/${aligner}_base.sorted.bam"
                    
                    # Convert ParaDISM's SAM to sorted BAM for analysis
                    samtools view -b "$paradism_sam" | samtools sort -o "$base_bam"
                    samtools index "$base_bam"

                    # Run iteration-by-iteration read mapping analysis.
                    analysis_dir="$aligner_dir/read_mapping_analysis"
                    mkdir -p "$analysis_dir"

                    # Sorted list of iteration directories.
                    while IFS= read -r iter_dir; do
                        iter_num="${iter_dir##*/iteration_}"
                        iter_fastq_dir="${iter_dir}/${paradism_prefix}_fastq"
                        if [[ -d "$iter_fastq_dir" ]]; then
                            python "${SCRIPT_DIR}/read_mapping_analysis.py" \
                                --aligner "$aligner" \
                                --fastq-r1 "$r1" \
                                --mapper-fastq-dir "$iter_fastq_dir" \
                                --direct-sam "$base_bam" \
                                --analysis-dir "$analysis_dir" \
                                --output-prefix "seed_${seed}_${aligner}_iter${iter_num}"
                        fi
                    done < <(find "$paradism_output" -maxdepth 1 -type d -name 'iteration_*' | sort -t_ -k2,2n)

                    # Run final summary analysis from final_outputs only (no fallback)
                    mapper_fastq_dir="${paradism_output}/final_outputs/${paradism_prefix}_fastq"
                    if [[ ! -d "$mapper_fastq_dir" ]]; then
                        echo "ERROR: Final ParaDISM FASTQ directory not found: $mapper_fastq_dir" >&2
                        exit 1
                    fi
                    echo "  Using final_outputs for aggregate summary analysis" >&2
                    python "${SCRIPT_DIR}/read_mapping_analysis.py" \
                        --aligner "$aligner" \
                        --fastq-r1 "$r1" \
                        --mapper-fastq-dir "$mapper_fastq_dir" \
                        --direct-sam "$base_bam" \
                        --analysis-dir "$analysis_dir" \
                        --output-prefix "seed_${seed}_${aligner}"

                    echo "--- Completed aligner: $aligner (seed $seed) ---"
                    echo "  ParaDISM time: ${formatted_time} (${paradism_duration} seconds)"
                ) &
                aligner_pids+=($!)
            done
            
            # Wait for all aligners to complete for this seed and fail if any failed
            aligner_failed=0
            for pid in "${aligner_pids[@]}"; do
                if ! wait "$pid"; then
                    aligner_failed=1
                fi
            done
            if [[ "$aligner_failed" -ne 0 ]]; then
                echo "ERROR: One or more aligner jobs failed for seed $seed" >&2
                exit 1
            fi
            echo "=============================="
            echo "Completed seed: $seed"
            echo "=============================="
        ) &
        seed_pids+=($!)
    done
    
    # Wait for all seeds in this batch to complete and fail if any failed
    batch_failed=0
    for pid in "${seed_pids[@]}"; do
        if ! wait "$pid"; then
            batch_failed=1
        fi
    done
    if [[ "$batch_failed" -ne 0 ]]; then
        echo "ERROR: One or more seed jobs failed in batch ${seed_array[$batch_start]} to ${seed_array[$((batch_end-1))]}" >&2
        exit 1
    fi

    echo "=============================="
    echo "Completed batch: ${seed_array[$batch_start]} to ${seed_array[$((batch_end-1))]}"
    echo "=============================="
done

# ------------------------------------------------------------------
# Final post-processing on full requested seed range
# ------------------------------------------------------------------
final_seed="${seed_array[$((total_seeds - 1))]}"
run_postprocessing "$final_seed" "$SIM_OUTPUT_BASE" "Final"
rm -rf "$TIMING_TMP_DIR"
echo "Timing temp files removed: $TIMING_TMP_DIR"
