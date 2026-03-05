#!/usr/bin/env bash
# Run ParaDISM for all clinical HTS samples and controls with separate output roots.
# Example:
#   bash hts_analysis/run_hts_and_controls.sh \
#     --reference /homes/dtzim01/ref.fa \
#     --control-data-dir /path/to/control_fastqs \
#     --out-root hts_analysis/HTS_rerun_ref1_fixed_all

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash hts_analysis/run_hts_and_controls.sh [options]

Options:
  --reference PATH            Reference FASTA (default: ref.fa)
  --aligner NAME              bowtie2|bwa-mem2|minimap2 (default: bowtie2)
  --threads N                 Threads for ParaDISM (default: 8)
  --iterations N              ParaDISM iterations (default: 10)
  --disable-quality-filters   Disable ParaDISM variant quality filters
                              (enabled by default)
  --qual-threshold N          Minimum QUAL for variant quality filters (default: 20)
  --dp-threshold N            Minimum DP for variant quality filters (default: 10)
  --af-threshold F            Minimum AF for variant quality filters (default: 0.05)
  --clinical-data-dir PATH    Directory containing clinical FASTQs
                              (default: /mnt/STORAGE-BioGeMT-01/pkd_data)
  --control-data-dir PATH     Directory containing control FASTQs (required if auto-detect fails)
  --out-root PATH             Output root
                              (default: hts_analysis/HTS_rerun_all)
  --search-depth N            Max depth for FASTQ discovery (default: 5)
  --parallel-samples N        Max samples to run in parallel (default: 1)
  --dry-run                   Resolve inputs and print commands without running ParaDISM
  -h, --help                  Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

REFERENCE="ref.fa"
ALIGNER="bowtie2"
THREADS=8
ITERATIONS=10
ADD_QUALITY_FILTERS=1
QUAL_THRESHOLD=20
DP_THRESHOLD=10
AF_THRESHOLD=0.05
CLINICAL_DATA_DIR="/mnt/STORAGE-BioGeMT-01/pkd_data"
CONTROL_DATA_DIR=""
OUT_ROOT="hts_analysis/HTS_rerun_all"
SEARCH_DEPTH=5
PARALLEL_SAMPLES=1
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --aligner)
            ALIGNER="$2"
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
        --disable-quality-filters)
            ADD_QUALITY_FILTERS=0
            shift
            ;;
        --qual-threshold)
            QUAL_THRESHOLD="$2"
            shift 2
            ;;
        --dp-threshold)
            DP_THRESHOLD="$2"
            shift 2
            ;;
        --af-threshold)
            AF_THRESHOLD="$2"
            shift 2
            ;;
        --clinical-data-dir)
            CLINICAL_DATA_DIR="$2"
            shift 2
            ;;
        --control-data-dir)
            CONTROL_DATA_DIR="$2"
            shift 2
            ;;
        --out-root)
            OUT_ROOT="$2"
            shift 2
            ;;
        --search-depth)
            SEARCH_DEPTH="$2"
            shift 2
            ;;
        --parallel-samples)
            PARALLEL_SAMPLES="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: Reference file not found: $REFERENCE"
    exit 1
fi

if [[ ! -d "$CLINICAL_DATA_DIR" ]]; then
    echo "Error: Clinical data directory not found: $CLINICAL_DATA_DIR"
    exit 1
fi

case "$ALIGNER" in
    bowtie2|bwa-mem2|minimap2)
        ;;
    *)
        echo "Error: Unsupported aligner '$ALIGNER' (expected: bowtie2|bwa-mem2|minimap2)"
        exit 1
        ;;
esac

if ! [[ "$PARALLEL_SAMPLES" =~ ^[0-9]+$ ]] || (( PARALLEL_SAMPLES <= 0 )); then
    echo "Error: --parallel-samples must be a positive integer (got: $PARALLEL_SAMPLES)"
    exit 1
fi

if [[ -z "$CONTROL_DATA_DIR" ]]; then
    for candidate in \
        "/mnt/STORAGE-BioGeMT-01/pkd_data_controls" \
        "/mnt/STORAGE-BioGeMT-01/pkd_data_control" \
        "/mnt/STORAGE-BioGeMT-01/pkd_data_cntrl" \
        "/mnt/STORAGE-BioGeMT-01/controls/pkd_data"; do
        if [[ -d "$candidate" ]]; then
            CONTROL_DATA_DIR="$candidate"
            break
        fi
    done
fi

if [[ -z "$CONTROL_DATA_DIR" || ! -d "$CONTROL_DATA_DIR" ]]; then
    echo "Error: Control data directory not found."
    echo "Pass it explicitly with: --control-data-dir /path/to/control_fastqs"
    exit 1
fi

CLINICAL_SAMPLES=(
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
    "HTS172"
    "IndexCHKPEPI00000968_HTS002"
)

declare -A CLINICAL_LANE_MAP=(
    ["HTS131"]="HTS131_L1 HTS131_L2"
    ["HTS132"]="HTS132_L1 HTS132_L2"
)

CONTROL_SAMPLES=(
    "HTS002"
    "HTS003"
    "HTS005"
    "HTS007"
    "HTS031"
    "HTS034"
    "HTS035"
    "HTS036"
    "HTS037"
    "HTS038"
    "HTS039"
    "HTS040"
    "HTS041"
    "HTS042"
    "HTS043"
    "HTS044"
    "HTS045"
    "HTS046"
    "HTS047"
    "HTS048"
    "HTS049"
    "HTS067"
    "HTS068"
    "HTS070"
)

find_read_pair() {
    local data_root="$1"
    local sample="$2"
    local depth="$3"

    local -a suffix_pairs=(
        "_R1.fq:_R2.fq"
        "_R1.fastq:_R2.fastq"
        "_R1.fq.gz:_R2.fq.gz"
        "_R1.fastq.gz:_R2.fastq.gz"
        "_1.fq:_2.fq"
        "_1.fastq:_2.fastq"
        "_1.fq.gz:_2.fq.gz"
        "_1.fastq.gz:_2.fastq.gz"
    )

    local pair r1_suffix r2_suffix r1 r2
    for pair in "${suffix_pairs[@]}"; do
        IFS=":" read -r r1_suffix r2_suffix <<< "$pair"
        while IFS= read -r r1; do
            r2="${r1%$r1_suffix}$r2_suffix"
            if [[ -f "$r2" ]]; then
                printf "%s\t%s\n" "$r1" "$r2"
                return 0
            fi
        done < <(find "$data_root" -maxdepth "$depth" -type f -name "${sample}${r1_suffix}" 2>/dev/null | sort)
    done

    return 1
}

append_fastq_content() {
    local src="$1"
    local dest="$2"

    if [[ "$src" == *.gz ]]; then
        gzip -dc "$src" >> "$dest"
    else
        cat "$src" >> "$dest"
    fi
}

resolve_pair_for_sample() {
    local group_label="$1"
    local data_root="$2"
    local sample="$3"
    local depth="$4"

    local lane_members lane lane_pair lane_r1 lane_r2
    local merged_dir merged_r1 merged_r2

    lane_members="${CLINICAL_LANE_MAP[$sample]-}"
    if [[ "$group_label" == "clinical" && -n "$lane_members" ]]; then
        merged_dir="$OUT_ROOT/_merged_fastq"
        merged_r1="$merged_dir/${sample}_R1.fq"
        merged_r2="$merged_dir/${sample}_R2.fq"

        if [[ "$DRY_RUN" -eq 0 ]]; then
            mkdir -p "$merged_dir"
            : > "$merged_r1"
            : > "$merged_r2"
        fi

        for lane in $lane_members; do
            if ! lane_pair="$(find_read_pair "$data_root" "$lane" "$depth")"; then
                return 1
            fi
            IFS=$'\t' read -r lane_r1 lane_r2 <<< "$lane_pair"
            if [[ "$DRY_RUN" -eq 0 ]]; then
                append_fastq_content "$lane_r1" "$merged_r1"
                append_fastq_content "$lane_r2" "$merged_r2"
            fi
        done

        printf "%s\t%s\n" "$merged_r1" "$merged_r2"
        return 0
    fi

    find_read_pair "$data_root" "$sample" "$depth"
}

CLINICAL_OUTPUT_BASE="$OUT_ROOT/HTS_output/HTS_${ALIGNER}_output"
CONTROL_OUTPUT_BASE="$OUT_ROOT/HTS_CNTRL/HTS_${ALIGNER}_output"
mkdir -p "$CLINICAL_OUTPUT_BASE/mapper_logs" "$CONTROL_OUTPUT_BASE/mapper_logs"

BATCH_LOG="$OUT_ROOT/batch_${ALIGNER}_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$OUT_ROOT"

TOTAL_PROCESSED=0
TOTAL_SUCCESS=0
TOTAL_FAILED=0
TOTAL_MISSING=0

process_sample() {
    local group_label="$1"
    local data_root="$2"
    local output_base="$3"
    local sample="$4"

    local pair r1_file r2_file sample_output sample_log
    local -a cmd

    echo "[$(date)] [${group_label}] Resolving reads for ${sample}" | tee -a "$BATCH_LOG"

    if [[ "$group_label" == "clinical" && -n "${CLINICAL_LANE_MAP[$sample]-}" ]]; then
        echo "[$(date)] [${group_label}] Merging lanes for ${sample}: ${CLINICAL_LANE_MAP[$sample]}" | tee -a "$BATCH_LOG"
    fi

    if ! pair="$(resolve_pair_for_sample "$group_label" "$data_root" "$sample" "$SEARCH_DEPTH")"; then
        echo "[$(date)] [${group_label}] Missing FASTQ pair for ${sample}; skipping" | tee -a "$BATCH_LOG"
        return 10
    fi

    IFS=$'\t' read -r r1_file r2_file <<< "$pair"
    sample_output="${output_base}/${sample}"
    sample_log="${output_base}/mapper_logs/${sample}.log"

    cmd=(
        python paradism.py
        --read1 "$r1_file"
        --read2 "$r2_file"
        --reference "$REFERENCE"
        --aligner "$ALIGNER"
        --threads "$THREADS"
        --output-dir "$sample_output"
        --iterations "$ITERATIONS"
    )

    if [[ "$ADD_QUALITY_FILTERS" -eq 1 ]]; then
        cmd+=(
            --add-quality-filters
            --qual-threshold "$QUAL_THRESHOLD"
            --dp-threshold "$DP_THRESHOLD"
            --af-threshold "$AF_THRESHOLD"
        )
    fi

    echo "[$(date)] [${group_label}] Processing ${sample}" | tee -a "$BATCH_LOG"
    echo "[$(date)] [${group_label}] R1=${r1_file}" | tee -a "$BATCH_LOG"
    echo "[$(date)] [${group_label}] R2=${r2_file}" | tee -a "$BATCH_LOG"
    echo "[$(date)] [${group_label}] CMD=${cmd[*]}" | tee -a "$BATCH_LOG"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        return 0
    fi

    if "${cmd[@]}" >> "$sample_log" 2>&1; then
        echo "[$(date)] [${group_label}] Completed ${sample}" | tee -a "$BATCH_LOG"
        return 0
    fi

    echo "[$(date)] [${group_label}] Failed ${sample}" | tee -a "$BATCH_LOG"
    return 20
}

run_group() {
    local group_label="$1"
    local data_root="$2"
    local output_base="$3"
    shift 3
    local samples=("$@")

    local processed=0
    local success=0
    local failed=0
    local missing=0
    local sample
    local -a pending_samples=()

    echo "[$(date)] Starting ${group_label} batch" | tee -a "$BATCH_LOG"
    echo "[$(date)] ${group_label} data directory: ${data_root}" | tee -a "$BATCH_LOG"
    echo "[$(date)] ${group_label} parallel samples: ${PARALLEL_SAMPLES}" | tee -a "$BATCH_LOG"

    run_pending_chunk() {
        local -a pids=()
        local -a pid_samples=()
        local idx pid rc

        for sample in "${pending_samples[@]}"; do
            process_sample "$group_label" "$data_root" "$output_base" "$sample" &
            pids+=("$!")
            pid_samples+=("$sample")
        done

        for idx in "${!pids[@]}"; do
            pid="${pids[$idx]}"
            if wait "$pid"; then
                rc=0
            else
                rc=$?
            fi

            case "$rc" in
                0)
                    success=$((success + 1))
                    ;;
                10)
                    missing=$((missing + 1))
                    ;;
                *)
                    failed=$((failed + 1))
                    ;;
            esac
        done
    }

    for sample in "${samples[@]}"; do
        processed=$((processed + 1))
        pending_samples+=("$sample")

        if (( ${#pending_samples[@]} >= PARALLEL_SAMPLES )); then
            run_pending_chunk
            pending_samples=()
        fi
    done

    if (( ${#pending_samples[@]} > 0 )); then
        run_pending_chunk
        pending_samples=()
    fi

    TOTAL_PROCESSED=$((TOTAL_PROCESSED + processed))
    TOTAL_SUCCESS=$((TOTAL_SUCCESS + success))
    TOTAL_FAILED=$((TOTAL_FAILED + failed))
    TOTAL_MISSING=$((TOTAL_MISSING + missing))

    echo "[$(date)] ${group_label} summary: total=${processed}, success=${success}, failed=${failed}, missing=${missing}" | tee -a "$BATCH_LOG"
}

echo "========================================="
echo "Running ParaDISM on HTS clinical + controls"
echo "Reference: $REFERENCE"
echo "Aligner: $ALIGNER"
echo "Threads: $THREADS"
echo "Iterations: $ITERATIONS"
echo "Quality filters enabled: $ADD_QUALITY_FILTERS"
if [[ "$ADD_QUALITY_FILTERS" -eq 1 ]]; then
    echo "Quality thresholds: QUAL>=$QUAL_THRESHOLD DP>=$DP_THRESHOLD AF>=$AF_THRESHOLD"
fi
echo "Clinical data: $CLINICAL_DATA_DIR"
echo "Control data: $CONTROL_DATA_DIR"
echo "Output root: $OUT_ROOT"
echo "Parallel samples: $PARALLEL_SAMPLES"
echo "Dry run: $DRY_RUN"
echo "Batch log: $BATCH_LOG"
echo "========================================="

echo "[$(date)] Configuration loaded" | tee -a "$BATCH_LOG"
echo "[$(date)] Clinical output base: $CLINICAL_OUTPUT_BASE" | tee -a "$BATCH_LOG"
echo "[$(date)] Control output base: $CONTROL_OUTPUT_BASE" | tee -a "$BATCH_LOG"

run_group "clinical" "$CLINICAL_DATA_DIR" "$CLINICAL_OUTPUT_BASE" "${CLINICAL_SAMPLES[@]}"
run_group "control" "$CONTROL_DATA_DIR" "$CONTROL_OUTPUT_BASE" "${CONTROL_SAMPLES[@]}"

echo "========================================="
echo "Final summary: processed=$TOTAL_PROCESSED success=$TOTAL_SUCCESS failed=$TOTAL_FAILED missing=$TOTAL_MISSING"
echo "Done. See $BATCH_LOG"
echo "========================================="
