#!/usr/bin/env bash
# Run the full HTS pipeline with bowtie2 only:
#   1) ParaDISM mapping/iterations for clinical + controls
#   2) RAW variant calling (final_outputs only)
#   3) Final filtering from RAW to per-sample final VCFs
#
# Example:
#   bash hts_analysis/run_full_hts_bowtie2.sh \
#     --reference /homes/dtzim01/ref.fa \
#     --out-root hts_analysis/HTS_rerun_bowtie2_full

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: bash hts_analysis/run_full_hts_bowtie2.sh [options]

Options:
  --reference PATH            Reference FASTA (default: ref.fa)
  --clinical-data-dir PATH    Clinical FASTQ dir
                              (default: /mnt/STORAGE-BioGeMT-01/pkd_data)
  --control-data-dir PATH     Control FASTQ dir (optional; auto-detect if omitted)
  --out-root PATH             Output root used by all steps
                              (default: hts_analysis/HTS_rerun_bowtie2_full)
  --threads N                 Threads for ParaDISM mapping (default: 8)
  --iterations N              ParaDISM iterations (default: 10)
  --search-depth N            FASTQ discovery max depth (default: 5)
  --parallel-samples N        Max samples to run in parallel (default: 5)
  --disable-quality-filters   Disable ParaDISM quality filters
                              (enabled by default)
  --qual-threshold N          QUAL threshold for ParaDISM filters (default: 20)
  --dp-threshold N            DP threshold for ParaDISM filters (default: 10)
  --af-threshold F            AF threshold for ParaDISM filters (default: 0.05)
  --dry-run                   Dry run mapping step only; skips variant steps
  --skip-mapping              Skip mapping step
  --skip-raw-calling          Skip raw variant calling step
  --skip-final-filter         Skip final filter step
  -h, --help                  Show this help
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

REFERENCE="ref.fa"
CLINICAL_DATA_DIR="/mnt/STORAGE-BioGeMT-01/pkd_data"
CONTROL_DATA_DIR=""
OUT_ROOT="hts_analysis/HTS_rerun_bowtie2_full"
THREADS=8
ITERATIONS=10
SEARCH_DEPTH=5
PARALLEL_SAMPLES=5
ADD_QUALITY_FILTERS=1
QUAL_THRESHOLD=20
DP_THRESHOLD=10
AF_THRESHOLD=0.05
DRY_RUN=0
SKIP_MAPPING=0
SKIP_RAW_CALLING=0
SKIP_FINAL_FILTER=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reference)
            REFERENCE="$2"
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
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --iterations)
            ITERATIONS="$2"
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
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --skip-mapping)
            SKIP_MAPPING=1
            shift
            ;;
        --skip-raw-calling)
            SKIP_RAW_CALLING=1
            shift
            ;;
        --skip-final-filter)
            SKIP_FINAL_FILTER=1
            shift
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

MAPPING_CMD=(
    bash hts_analysis/run_hts_and_controls.sh
    --reference "$REFERENCE"
    --aligner bowtie2
    --threads "$THREADS"
    --iterations "$ITERATIONS"
    --clinical-data-dir "$CLINICAL_DATA_DIR"
    --out-root "$OUT_ROOT"
    --search-depth "$SEARCH_DEPTH"
    --parallel-samples "$PARALLEL_SAMPLES"
)

if [[ -n "$CONTROL_DATA_DIR" ]]; then
    MAPPING_CMD+=(--control-data-dir "$CONTROL_DATA_DIR")
fi

if [[ "$ADD_QUALITY_FILTERS" -eq 0 ]]; then
    MAPPING_CMD+=(--disable-quality-filters)
else
    MAPPING_CMD+=(
        --qual-threshold "$QUAL_THRESHOLD"
        --dp-threshold "$DP_THRESHOLD"
        --af-threshold "$AF_THRESHOLD"
    )
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    MAPPING_CMD+=(--dry-run)
fi

RAW_CALL_CMD=(
    bash hts_analysis/call_variants_raw_hts.sh
    --run-root "$OUT_ROOT"
    --reference "$REFERENCE"
    --aligner bowtie2
)

FINAL_FILTER_CMD=(
    bash hts_analysis/filter_variants_hts.sh
    --run-root "$OUT_ROOT"
    --aligner bowtie2
)

echo "=========================================="
echo "HTS full pipeline (bowtie2 only)"
echo "Project root: $PROJECT_ROOT"
echo "Reference: $REFERENCE"
echo "Out root: $OUT_ROOT"
echo "Threads: $THREADS"
echo "Iterations: $ITERATIONS"
echo "Parallel samples: $PARALLEL_SAMPLES"
echo "Dry run: $DRY_RUN"
echo "Skip mapping: $SKIP_MAPPING"
echo "Skip raw calling: $SKIP_RAW_CALLING"
echo "Skip final filter: $SKIP_FINAL_FILTER"
echo "=========================================="

if [[ "$SKIP_MAPPING" -eq 0 ]]; then
    echo "[1/3] Mapping + ParaDISM iterations (bowtie2)"
    "${MAPPING_CMD[@]}"
else
    echo "[1/3] Skipped mapping step"
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "Dry run enabled: skipping raw calling and final filtering."
    exit 0
fi

if [[ "$SKIP_RAW_CALLING" -eq 0 ]]; then
    echo "[2/3] RAW variant calling (final outputs only)"
    "${RAW_CALL_CMD[@]}"
else
    echo "[2/3] Skipped raw variant calling"
fi

if [[ "$SKIP_FINAL_FILTER" -eq 0 ]]; then
    echo "[3/3] Final variant filtering"
    "${FINAL_FILTER_CMD[@]}"
else
    echo "[3/3] Skipped final filtering"
fi

echo ""
echo "Complete."
echo "Run root: $OUT_ROOT"
echo "Clinical root: $OUT_ROOT/HTS_output/HTS_bowtie2_output"
echo "Control root:  $OUT_ROOT/HTS_CNTRL/HTS_bowtie2_output"
echo "RAW summary:   $OUT_ROOT/variant_calling_raw_summary_bowtie2.tsv"
echo "Final summary: $OUT_ROOT/variant_filtering_summary_bowtie2.tsv"
