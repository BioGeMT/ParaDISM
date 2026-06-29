#!/usr/bin/env bash
# Run ParaDISM on GIAB HG002 with G60 threshold.
# Uses bowtie2, minalt 5, qual filtered
# Run from anywhere: bash benchmark/giab/run_giab.sh

set -euo pipefail

CALLER_CWD="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

READS_DIR="${READS_DIR:-$SCRIPT_DIR/giab_hg002_reads}"
REFERENCE="${REFERENCE:-$PROJECT_ROOT/benchmark/references/pkd1_panel.fa}"
PARADISM="$PROJECT_ROOT/paradism.py"
THREADS="${THREADS:-8}"
WORKERS="${WORKERS:-1}"
ITERATIONS="${ITERATIONS:-1}"
MIN_ALT_COUNT="${MIN_ALT_COUNT:-5}"
OUTPUT_G60="${OUTPUT_G60:-$SCRIPT_DIR/giab_hg002_output_bowtie2_G60_min5_qfilters}"

usage() {
    cat <<'EOF'
Usage:
  bash benchmark/giab/run_giab.sh [options]

Options:
  --reads-dir DIR       Directory containing HG002_R1/HG002_R2 FASTQ files
  --reference FILE      ParaDISM reference FASTA (default: benchmark/references/pkd1_panel.fa)
  --output-dir DIR      Output directory (default: benchmark/giab/giab_hg002_output_bowtie2_G60_min5_qfilters)
  --threads N           Threads for ParaDISM (default: 8 or env THREADS)
  --workers N           Read-assignment worker processes (default: 1 or env WORKERS)
  --iterations N        ParaDISM iterations (default: 1 or env ITERATIONS)
  --min-alt-count N     ParaDISM --min-alternate-count (default: 5 or env MIN_ALT_COUNT)
  -h, --help            Show help
EOF
}

to_abs_path() {
    local path="$1"
    if [[ "$path" == /* ]]; then
        printf "%s\n" "$path"
    else
        printf "%s\n" "$CALLER_CWD/$path"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --reads-dir)
            READS_DIR="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_G60="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --workers)
            WORKERS="$2"
            shift 2
            ;;
        --iterations)
            ITERATIONS="$2"
            shift 2
            ;;
        --min-alt-count)
            MIN_ALT_COUNT="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

READS_DIR="$(to_abs_path "$READS_DIR")"
REFERENCE="$(to_abs_path "$REFERENCE")"
OUTPUT_G60="$(to_abs_path "$OUTPUT_G60")"

if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: reference FASTA not found: $REFERENCE" >&2
    exit 1
fi

# Find merged or downsampled reads.
R1_MERGED=""
R2_MERGED=""
for extension in fq.gz fastq.gz fq fastq; do
    candidate_r1="${READS_DIR}/HG002_R1.${extension}"
    candidate_r2="${READS_DIR}/HG002_R2.${extension}"
    if [[ -f "$candidate_r1" && -f "$candidate_r2" ]]; then
        R1_MERGED="$candidate_r1"
        R2_MERGED="$candidate_r2"
        break
    fi
done

if [[ -z "$R1_MERGED" ]]; then
    echo "Error: HG002 read pair not found in $READS_DIR"
    echo "Expected: HG002_R1/HG002_R2 with extension .fq.gz, .fastq.gz, .fq, or .fastq"
    exit 1
fi

echo "=========================================="
echo "ParaDISM GIAB HG002 - G60"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Read1: $R1_MERGED"
echo "  Read2: $R2_MERGED"
echo "  Reference: $REFERENCE"
echo "  Threads: $THREADS"
echo "  Workers: $WORKERS"
echo "  Iterations: $ITERATIONS"
echo "  Min-alternate-count: $MIN_ALT_COUNT"
echo ""

# Use the pinned environment if available, but also support already-active envs.
if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate paradism 2>/dev/null || true
fi

# Run G60 (recommended threshold)
if [[ -d "$OUTPUT_G60/final_outputs" ]]; then
    echo "G60 already complete, skipping..."
else
    echo "=========================================="
    echo "Running G60 threshold..."
    echo "=========================================="
    python "$PARADISM" \
        --read1 "$R1_MERGED" \
        --read2 "$R2_MERGED" \
        --reference "$REFERENCE" \
        --aligner bowtie2 \
        --threads "$THREADS" \
        --workers "$WORKERS" \
        --iterations "$ITERATIONS" \
        --min-alternate-count "$MIN_ALT_COUNT" \
        --add-quality-filters \
        --qual-threshold 20 \
        --dp-threshold 10 \
        --af-threshold 0.05 \
        --threshold "G,60,60" \
        --output-dir "$OUTPUT_G60"
fi

echo ""
echo "=========================================="
echo "Complete!"
echo "=========================================="
echo ""
echo "Output directories:"
echo "  G60: $OUTPUT_G60"
echo ""
echo "Next: Run variant calling"
echo "  bash benchmark/giab/run_variant_calling_to_vcf_out.sh --run-dir \"$OUTPUT_G60\" --out-dir benchmark/giab/vcf_out_full --threads \"$THREADS\""
