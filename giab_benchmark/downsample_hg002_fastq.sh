#!/usr/bin/env bash
# Build a HG002 subset FASTQ pair (no ParaDISM run), optimized for huge inputs.
#
# Strategy:
# - Uses pre-split shard files under giab_hg002_reads (D1_*_R1_*.fastq / R2).
# - Distributes target pairs across all shards by shard size (bytes), so all shards contribute.
# - Streams each shard and copies only its allocated number of pairs.
# - Avoids scanning the full merged HG002 FASTQs.
#
# Usage:
#   ./giab_benchmark/downsample_hg002_fastq.sh [options]
#
# Options:
#   --coverage N          Target coverage (default: 10)
#   --bp-per-pair N       Bases per pair (default: 500 for 2x250)
#   --genome-bp N         Genome size in bp (default: 3137300923)
#   --progress-every N    Progress interval in written pairs (default: 5000000)
#   --out-dir PATH        Output directory (default: auto in giab_hg002_reads)
#   --force               Overwrite existing output FASTQs
#   -h, --help            Show help

set -euo pipefail

print_help() {
  sed -n '1,25p' "$0"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
READS_DIR="$SCRIPT_DIR/giab_hg002_reads"

COVERAGE=10
BP_PER_PAIR=500
GENOME_BP=3137300923
PROGRESS_EVERY=5000000
FORCE=0
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --coverage)
      COVERAGE="$2"
      shift 2
      ;;
    --bp-per-pair)
      BP_PER_PAIR="$2"
      shift 2
      ;;
    --genome-bp)
      GENOME_BP="$2"
      shift 2
      ;;
    --progress-every)
      PROGRESS_EVERY="$2"
      shift 2
      ;;
    --out-dir)
      OUT_DIR="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    *)
      echo "Error: unknown option: $1" >&2
      print_help >&2
      exit 1
      ;;
  esac
done

for numeric_value in "$COVERAGE" "$BP_PER_PAIR" "$GENOME_BP" "$PROGRESS_EVERY"; do
  if ! [[ "$numeric_value" =~ ^[0-9]+$ ]]; then
    echo "Error: all numeric options must be positive integers." >&2
    exit 1
  fi
done

if (( BP_PER_PAIR == 0 )); then
  echo "Error: --bp-per-pair must be > 0" >&2
  exit 1
fi

TARGET_PAIRS=$(( (COVERAGE * GENOME_BP + BP_PER_PAIR - 1) / BP_PER_PAIR ))
if (( TARGET_PAIRS == 0 )); then
  echo "Error: computed target pairs is 0." >&2
  exit 1
fi

mapfile -t read1_shards < <(find "$READS_DIR" -maxdepth 1 -type f -name 'D1_*_R1_*.fastq' | LC_ALL=C sort)
if (( ${#read1_shards[@]} == 0 )); then
  echo "Error: no shard FASTQs found in $READS_DIR (expected D1_*_R1_*.fastq)." >&2
  exit 1
fi

declare -a read2_shards shard_bytes shard_quotas shard_remainders
total_shard_bytes=0

for read1_shard in "${read1_shards[@]}"; do
  read2_shard="${read1_shard/_R1_/_R2_}"
  if [[ ! -f "$read2_shard" ]]; then
    echo "Error: missing paired shard for $read1_shard" >&2
    echo "Expected: $read2_shard" >&2
    exit 1
  fi
  read2_shards+=("$read2_shard")

  bytes="$(stat -c '%s' "$read1_shard")"
  if ! [[ "$bytes" =~ ^[0-9]+$ ]]; then
    echo "Error: failed to read size for shard: $read1_shard" >&2
    exit 1
  fi
  shard_bytes+=("$bytes")
  total_shard_bytes=$(( total_shard_bytes + bytes ))
done

if (( total_shard_bytes == 0 )); then
  echo "Error: total shard bytes is 0." >&2
  exit 1
fi

if [[ -z "$OUT_DIR" ]]; then
  OUT_DIR="${READS_DIR}/subset_genome${COVERAGE}x_bp${BP_PER_PAIR}_balanced"
fi

OUTPUT_READ1="${OUT_DIR}/HG002_R1.fastq"
OUTPUT_READ2="${OUT_DIR}/HG002_R2.fastq"
META_OUT="${OUT_DIR}/downsample_metadata.txt"
SELECTED_SHARDS_OUT="${OUT_DIR}/selected_shards.txt"

mkdir -p "$OUT_DIR"

if [[ -s "$OUTPUT_READ1" || -s "$OUTPUT_READ2" ]]; then
  if (( FORCE == 0 )); then
    echo "Error: output exists. Use --force to overwrite." >&2
    echo "  $OUTPUT_READ1" >&2
    echo "  $OUTPUT_READ2" >&2
    exit 1
  fi
fi

if (( FORCE == 1 )); then
  rm -f "$OUTPUT_READ1" "$OUTPUT_READ2" "$META_OUT" "$SELECTED_SHARDS_OUT"
fi

awk_bin="$(command -v mawk || true)"
if [[ -z "$awk_bin" ]]; then
  awk_bin="$(command -v gawk || true)"
fi
if [[ -z "$awk_bin" ]]; then
  awk_bin="$(command -v awk || true)"
fi
if [[ -z "$awk_bin" ]]; then
  echo "Error: awk not found in PATH." >&2
  exit 1
fi

quota_sum=0
for i in "${!read1_shards[@]}"; do
  bytes="${shard_bytes[$i]}"
  quota=$(( TARGET_PAIRS * bytes / total_shard_bytes ))
  remainder=$(( TARGET_PAIRS * bytes % total_shard_bytes ))
  shard_quotas[$i]="$quota"
  shard_remainders[$i]="$remainder"
  quota_sum=$(( quota_sum + quota ))
done

leftover=$(( TARGET_PAIRS - quota_sum ))
if (( leftover > 0 )); then
  mapfile -t ranked_indices < <(
    for i in "${!read1_shards[@]}"; do
      printf "%s\t%s\n" "$i" "${shard_remainders[$i]}"
    done | sort -t $'\t' -k2,2nr -k1,1n | awk -F $'\t' '{print $1}'
  )

  for (( j=0; j<leftover; j++ )); do
    idx="${ranked_indices[$j]}"
    shard_quotas[$idx]=$(( shard_quotas[$idx] + 1 ))
  done
fi

temp_read1="${OUTPUT_READ1}.tmp.$$"
temp_read2="${OUTPUT_READ2}.tmp.$$"
rm -f "$temp_read1" "$temp_read2"

cleanup_tmp() {
  rm -f "$temp_read1" "$temp_read2"
}
trap cleanup_tmp EXIT

remaining_pairs="$TARGET_PAIRS"
written_pairs_total=0
selected_shards_count=0

echo "Subset configuration:"
echo "  Coverage:          ${COVERAGE}x"
echo "  Genome bp:         $GENOME_BP"
echo "  BP per pair:       $BP_PER_PAIR"
echo "  Target pairs:      $TARGET_PAIRS"
echo "  Shard files found: ${#read1_shards[@]}"
echo "  Strategy:          balanced_across_all_shards"
echo "  Output R1:         $OUTPUT_READ1"
echo "  Output R2:         $OUTPUT_READ2"
echo "  Awk binary:        $awk_bin"

for i in "${!read1_shards[@]}"; do
  if (( remaining_pairs == 0 )); then
    break
  fi

  read1_shard="${read1_shards[$i]}"
  read2_shard="${read2_shards[$i]}"
  requested_pairs="${shard_quotas[$i]}"

  if (( i == ${#read1_shards[@]} - 1 )); then
    requested_pairs="$remaining_pairs"
  fi
  if (( requested_pairs > remaining_pairs )); then
    requested_pairs="$remaining_pairs"
  fi
  if (( requested_pairs == 0 )); then
    continue
  fi

  copied_pairs="$(
    "$awk_bin" \
      -v read1_path="$read1_shard" \
      -v read2_path="$read2_shard" \
      -v output_read1="$temp_read1" \
      -v output_read2="$temp_read2" \
      -v max_pairs="$requested_pairs" \
      -v progress_every="$PROGRESS_EVERY" \
      -v already_written="$written_pairs_total" '
BEGIN{
  copied=0;
  while (copied < max_pairs && (getline header1 < read1_path) > 0) {
    if ((getline sequence1 < read1_path) <= 0 || (getline plus1 < read1_path) <= 0 || (getline quality1 < read1_path) <= 0) exit 2;
    if ((getline header2 < read2_path) <= 0 || (getline sequence2 < read2_path) <= 0 || (getline plus2 < read2_path) <= 0 || (getline quality2 < read2_path) <= 0) exit 3;

    printf "%s\n%s\n%s\n%s\n", header1, sequence1, plus1, quality1 >> output_read1;
    printf "%s\n%s\n%s\n%s\n", header2, sequence2, plus2, quality2 >> output_read2;

    copied++;
    total_written = already_written + copied;
    if (progress_every > 0 && total_written % progress_every == 0) {
      printf "  written_pairs=%d\n", total_written > "/dev/stderr";
    }
  }
  print copied;
}' /dev/null
  )"

  if ! [[ "$copied_pairs" =~ ^[0-9]+$ ]]; then
    echo "Error: failed to parse copied pairs for shard $read1_shard" >&2
    exit 1
  fi

  remaining_pairs=$(( remaining_pairs - copied_pairs ))
  written_pairs_total=$(( written_pairs_total + copied_pairs ))

  if (( copied_pairs > 0 )); then
    selected_shards_count=$(( selected_shards_count + 1 ))
    printf "%s\trequested_pairs=%s\tcopied_pairs=%s\n" "$(basename "$read1_shard")" "$requested_pairs" "$copied_pairs" >> "$SELECTED_SHARDS_OUT"
    echo "  shard: $(basename "$read1_shard") requested_pairs=$requested_pairs copied_pairs=$copied_pairs remaining_pairs=$remaining_pairs"
  fi
done

if (( remaining_pairs > 0 )); then
  echo "Error: not enough pairs in available shard files." >&2
  echo "  missing_pairs=$remaining_pairs" >&2
  exit 1
fi

mv -f "$temp_read1" "$OUTPUT_READ1"
mv -f "$temp_read2" "$OUTPUT_READ2"

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "strategy=balanced_across_all_shards"
  echo "reads_dir=$READS_DIR"
  echo "coverage=$COVERAGE"
  echo "genome_bp=$GENOME_BP"
  echo "bp_per_pair=$BP_PER_PAIR"
  echo "target_pairs=$TARGET_PAIRS"
  echo "written_pairs=$written_pairs_total"
  echo "progress_every=$PROGRESS_EVERY"
  echo "total_shard_bytes=$total_shard_bytes"
  echo "selected_shards_count=$selected_shards_count"
  echo "read1_out=$OUTPUT_READ1"
  echo "read2_out=$OUTPUT_READ2"
  echo "selected_shards_file=$SELECTED_SHARDS_OUT"
} > "$META_OUT"

trap - EXIT

echo "Subset creation complete."
echo "  Written pairs: $written_pairs_total"
echo "  R1: $OUTPUT_READ1"
echo "  R2: $OUTPUT_READ2"
echo "  Metadata: $META_OUT"
echo "  Shards list: $SELECTED_SHARDS_OUT"
