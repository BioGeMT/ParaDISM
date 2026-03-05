#!/usr/bin/env bash
# Export coverage metrics from IGV-ready per-sample BAM files.
# Produces:
# - per_contig_coverage.tsv  (one row per sample-contig)
# - per_sample_coverage.tsv  (one row per sample, aggregated across contigs)
#
# Example:
#   bash hts_analysis/export_coverage_metrics_hts.sh \
#     --igv-dir hts_analysis/HTS_rerun_bowtie2_t4_p10/igv_per_sample

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash hts_analysis/export_coverage_metrics_hts.sh [options]

Options:
  --igv-dir PATH      IGV per-sample directory with clinical/ and control/
                      (default: hts_analysis/HTS_rerun_bowtie2_t4_p10/igv_per_sample)
  --out-dir PATH      Output directory for coverage TSVs
                      (default: <igv-dir>/coverage_metrics)
  --min-mapq N        Minimum mapping quality passed to samtools coverage (default: 0)
  --min-baseq N       Minimum base quality passed to samtools coverage (default: 0)
  -h, --help          Show this help
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

IGV_DIR="hts_analysis/HTS_rerun_bowtie2_t4_p10/igv_per_sample"
OUT_DIR=""
MIN_MAPQ=0
MIN_BASEQ=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --igv-dir)
            IGV_DIR="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --min-mapq)
            MIN_MAPQ="$2"
            shift 2
            ;;
        --min-baseq)
            MIN_BASEQ="$2"
            shift 2
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

if [[ -z "$OUT_DIR" ]]; then
    OUT_DIR="${IGV_DIR}/coverage_metrics"
fi

if ! [[ "$MIN_MAPQ" =~ ^[0-9]+$ ]]; then
    echo "Error: --min-mapq must be a non-negative integer (got: $MIN_MAPQ)"
    exit 1
fi
if ! [[ "$MIN_BASEQ" =~ ^[0-9]+$ ]]; then
    echo "Error: --min-baseq must be a non-negative integer (got: $MIN_BASEQ)"
    exit 1
fi

SAMTOOLS="$(command -v samtools || true)"
if [[ -z "$SAMTOOLS" && -x "$HOME/miniconda3/envs/paradism_env/bin/samtools" ]]; then
    SAMTOOLS="$HOME/miniconda3/envs/paradism_env/bin/samtools"
fi
if [[ -z "$SAMTOOLS" ]]; then
    echo "Error: samtools not found in PATH and not found in paradism_env."
    exit 1
fi

CLINICAL_DIR="${IGV_DIR}/clinical"
CONTROL_DIR="${IGV_DIR}/control"
for d in "$CLINICAL_DIR" "$CONTROL_DIR"; do
    if [[ ! -d "$d" ]]; then
        echo "Error: expected directory not found: $d"
        exit 1
    fi
done

mkdir -p "$OUT_DIR"

PER_CONTIG_TSV="${OUT_DIR}/per_contig_coverage.tsv"
PER_SAMPLE_TSV="${OUT_DIR}/per_sample_coverage.tsv"
STATUS_TSV="${OUT_DIR}/coverage_status.tsv"

echo -e "group\tsample\tcontig\tstart\tend\tlength\tnumreads\tcovbases\tcoverage_pct\tmean_depth\tmean_baseq\tmean_mapq\tbam_path" > "$PER_CONTIG_TSV"
echo -e "group\tsample\tcontigs\ttotal_bases\tnumreads\tcovbases\tcoverage_pct\tmean_depth\tmean_baseq\tmean_mapq\tbam_path" > "$PER_SAMPLE_TSV"
echo -e "group\tsample\tok\tbam_path\treason" > "$STATUS_TSV"

TOTAL=0
OK=0
FAIL=0

process_group() {
    local group_label="$1"
    local group_dir="$2"
    local bam sample tmp_cov reason

    while IFS= read -r bam; do
        sample="$(basename "$bam" .bam)"
        TOTAL=$((TOTAL + 1))
        tmp_cov="$(mktemp)"
        reason=""

        if "$SAMTOOLS" coverage -H -q "$MIN_MAPQ" -Q "$MIN_BASEQ" "$bam" > "$tmp_cov"; then
            if [[ -s "$tmp_cov" ]]; then
                awk -v grp="$group_label" -v smp="$sample" -v bam_path="$bam" '
                    BEGIN { OFS="\t" }
                    NF >= 9 {
                        len = $3 - $2 + 1
                        printf "%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                            grp, smp, $1, $2, $3, len, $4, $5, $6, $7, $8, $9, bam_path
                    }
                ' "$tmp_cov" >> "$PER_CONTIG_TSV"

                awk -v grp="$group_label" -v smp="$sample" -v bam_path="$bam" '
                    BEGIN {
                        OFS="\t"
                        contigs=0; total_len=0; total_reads=0; total_covbases=0
                        depth_weighted=0; baseq_weighted=0; mapq_weighted=0
                    }
                    NF >= 9 {
                        len = $3 - $2 + 1
                        contigs += 1
                        total_len += len
                        total_reads += $4
                        total_covbases += $5
                        depth_weighted += ($7 * len)
                        baseq_weighted += ($8 * $5)
                        mapq_weighted += ($9 * $4)
                    }
                    END {
                        if (total_len > 0) {
                            cov_pct = (100.0 * total_covbases) / total_len
                            mean_depth = depth_weighted / total_len
                        } else {
                            cov_pct = 0
                            mean_depth = 0
                        }
                        if (total_covbases > 0) {
                            mean_baseq = baseq_weighted / total_covbases
                        } else {
                            mean_baseq = 0
                        }
                        if (total_reads > 0) {
                            mean_mapq = mapq_weighted / total_reads
                        } else {
                            mean_mapq = 0
                        }

                        printf "%s\t%s\t%d\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%s\n",
                            grp, smp, contigs, total_len, total_reads, total_covbases,
                            cov_pct, mean_depth, mean_baseq, mean_mapq, bam_path
                    }
                ' "$tmp_cov" >> "$PER_SAMPLE_TSV"

                OK=$((OK + 1))
                printf "%s\t%s\t1\t%s\t%s\n" "$group_label" "$sample" "$bam" "" >> "$STATUS_TSV"
            else
                FAIL=$((FAIL + 1))
                reason="no_coverage_rows"
                printf "%s\t%s\t0\t%s\t%s\n" "$group_label" "$sample" "$bam" "$reason" >> "$STATUS_TSV"
            fi
        else
            FAIL=$((FAIL + 1))
            reason="samtools_coverage_failed"
            printf "%s\t%s\t0\t%s\t%s\n" "$group_label" "$sample" "$bam" "$reason" >> "$STATUS_TSV"
        fi

        rm -f "$tmp_cov"
    done < <(find "$group_dir" -maxdepth 1 -type f -name '*.bam' | sort)
}

echo "=========================================="
echo "Exporting coverage metrics"
echo "IGV dir: $IGV_DIR"
echo "Out dir: $OUT_DIR"
echo "samtools: $SAMTOOLS"
echo "Min MAPQ: $MIN_MAPQ"
echo "Min BASEQ: $MIN_BASEQ"
echo "=========================================="

process_group "clinical" "$CLINICAL_DIR"
process_group "control" "$CONTROL_DIR"

echo ""
echo "Done."
echo "Samples processed: $TOTAL"
echo "Coverage success: $OK"
echo "Coverage failed: $FAIL"
echo "Per-contig TSV: $PER_CONTIG_TSV"
echo "Per-sample TSV: $PER_SAMPLE_TSV"
echo "Status TSV: $STATUS_TSV"
