#!/usr/bin/env bash
# Call variants from GNAQ ParaDISM outputs.
# Runs FreeBayes per-gene on final BAMs, merges into one VCF per sample.
#
# Usage:
#   bash gnaq_analysis/call_variants.sh
#   bash gnaq_analysis/call_variants.sh --samples SRR5602384 SRR5602389
#   bash gnaq_analysis/call_variants.sh --run-dir gnaq_analysis/bwa-mem2_160_output

set -euo pipefail

CALLER_CWD="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

REFERENCE="$SCRIPT_DIR/gnaq-gnaqp_ref.fa"
RUN_DIR="$SCRIPT_DIR/bwa-mem2_160_output"
THREADS=8
MIN_ALT_COUNT=""
SAMPLES=()

usage() {
    cat <<'EOF'
Usage:
  bash gnaq_analysis/call_variants.sh [options]

Options:
  --run-dir DIR          ParaDISM run directory (default: gnaq_analysis/bwa-mem2_160_output)
  --reference FILE       Reference fasta (default: gnaq_analysis/gnaq-gnaqp_ref.fa)
  --samples S1 S2 ...    Specific samples to process (default: all subdirs in run-dir)
  --threads N            Threads for samtools sort (default: 8)
  --min-alt-count N      FreeBayes --min-alternate-count (default: not set)
  -h, --help             Show help
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
        --run-dir)
            RUN_DIR="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --samples)
            shift
            while [[ $# -gt 0 && ! "$1" == --* ]]; do
                SAMPLES+=("$1")
                shift
            done
            ;;
        --threads)
            THREADS="$2"
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

RUN_DIR="$(to_abs_path "$RUN_DIR")"
REFERENCE="$(to_abs_path "$REFERENCE")"

if [[ ! -d "$RUN_DIR" ]]; then
    echo "Error: run directory not found: $RUN_DIR" >&2
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "Error: reference not found: $REFERENCE" >&2
    exit 1
fi

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS <= 0 )); then
    echo "Error: --threads must be a positive integer (got: $THREADS)" >&2
    exit 1
fi

# Activate conda environment
if command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate paradism_env 2>/dev/null || true
fi

# Discover genes from reference
GENES=($(grep "^>" "$REFERENCE" | sed 's/^>//'))

# Discover samples if not specified
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    for d in "$RUN_DIR"/*/; do
        sample="$(basename "$d")"
        if [[ "$sample" != "mapper_logs" && -d "$d/final_outputs" ]]; then
            SAMPLES+=("$sample")
        fi
    done
fi

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "Error: no samples found in $RUN_DIR" >&2
    exit 1
fi

# FreeBayes args
FREEBAYES_ARGS=(--ploidy 2)
if [[ -n "$MIN_ALT_COUNT" ]]; then
    FREEBAYES_ARGS+=(--min-alternate-count "$MIN_ALT_COUNT")
fi

# Quality filter (same as GIAB benchmark)
FILTER_EXPR='INFO/AO>=3 && QUAL>=10 && INFO/SAF>=1 && INFO/SAR>=1'

echo "=========================================="
echo "GNAQ Variant Calling"
echo "=========================================="
echo "Run dir:    $RUN_DIR"
echo "Reference:  $REFERENCE"
echo "Genes:      ${GENES[*]}"
echo "Samples:    ${SAMPLES[*]}"
echo "FreeBayes:  ${FREEBAYES_ARGS[*]}"
echo "Filter:     $FILTER_EXPR"
echo ""

for sample in "${SAMPLES[@]}"; do
    echo "--- $sample ---"

    bam_dir="$RUN_DIR/$sample/final_outputs/${sample}_bam"
    out_dir="$RUN_DIR/$sample/variant_calling"
    mkdir -p "$out_dir/per_gene"

    if [[ ! -d "$bam_dir" ]]; then
        echo "  BAM directory not found: $bam_dir — skipping"
        continue
    fi

    gene_vcfs=()

    for gene in "${GENES[@]}"; do
        gene_bam="${bam_dir}/${sample}_${gene}.sorted.bam"
        gene_vcf="${out_dir}/per_gene/${gene}.vcf"
        gene_vcfgz="${gene_vcf}.gz"

        if [[ ! -f "$gene_bam" ]]; then
            echo "  BAM not found for ${gene}, skipping"
            continue
        fi

        echo "  FreeBayes: ${gene}"
        freebayes --bam "$gene_bam" --fasta-reference "$REFERENCE" \
            "${FREEBAYES_ARGS[@]}" > "$gene_vcf"

        bcftools sort -Oz -o "$gene_vcfgz" "$gene_vcf"
        bcftools index -f "$gene_vcfgz"
        gene_vcfs+=("$gene_vcfgz")
    done

    if [[ ${#gene_vcfs[@]} -eq 0 ]]; then
        echo "  No gene BAMs found; skipping merge"
        continue
    fi

    # Merge per-gene VCFs into one raw VCF per sample
    sample_raw="${out_dir}/${sample}.raw.vcf.gz"
    sample_filtered="${out_dir}/${sample}.filtered.vcf.gz"
    echo "  Merging ${#gene_vcfs[@]} per-gene VCFs..."
    bcftools concat -a -O v "${gene_vcfs[@]}" | \
        bcftools sort -Oz -o "$sample_raw"
    bcftools index -f "$sample_raw"

    # Apply quality filters
    echo "  Filtering: $FILTER_EXPR"
    bcftools view -i "$FILTER_EXPR" -Oz -o "$sample_filtered" "$sample_raw"
    bcftools index -f "$sample_filtered"

    # Summary
    raw_count=$(bcftools view -H "$sample_raw" | wc -l)
    filtered_count=$(bcftools view -H "$sample_filtered" | wc -l)
    echo "  Raw variants: $raw_count"
    echo "  After filtering: $filtered_count"
    echo "  Output: $sample_filtered"
    echo ""
done

echo "Done!"
