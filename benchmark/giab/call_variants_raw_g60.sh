#!/usr/bin/env bash
# Call variants from a ParaDISM run directory.
# Produces raw VCFs first, then filtered simple biallelic A/C/G/T SNP VCFs
# so raw vs filtered can be compared directly.
#
# Run from anywhere:
#   bash benchmark/giab/call_variants_raw_g60.sh \
#     --input-dir benchmark/giab/giab_hg002_output_bowtie2_G60_min5_qfilters \
#     --output-dir benchmark/giab/giab_hg002_output_bowtie2_G60_min5_qfilters/variant_calling_custom

set -euo pipefail

CALLER_CWD="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

REFERENCE="$PROJECT_ROOT/benchmark/references/pkd1_panel.fa"
MIN_ALT_COUNT="${MIN_ALT_COUNT:-}"
ATOMIZE_SCRIPT="$SCRIPT_DIR/atomize_equal_length_substitutions.py"
THREADS="${THREADS:-8}"
INPUT_DIR=""
OUTPUT_DIR=""
BAM_PREFIX=""

usage() {
    cat <<'EOF'
Usage:
  bash benchmark/giab/call_variants_raw_g60.sh --input-dir DIR --output-dir DIR [options]

Required:
  --input-dir DIR       ParaDISM run output directory (contains final_outputs and iteration_1)
  --output-dir DIR      Destination directory for variant calling outputs

Optional:
  --reference FILE      Reference fasta (default: benchmark/references/pkd1_panel.fa)
  --bam-prefix PREFIX   Prefix used in per-gene BAM names (default: basename of --input-dir)
  --threads N           Threads for samtools sort (default: 8 or env THREADS)
  -h, --help            Show help

Environment:
  MIN_ALT_COUNT         If set, passed to FreeBayes as --min-alternate-count
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --reference)
            REFERENCE="$2"
            shift 2
            ;;
        --bam-prefix)
            BAM_PREFIX="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
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

to_abs_path() {
    local path="$1"
    if [[ "$path" == /* ]]; then
        printf "%s\n" "$path"
    else
        printf "%s\n" "$CALLER_CWD/$path"
    fi
}

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: --input-dir and --output-dir are required." >&2
    usage >&2
    exit 1
fi

INPUT_DIR="$(to_abs_path "$INPUT_DIR")"
OUTPUT_DIR="$(to_abs_path "$OUTPUT_DIR")"
if [[ "$REFERENCE" != /* ]]; then
    REFERENCE="$(to_abs_path "$REFERENCE")"
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: input directory not found: $INPUT_DIR" >&2
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

if [[ -z "$BAM_PREFIX" ]]; then
    BAM_PREFIX="$(basename "$INPUT_DIR")"
fi

if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate paradism 2>/dev/null || true
fi

for required_tool in freebayes bcftools samtools; do
    if ! command -v "$required_tool" >/dev/null 2>&1; then
        echo "Error: required tool not found on PATH: $required_tool" >&2
        echo "Install the pinned environment with: conda env create -f environment.yml" >&2
        exit 1
    fi
done

GENES=()
while IFS= read -r gene_name; do
    GENES+=("$gene_name")
done < <(grep "^>" "$REFERENCE" | sed 's/^>//')

FREEBAYES_ARGS=(--ploidy 2)
SNP_FILTER_ARGS=(-v snps -m2 -M2)
SNP_ACGT_EXPR='REF~"^[ACGT]$" && ALT~"^[ACGT]$"'
if [[ -n "$MIN_ALT_COUNT" ]]; then
    if ! [[ "$MIN_ALT_COUNT" =~ ^[0-9]+$ ]]; then
        echo "Error: MIN_ALT_COUNT must be a non-negative integer (got: $MIN_ALT_COUNT)" >&2
        exit 1
    fi
    FREEBAYES_ARGS+=(--min-alternate-count "$MIN_ALT_COUNT")
fi

if [[ ! -f "$ATOMIZE_SCRIPT" ]]; then
    echo "Error: atomize script not found: $ATOMIZE_SCRIPT" >&2
    exit 1
fi

call_raw_from_gene_bams() {
    local bam_dir=$1
    local out_dir=$2
    local prefix=$3

    echo "Calling RAW variants per-gene: $bam_dir"
    mkdir -p "$out_dir/per_gene"

    local gene_vcfs_raw=()
    local gene_vcfs_filtered=()
    for gene in "${GENES[@]}"; do
        local gene_bam="${bam_dir}/${prefix}_${gene}.sorted.bam"
        local gene_vcf_raw="${out_dir}/per_gene/${gene}.raw.vcf"
        local gene_vcfgz_raw="${gene_vcf_raw}.gz"
        local gene_vcf_filtered="${out_dir}/per_gene/${gene}.simple_snps_acgt.vcf"
        local gene_vcfgz_filtered="${gene_vcf_filtered}.gz"
        if [[ -f "$gene_bam" ]]; then
            echo "  FreeBayes raw: $gene"
            freebayes --bam "$gene_bam" --fasta-reference "$REFERENCE" \
                "${FREEBAYES_ARGS[@]}" > "$gene_vcf_raw"
            bcftools sort -Oz -o "$gene_vcfgz_raw" "$gene_vcf_raw"
            bcftools index -f "$gene_vcfgz_raw"
            gene_vcfs_raw+=("$gene_vcfgz_raw")

            bcftools norm -m -any -O v "$gene_vcf_raw" | \
                python3 "$ATOMIZE_SCRIPT" | \
                bcftools view "${SNP_FILTER_ARGS[@]}" -i "$SNP_ACGT_EXPR" -O v -o "$gene_vcf_filtered"
            bcftools sort -Oz -o "$gene_vcfgz_filtered" "$gene_vcf_filtered"
            bcftools index -f "$gene_vcfgz_filtered"
            gene_vcfs_filtered+=("$gene_vcfgz_filtered")
        fi
    done

    if [[ ${#gene_vcfs_raw[@]} -eq 0 ]]; then
        echo "Error: no expected per-gene BAMs found in $bam_dir" >&2
        return 1
    fi

    echo "  Merging ${#gene_vcfs_raw[@]} per-gene raw VCFs..."
    bcftools concat -a -O v -o "${out_dir}/variants_raw.vcf" "${gene_vcfs_raw[@]}"
    bcftools sort -Oz -o "${out_dir}/variants_raw.vcf.gz" "${out_dir}/variants_raw.vcf"
    bcftools index -f "${out_dir}/variants_raw.vcf.gz"

    echo "  Merging ${#gene_vcfs_filtered[@]} per-gene filtered VCFs..."
    bcftools concat -a -O v -o "${out_dir}/variants_simple_snps_acgt.vcf" "${gene_vcfs_filtered[@]}"
    bcftools sort -Oz -o "${out_dir}/variants_simple_snps_acgt.vcf.gz" "${out_dir}/variants_simple_snps_acgt.vcf"
    bcftools index -f "${out_dir}/variants_simple_snps_acgt.vcf.gz"
}

call_raw_from_alignment() {
    local alignment_file=$1
    local out_dir=$2

    echo "Calling RAW variants from original alignment: $alignment_file"
    mkdir -p "$out_dir"

    if [[ ! -f "$alignment_file" ]]; then
        echo "  Alignment file not found: $alignment_file"
        return
    fi

    local sorted_bam="${out_dir}/mapped_reads.sorted.bam"
    if [[ ! -f "$sorted_bam" ]]; then
        echo "  Sorting base-alignment SAM/BAM..."
        samtools sort -@ "$THREADS" -o "$sorted_bam" "$alignment_file"
        samtools index "$sorted_bam"
    fi

    echo "  FreeBayes raw..."
    freebayes --bam "$sorted_bam" --fasta-reference "$REFERENCE" \
        "${FREEBAYES_ARGS[@]}" > "${out_dir}/variants_raw.vcf"

    bcftools sort -Oz -o "${out_dir}/variants_raw.vcf.gz" "${out_dir}/variants_raw.vcf"
    bcftools index -f "${out_dir}/variants_raw.vcf.gz"

    echo "  Filtering to simple biallelic A/C/G/T SNPs..."
    bcftools norm -m -any -O v "${out_dir}/variants_raw.vcf" | \
        python3 "$ATOMIZE_SCRIPT" | \
        bcftools view "${SNP_FILTER_ARGS[@]}" -i "$SNP_ACGT_EXPR" \
            -O v -o "${out_dir}/variants_simple_snps_acgt.vcf"
    bcftools sort -Oz -o "${out_dir}/variants_simple_snps_acgt.vcf.gz" "${out_dir}/variants_simple_snps_acgt.vcf"
    bcftools index -f "${out_dir}/variants_simple_snps_acgt.vcf.gz"
}

echo "=========================================="
echo "Variant Calling (raw + simple SNP filtered)"
echo "=========================================="
echo ""
echo "Input dir: $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "BAM prefix: $BAM_PREFIX"
echo "Reference: $REFERENCE"
echo "Threads: $THREADS"
echo "FreeBayes args: ${FREEBAYES_ARGS[*]}"
echo "SNP filter args: ${SNP_FILTER_ARGS[*]}"
echo "SNP ACGT expr: $SNP_ACGT_EXPR"
echo "Atomize script: $ATOMIZE_SCRIPT"
echo ""

PARADISM_BAM_DIR="$INPUT_DIR/final_outputs/${BAM_PREFIX}_bam"
PARADISM_OUT_DIR="$OUTPUT_DIR/paradism_raw"
BASE_ALIGNMENT="$INPUT_DIR/iteration_1/mapped_reads.sam"
if [[ ! -f "$BASE_ALIGNMENT" ]]; then
    BASE_ALIGNMENT="$INPUT_DIR/iteration_1/mapped_reads.bam"
fi
BASE_OUT_DIR="$OUTPUT_DIR/basealigner_raw"

if [[ ! -d "$PARADISM_BAM_DIR" ]]; then
    echo "Error: ParaDISM BAM directory not found: $PARADISM_BAM_DIR" >&2
    echo "The ParaDISM run is incomplete; do not start GIAB post-processing." >&2
    exit 1
fi

if [[ ! -f "$BASE_ALIGNMENT" ]]; then
    echo "Error: base aligner SAM/BAM not found in: $INPUT_DIR/iteration_1" >&2
    echo "The ParaDISM run is incomplete; do not start GIAB post-processing." >&2
    exit 1
fi

echo "=== ParaDISM (raw per-gene) ==="
call_raw_from_gene_bams "$PARADISM_BAM_DIR" \
                        "$PARADISM_OUT_DIR" \
                        "$BAM_PREFIX"

echo "=== Base Aligner (raw) ==="
call_raw_from_alignment "$BASE_ALIGNMENT" "$BASE_OUT_DIR"

echo ""
echo "Done!"
