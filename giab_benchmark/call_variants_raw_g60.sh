#!/usr/bin/env bash
# Call variants from a ParaDISM run directory.
# Produces raw VCFs first, then filtered simple biallelic A/C/G/T SNP VCFs
# so raw vs filtered can be compared directly.
#
# Run from anywhere:
#   bash giab_benchmark/call_variants_raw_g60.sh \
#     --input-dir giab_benchmark/giab_hg002_output_bowtie2_G60_min5_qfilters \
#     --output-dir giab_benchmark/giab_hg002_output_bowtie2_G60_min5_qfilters/variant_calling_custom

set -euo pipefail

CALLER_CWD="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

REFERENCE="$SCRIPT_DIR/ref.fa"
MIN_ALT_COUNT="${MIN_ALT_COUNT:-}"
ATOMIZE_SCRIPT="$SCRIPT_DIR/atomize_equal_length_substitutions.py"
THREADS="${THREADS:-8}"
INPUT_DIR=""
OUTPUT_DIR=""
BAM_PREFIX=""

usage() {
    cat <<'EOF'
Usage:
  bash giab_benchmark/call_variants_raw_g60.sh --input-dir DIR --output-dir DIR [options]

Required:
  --input-dir DIR       ParaDISM run output directory (contains final_outputs and/or iteration_1)
  --output-dir DIR      Destination directory for variant calling outputs

Optional:
  --reference FILE      Reference fasta (default: giab_benchmark/ref.fa)
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
if [[ "$REFERENCE" != "$SCRIPT_DIR/"* ]]; then
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

source ~/miniconda3/etc/profile.d/conda.sh
conda activate paradism_env

GENES=($(grep "^>" "$REFERENCE" | sed 's/^>//'))

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
        echo "  No gene BAMs found; skipping merge"
        return
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

call_raw_from_sam() {
    local sam_file=$1
    local out_dir=$2

    echo "Calling RAW variants from original SAM: $sam_file"
    mkdir -p "$out_dir"

    if [[ ! -f "$sam_file" ]]; then
        echo "  SAM file not found: $sam_file"
        return
    fi

    local sorted_bam="${out_dir}/mapped_reads.sorted.bam"
    if [[ ! -f "$sorted_bam" ]]; then
        echo "  Converting SAM to sorted BAM..."
        samtools view -bS "$sam_file" | samtools sort -@ "$THREADS" -o "$sorted_bam"
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
BASE_SAM="$INPUT_DIR/iteration_1/mapped_reads.sam"
BASE_OUT_DIR="$OUTPUT_DIR/basealigner_raw"

if [[ -d "$PARADISM_BAM_DIR" ]]; then
    echo "=== ParaDISM (raw per-gene) ==="
    call_raw_from_gene_bams "$PARADISM_BAM_DIR" \
                            "$PARADISM_OUT_DIR" \
                            "$BAM_PREFIX"
else
    echo "ParaDISM BAM directory not found, skipping: $PARADISM_BAM_DIR"
fi

if [[ -f "$BASE_SAM" ]]; then
    echo "=== Base Aligner (raw) ==="
    call_raw_from_sam "$BASE_SAM" "$BASE_OUT_DIR"
else
    echo "Base aligner SAM not found, skipping: $BASE_SAM"
fi

echo ""
echo "Done!"
