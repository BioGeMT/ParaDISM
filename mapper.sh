#!/bin/bash

SCRIPTS_DIR="./scripts"
OUTPUT_DIR="./output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Starting pipeline..."

# Initialize variables
R1=""
R2=""
REF=""

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r1)
      if [[ -z "$2" ]]; then
        echo "ERROR: -r1 requires a filename argument"
        exit 1
      fi
      R1="$2"
      shift 2
      ;;
    -r2)
      if [[ -z "$2" ]]; then
        echo "ERROR: -r2 requires a filename argument"
        exit 1
      fi
      R2="$2"
      shift 2
      ;;
    -ref)
      if [[ -z "$2" ]]; then
        echo "ERROR: -ref requires a filename argument"
        exit 1
      fi
      REF="$2"
      shift 2
      ;;
    -*)
      echo "ERROR: Unknown option '$1'"
      echo "Valid options: -r1, -r2, -ref"
      exit 1
      ;;
    *)
      echo "ERROR: Unexpected argument '$1' - options must start with '-'"
      exit 1
      ;;
  esac
done

# Validate required arguments
if [[ -z "$R1" || -z "$R2" || -z "$REF" ]]; then
  echo "Usage: $0 -r1 <R1> -r2 <R2> -ref <REF>"
  echo "ERROR: Missing required arguments (-r1, -r2, or -ref)."
  exit 1
fi

echo "Processing MSA..."
mafft --auto "$REF" > "$OUTPUT_DIR/ref_seq_msa.aln"

echo "Building Bowtie2 index..."
bowtie2-build "$REF" "$OUTPUT_DIR/PKD_index"

echo "Aligning reads..."
bowtie2 --very-sensitive --end-to-end -x "$OUTPUT_DIR/PKD_index" -1 "$R1" -2 "$R2" -S "$OUTPUT_DIR/mapped_reads.sam"

echo "Mapping reference sequences to MSA..."
python "$SCRIPTS_DIR/ref_2_msa.py" --reference_fasta "$REF" --msa_file "$OUTPUT_DIR/ref_seq_msa.aln" --output "$OUTPUT_DIR/ref_seq_msa.tsv"

echo "Mapping reads to reference sequences..."
python "$SCRIPTS_DIR/read_2_gene.py" --sam "$OUTPUT_DIR/mapped_reads.sam" --output "$OUTPUT_DIR/mapped_reads.tsv"

echo "Mapping reads..."
python "$SCRIPTS_DIR/mapper_algo.py" --read_map "$OUTPUT_DIR/mapped_reads.tsv" --msa "$OUTPUT_DIR/ref_seq_msa.tsv" --output_file "$OUTPUT_DIR/unique_mappings.tsv"

echo "Writing files..."
python "$SCRIPTS_DIR/output.py" --tsv "$OUTPUT_DIR//unique_mappings.tsv" --r1 "$R1" --r2 "$R2" --ref "$REF" --fastq-dir "$OUTPUT_DIR/fastq" --bam-dir "$OUTPUT_DIR/bam"

echo "Pipeline complete. Outputs saved to: $OUTPUT_DIR"