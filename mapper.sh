#!/bin/bash

SCRIPTS_DIR="./scripts"
OUTPUT_DIR="./output_new"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "Starting pipeline..."

# Initialize variables
R1=""
R2=""
REF=""
FASTQ_OUT=0
BAM_OUT=0

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
    -out)
      if [[ -z "$2" ]]; then
        echo "ERROR: --out requires an output mode argument"
        exit 1
      fi
      IFS=',' read -ra OUT_MODES <<< "$2"
      for mode in "${OUT_MODES[@]}"; do
        case "$mode" in
          fastq)
            FASTQ_OUT=1
            ;;
          bam)
            BAM_OUT=1
            ;;
          *)
            echo "ERROR: Invalid output mode '$mode'. Valid modes are 'fastq' and 'bam'."
            exit 1
            ;;
        esac
      done
      shift 2
      ;;
    -*)
      echo "ERROR: Unknown option '$1'"
      echo "Valid options: -r1, -r2, -ref, -out"
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
  echo "Usage: $0 -r1 <R1> -r2 <R2> -ref <REF> [--out <fastq,bam>]"
  echo "ERROR: Missing required arguments (-r1, -r2, or -ref)."
  exit 1
fi

# Default to both outputs if none specified
if [[ $FASTQ_OUT -eq 0 && $BAM_OUT -eq 0 ]]; then
  FASTQ_OUT=1
  BAM_OUT=1
fi

echo "Processing MSA..."
mafft --auto "$REF" > "$OUTPUT_DIR/ref_seq_msa.aln"

echo "Building Bowtie2 index..."
bowtie2-build "$REF" "$OUTPUT_DIR/PKD_index"

echo "Aligning reads..."
bowtie2 --very-sensitive --end-to-end -x "$OUTPUT_DIR/PKD_index" -1 "$R1" -2 "$R2" -S "$OUTPUT_DIR/mapped_reads.sam"

echo "Mapping reference sequences to MSA..."
python3 "$SCRIPTS_DIR/ref_2_msa.py" --reference_fasta "$REF" --msa_file "$OUTPUT_DIR/ref_seq_msa.aln" --output "$OUTPUT_DIR/ref_seq_msa.tsv"

echo "Mapping reads to reference sequences..."
python3 "$SCRIPTS_DIR/read_2_gene.py" --sam "$OUTPUT_DIR/mapped_reads.sam" --output "$OUTPUT_DIR/mapped_reads.tsv"

echo "Mapping reads..."
python3 "$SCRIPTS_DIR/mapper_algo.py" --read_map "$OUTPUT_DIR/mapped_reads.tsv" --msa "$OUTPUT_DIR/ref_seq_msa.tsv" --output_file "$OUTPUT_DIR/unique_mappings.tsv"

# Conditional execution of output steps
if [[ $FASTQ_OUT -eq 1 ]]; then
  echo "Writing FASTQ files..."
  python3 "$SCRIPTS_DIR/fq_output.py" --tsv "$OUTPUT_DIR/unique_mappings.tsv" --r1 "$R1" --r2 "$R2" -o "$OUTPUT_DIR/unique_mappings_fastq_files"
fi

if [[ $BAM_OUT -eq 1 ]]; then
  echo "Writing BAM files..."
  python3 "$SCRIPTS_DIR/bam_output.py" --tsv "$OUTPUT_DIR/unique_mappings.tsv" --r1 "$R1" --r2 "$R2" --ref "$REF" --output "$OUTPUT_DIR/unique_mappings_bam_files"
fi

echo "Pipeline complete. Outputs saved to: $OUTPUT_DIR"
