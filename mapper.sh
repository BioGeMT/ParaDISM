#!/bin/bash

# If no arguments provided, launch interactive mode
if [ $# -eq 0 ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    exec "$SCRIPT_DIR/mapper_interactive.sh"
fi

SCRIPTS_DIR="${SCRIPTS_DIR:-./scripts}"
OUTPUT_DIR="${OUTPUT_DIR:-./output}"
ALIGNER="${ALIGNER:-bowtie2}"
THREADS="${THREADS:-4}"
MINIMAP2_PROFILE="${MINIMAP2_PROFILE:-short}"
MINIMAP2_PROFILE_SET=0

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
    --read1)
      if [[ -z "$2" ]]; then
        echo "ERROR: --read1 requires a filename argument"
        exit 1
      fi
      R1="$2"
      shift 2
      ;;
    --read2)
      if [[ -z "$2" ]]; then
        echo "ERROR: --read2 requires a filename argument"
        exit 1
      fi
      R2="$2"
      shift 2
      ;;
    --reference)
      if [[ -z "$2" ]]; then
        echo "ERROR: --reference requires a filename argument"
        exit 1
      fi
      REF="$2"
      shift 2
      ;;
    --aligner)
      if [[ -z "$2" ]]; then
        echo "ERROR: --aligner requires an argument (bowtie2, bwa-mem2, or minimap2)"
        exit 1
      fi
      ALIGNER="$2"
      shift 2
      ;;
    --threads)
      if [[ -z "$2" ]]; then
        echo "ERROR: --threads requires an integer argument"
        exit 1
      fi
      THREADS="$2"
      shift 2
      ;;
    --minimap2-profile)
      if [[ -z "$2" ]]; then
        echo "ERROR: --minimap2-profile requires an argument (short, pacbio, or nanopore)"
        exit 1
      fi
      MINIMAP2_PROFILE="$2"
      MINIMAP2_PROFILE_SET=1
      shift 2
      ;;
    --*)
      echo "ERROR: Unknown option '$1'"
      echo "Valid options: --read1, --read2, --reference, --aligner, --threads, --minimap2-profile"
      exit 1
      ;;
    *)
      echo "ERROR: Unexpected argument '$1' - options must start with '--'"
      exit 1
      ;;
  esac
done

# Validate required arguments
if [[ -z "$R1" || -z "$R2" || -z "$REF" ]]; then
  echo "Usage: $0 --read1 <READ1> --read2 <READ2> --reference <REFERENCE> [--aligner <bowtie2|bwa-mem2|minimap2>] [--threads <N>] [--minimap2-profile <short|pacbio|nanopore>]"
  echo "ERROR: Missing required arguments (--read1, --read2, or --reference)."
  exit 1
fi

case "${ALIGNER,,}" in
  bowtie2|bwa-mem2|minimap2)
    ALIGNER="${ALIGNER,,}"
    ;;
  *)
    echo "ERROR: Unsupported aligner '$ALIGNER'. Choose from bowtie2, bwa-mem2, minimap2."
    exit 1
    ;;
esac

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "ERROR: --threads must be a positive integer."
  exit 1
fi

MINIMAP2_PROFILE="${MINIMAP2_PROFILE,,}"
if [[ "$ALIGNER" == "minimap2" ]]; then
  case "$MINIMAP2_PROFILE" in
    short|pacbio|nanopore)
      ;;
    *)
      echo "ERROR: Unsupported --minimap2-profile '$MINIMAP2_PROFILE'. Choose short, pacbio, or nanopore."
      exit 1
      ;;
  esac
elif [[ "$MINIMAP2_PROFILE_SET" -eq 1 ]]; then
  echo "WARNING: --minimap2-profile is ignored when --aligner is not minimap2."
fi

echo "Processing MSA..."
mafft --auto "$REF" > "$OUTPUT_DIR/ref_seq_msa.aln"

case "$ALIGNER" in
  bowtie2)
    echo "Building Bowtie2 index..."
    bowtie2-build "$REF" "$OUTPUT_DIR/PKD1_index"

    echo "Aligning reads with Bowtie2..."
    bowtie2 -p "$THREADS" -x "$OUTPUT_DIR/PKD1_index" -1 "$R1" -2 "$R2" -S "$OUTPUT_DIR/mapped_reads.sam"
    ;;
  bwa-mem2)
    echo "Building BWA-MEM2 index..."
    bwa-mem2 index -p "$OUTPUT_DIR/PKD1_index" "$REF"

    echo "Aligning reads with BWA-MEM2..."
    bwa-mem2 mem -t "$THREADS" "$OUTPUT_DIR/PKD1_index" "$R1" "$R2" > "$OUTPUT_DIR/mapped_reads.sam"
    ;;
  minimap2)
    case "$MINIMAP2_PROFILE" in
      short)
        MINIMAP2_INDEX_PRESET="sr"
        MINIMAP2_ALIGN_PRESET="sr"
        ;;
      pacbio)
        MINIMAP2_INDEX_PRESET="map-pb"
        MINIMAP2_ALIGN_PRESET="map-pb"
        ;;
      nanopore)
        MINIMAP2_INDEX_PRESET="map-ont"
        MINIMAP2_ALIGN_PRESET="map-ont"
        ;;
    esac

    echo "Building minimap2 ($MINIMAP2_PROFILE) index..."
    minimap2 -x "$MINIMAP2_INDEX_PRESET" -d "$OUTPUT_DIR/PKD1_index.mmi" "$REF"

    echo "Aligning reads with minimap2 ($MINIMAP2_PROFILE)..."
    minimap2 -ax "$MINIMAP2_ALIGN_PRESET" --MD -t "$THREADS" "$OUTPUT_DIR/PKD1_index.mmi" "$R1" "$R2" > "$OUTPUT_DIR/mapped_reads.sam"
    ;;
esac

echo "Mapping reference sequences to MSA..."
python "$SCRIPTS_DIR/ref_2_msa.py" --reference_fasta "$REF" --msa_file "$OUTPUT_DIR/ref_seq_msa.aln" --output "$OUTPUT_DIR/ref_seq_msa.tsv"

echo "Mapping reads to reference sequences..."
python "$SCRIPTS_DIR/read_2_gene.py" --sam "$OUTPUT_DIR/mapped_reads.sam" --output "$OUTPUT_DIR/mapped_reads.tsv"

echo "Refining reads..."
python "$SCRIPTS_DIR/mapper_algo_snp_only.py" --read_map "$OUTPUT_DIR/mapped_reads.tsv" --msa "$OUTPUT_DIR/ref_seq_msa.tsv" --output_file "$OUTPUT_DIR/unique_mappings.tsv"

echo "Writing output files..."
python "$SCRIPTS_DIR/output.py" --tsv "$OUTPUT_DIR/unique_mappings.tsv" --r1 "$R1" --r2 "$R2" --ref "$REF" --fastq-dir "$OUTPUT_DIR/fastq" --bam-dir "$OUTPUT_DIR/bam"

echo "Pipeline complete. Outputs saved to: $OUTPUT_DIR"
