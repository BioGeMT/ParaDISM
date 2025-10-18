#!/bin/bash

set -e  # Exit on error
trap 'tput cnorm; kill $(jobs -p) 2>/dev/null' EXIT  # Restore cursor and kill all background jobs on exit

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

# Create log file with timestamp
LOG_FILE="$OUTPUT_DIR/pipeline_$(date '+%Y%m%d_%H%M%S').log"

# Spinner characters
SPINNER_CHARS='⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'

# Colors
CYAN='\033[0;36m'
RED='\033[0;31m'
RESET='\033[0m'

# Function to show spinner while command runs
spinner() {
    local pid=$1
    local message=$2
    local spin_i=0

    # Hide cursor
    tput civis

    # Print message with spinner
    while kill -0 $pid 2>/dev/null; do
        local char=${SPINNER_CHARS:spin_i++%${#SPINNER_CHARS}:1}
        printf "\r%s %s" "$char" "$message"
        sleep 0.1
    done

    wait $pid
    local exit_code=$?

    # Show cursor
    tput cnorm

    if [ $exit_code -eq 0 ]; then
        printf "\r${CYAN}✓${RESET} ${CYAN}%s${RESET}\n" "$message"
    else
        printf "\r${RED}✗${RESET} ${RED}%s${RESET}\n" "$message"
        return $exit_code
    fi
}

# Function to run command with spinner
run_with_spinner() {
    local message=$1
    shift

    # Append to pipeline log with timestamp and separator
    echo "" >> "$LOG_FILE"
    echo "================================================" >> "$LOG_FILE"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" >> "$LOG_FILE"
    echo "================================================" >> "$LOG_FILE"

    # Run command in background, redirect output to log
    # Use eval to properly handle shell redirections in the command
    eval "$@" >> "$LOG_FILE" 2>&1 &
    local pid=$!

    # Show spinner
    if ! spinner $pid "$message"; then
        echo "Error occurred. Check log: $LOG_FILE"
        echo "Last 20 lines:"
        tail -20 "$LOG_FILE"
        return 1
    fi

    return 0
}

# Function to show progress output with spinner and completion message
run_with_progress() {
    local message=$1
    shift

    {
        echo ""
        echo "================================================"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message"
        echo "================================================"
    } >> "$LOG_FILE"

    tput civis

    # Print initial two lines: spinner + blank line for counter/progress
    printf "%s %s\n" "⠋" "$message"
    echo ""

    # Start spinner animation in background
    local spin_i=0
    (
        while true; do
            tput cuu 2  # Move cursor up 2 lines (to spinner line)
            tput cr     # Move to beginning of line
            local char=${SPINNER_CHARS:spin_i++%${#SPINNER_CHARS}:1}
            printf "%s %s" "$char" "$message"
            tput el     # Clear to end of line
            tput cud 2  # Move cursor down 2 lines (back to progress/output line)
            tput cr     # Move to beginning of line
            sleep 0.1
        done
    ) &
    local spinner_pid=$!

    set +e
    # Redirect stderr to log only, stdout goes to tee (screen + log)
    "$@" 2>>"$LOG_FILE" | tee -a "$LOG_FILE"
    local exit_code=${PIPESTATUS[0]}

    # Stop spinner
    kill -9 $spinner_pid 2>/dev/null
    wait $spinner_pid 2>/dev/null

    # Python script ends with newline, so we're now 1 line below counter line
    # Need to move up 3 lines total to reach spinner line
    tput cuu 3           # Move up 3 lines to spinner line
    tput cr              # Go to beginning of line
    if [ $exit_code -eq 0 ]; then
        printf "${CYAN}✓${RESET} ${CYAN}%s${RESET}\n" "$message"
    else
        printf "${RED}✗${RESET} ${RED}%s${RESET}\n" "$message"
    fi
    # Now cursor is on counter line after the newline from printf
    tput el              # Clear the counter/progress line
    tput cud 1           # Move down to blank line
    tput cr              # Go to beginning
    tput el              # Clear the blank line
    tput cuu 1           # Move back up to be ready for next command
    tput cr              # Go to beginning

    if [ $exit_code -ne 0 ]; then
        echo "Error occurred. Check log: $LOG_FILE"
        tail -20 "$LOG_FILE"
    fi

    # Show cursor again
    tput cnorm

    set -e

    if [ $exit_code -eq 0 ]; then
        return 0
    else
        return $exit_code
    fi
}

# Initialize variables
R1=""
R2=""
REF=""
SAM=""

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
    --sam)
      if [[ -z "$2" ]]; then
        echo "ERROR: --sam requires a filename argument"
        exit 1
      fi
      SAM="$2"
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
    --output-dir)
      if [[ -z "$2" ]]; then
        echo "ERROR: --output-dir requires a directory path argument"
        exit 1
      fi
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --*)
      echo "ERROR: Unknown option '$1'"
      echo "Valid options: --read1, --read2, --reference, --sam, --aligner, --threads, --minimap2-profile, --output-dir"
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
  echo "Usage: $0 --read1 <READ1> --read2 <READ2> --reference <REFERENCE> [--sam <SAM_FILE>] [--aligner <bowtie2|bwa-mem2|minimap2>] [--threads <N>] [--minimap2-profile <short|pacbio|nanopore>] [--output-dir <DIR>]"
  echo "ERROR: Missing required arguments (--read1, --read2, or --reference)."
  exit 1
fi

# Check if SAM file is provided
if [[ -n "$SAM" ]]; then
  if [[ ! -f "$SAM" ]]; then
    echo "ERROR: SAM file not found: $SAM"
    exit 1
  fi

  echo "Validating SAM file..."

  # Check if SAM file has valid header
  if ! samtools view -H "$SAM" &>/dev/null; then
    echo "ERROR: Invalid SAM file or missing header"
    exit 1
  fi

  # Check if SAM file has aligned reads
  ALIGNED_COUNT=$(samtools view -c -F 4 "$SAM" 2>/dev/null || echo "0")
  if [[ "$ALIGNED_COUNT" -eq 0 ]]; then
    echo "ERROR: SAM file contains no aligned reads"
    exit 1
  fi

  # Check if SAM file has MD tag (required for get_reference_sequence())
  HAS_MD=$(samtools view "$SAM" 2>/dev/null | head -n 100 | grep -c "MD:Z:" || echo "0")
  if [[ "$HAS_MD" -eq 0 ]]; then
    echo "ERROR: SAM file is missing required MD tag"
    exit 1
  fi

  echo "SAM validation passed ($ALIGNED_COUNT aligned reads)"
  echo "Using provided SAM file: $SAM"
  echo "Note: --aligner, --threads, and --minimap2-profile options will be ignored."
fi

# Validate aligner options only if SAM is not provided
if [[ -z "$SAM" ]]; then
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
fi

# Print pipeline header only if not called from interactive script
if [[ -z "$MAPPER_CALLED_FROM_INTERACTIVE" ]]; then
  echo -e "${CYAN}Running mapper pipeline...${RESET}"
  echo ""
fi

run_with_spinner "Processing MSA" "mafft --auto '$REF' > '$OUTPUT_DIR/ref_seq_msa.aln'"

# Perform alignment or use provided SAM file
if [[ -n "$SAM" ]]; then
  run_with_spinner "Copying provided SAM file" "cp '$SAM' '$OUTPUT_DIR/mapped_reads.sam'"
else
  case "$ALIGNER" in
    bowtie2)
      run_with_spinner "Building Bowtie2 index" "bowtie2-build '$REF' '$OUTPUT_DIR/ref_index'"
      run_with_spinner "Aligning reads with Bowtie2" "bowtie2 -p $THREADS -x '$OUTPUT_DIR/ref_index' -1 '$R1' -2 '$R2' -S '$OUTPUT_DIR/mapped_reads.sam'"
      ;;
    bwa-mem2)
      run_with_spinner "Building BWA-MEM2 index" "bwa-mem2 index -p '$OUTPUT_DIR/ref_index' '$REF'"
      run_with_spinner "Aligning reads with BWA-MEM2" "bwa-mem2 mem -t $THREADS '$OUTPUT_DIR/ref_index' '$R1' '$R2' > '$OUTPUT_DIR/mapped_reads.sam'"
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

      run_with_spinner "Building minimap2 ($MINIMAP2_PROFILE) index" "minimap2 -x $MINIMAP2_INDEX_PRESET -d '$OUTPUT_DIR/ref_index.mmi' '$REF'"
      run_with_spinner "Aligning reads with minimap2 ($MINIMAP2_PROFILE)" "minimap2 -ax $MINIMAP2_ALIGN_PRESET --MD -t $THREADS '$OUTPUT_DIR/ref_index.mmi' '$R1' '$R2' > '$OUTPUT_DIR/mapped_reads.sam'"
      ;;
  esac
fi

run_with_spinner "Mapping reference sequences to MSA" "python '$SCRIPTS_DIR/ref_2_msa.py' --reference_fasta '$REF' --msa_file '$OUTPUT_DIR/ref_seq_msa.aln' --output '$OUTPUT_DIR/ref_seq_msa.tsv'"

run_with_progress "Mapping reads to reference sequences" python "$SCRIPTS_DIR/read_2_gene.py" --sam "$OUTPUT_DIR/mapped_reads.sam" --output "$OUTPUT_DIR/mapped_reads.tsv" --fastq "$R1"

run_with_progress "Refining mapping" python "$SCRIPTS_DIR/mapper_algo_snp_only.py" --read_map "$OUTPUT_DIR/mapped_reads.tsv" --msa "$OUTPUT_DIR/ref_seq_msa.tsv" --output_file "$OUTPUT_DIR/unique_mappings.tsv"

run_with_progress "Writing output files" python "$SCRIPTS_DIR/output.py" --tsv "$OUTPUT_DIR/unique_mappings.tsv" --r1 "$R1" --r2 "$R2" --ref "$REF" --fastq-dir "$OUTPUT_DIR/fastq" --bam-dir "$OUTPUT_DIR/bam" --aligner "$ALIGNER" --threads "$THREADS" --minimap2-profile "$MINIMAP2_PROFILE"

# Brief delay to let output buffers flush
sleep 0.2

# Cleanup intermediate files
echo "Cleaning up intermediate files..."
rm -f "$OUTPUT_DIR"/ref_index.*
rm -f "$OUTPUT_DIR"/ref_seq_msa.aln
rm -f "$OUTPUT_DIR"/ref_seq_msa.tsv
rm -f "$OUTPUT_DIR"/mapped_reads.tsv

echo "Pipeline complete. Outputs saved to: $OUTPUT_DIR"
