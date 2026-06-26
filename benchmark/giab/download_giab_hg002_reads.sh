#!/usr/bin/env bash
set -euo pipefail

# Create directory for GIAB HG002 reads
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$SCRIPT_DIR/giab_hg002_reads"
cd "$SCRIPT_DIR/giab_hg002_reads"

# Base URL for GIAB HG002 reads
BASE_URL="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/"

# Download L001 lane files (001-017)
for i in {1..17}; do
    chunk=$(printf "%03d" $i)
    echo "Downloading L001 chunk ${chunk}..."
    wget -nc "${BASE_URL}D1_S1_L001_R1_${chunk}.fastq.gz"
    wget -nc "${BASE_URL}D1_S1_L001_R2_${chunk}.fastq.gz"
done

# Download L002 lane files (001-017)
for i in {1..17}; do
    chunk=$(printf "%03d" $i)
    echo "Downloading L002 chunk ${chunk}..."
    wget -nc "${BASE_URL}D1_S1_L002_R1_${chunk}.fastq.gz"
    wget -nc "${BASE_URL}D1_S1_L002_R2_${chunk}.fastq.gz"
done

echo "Creating merged FASTQs expected by run_giab.sh..."
cat D1_S1_L00*_R1_*.fastq.gz > HG002_R1.fq.gz
cat D1_S1_L00*_R2_*.fastq.gz > HG002_R2.fq.gz

echo "Done."
echo "  $SCRIPT_DIR/giab_hg002_reads/HG002_R1.fq.gz"
echo "  $SCRIPT_DIR/giab_hg002_reads/HG002_R2.fq.gz"
