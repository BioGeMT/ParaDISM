#!/bin/bash

# Create directory for GIAB HG002 reads
mkdir -p giab_hg002_reads
cd giab_hg002_reads

# Base URL for GIAB HG002 reads
BASE_URL="https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/reads/"

echo "Downloading GIAB HG002 reads..."
echo "This will download 34 files (17 pairs from L001 + 17 pairs from L002)"
echo ""

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

echo ""
echo "Downloads complete! Unzipping files..."
gunzip -f *.fastq.gz

echo ""
echo "Concatenating files by lane..."
# Concatenate L001 lane files
cat D1_S1_L001_R1_*.fastq > L001_R1.fastq
cat D1_S1_L001_R2_*.fastq > L001_R2.fastq
# Concatenate L002 lane files
cat D1_S1_L002_R1_*.fastq > L002_R1.fastq
cat D1_S1_L002_R2_*.fastq > L002_R2.fastq

echo ""
echo "Combining all lanes into single files..."
# Combine all lanes
cat L001_R1.fastq L002_R1.fastq > HG002_R1.fastq
cat L001_R2.fastq L002_R2.fastq > HG002_R2.fastq

echo ""
echo "Done! Final files:"
echo "  Lane-specific (4 files):"
echo "    L001_R1.fastq"
echo "    L001_R2.fastq"
echo "    L002_R1.fastq"
echo "    L002_R2.fastq"
echo ""
echo "  Combined (2 files):"
echo "    HG002_R1.fastq"
echo "    HG002_R2.fastq"
echo ""
echo "Individual chunk files are still in giab_hg002_reads/ if needed."


