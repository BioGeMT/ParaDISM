#!/bin/bash

# Create directory for GIAB HG002 reads
mkdir -p giab_hg002_reads
cd giab_hg002_reads

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