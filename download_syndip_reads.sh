#!/bin/bash

# Create directory for SynDip reads
mkdir -p syndip_reads
cd syndip_reads

# Download SynDip reads
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/005/ERR1341795/ERR1341795_1.fastq.gz
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/006/ERR1341796/ERR1341796_2.fastq.gz
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/006/ERR1341796/ERR1341796_1.fastq.gz
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/003/ERR1341793/ERR1341793_2.fastq.gz
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/003/ERR1341793/ERR1341793_1.fastq.gz
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/004/ERR1341794/ERR1341794_2.fastq.gz
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/005/ERR1341795/ERR1341795_2.fastq.gz
wget -nc https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR134/004/ERR1341794/ERR1341794_1.fastq.gz

echo "Downloads complete! Unzipping files..."
gunzip -f *.fastq.gz

echo "Done! Unzipped FASTQ files are in syndip_reads/"
