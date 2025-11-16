#!/usr/bin/env bash
set -euo pipefail

mkdir -p gnaq_reads
cd gnaq_reads

echo "Downloading GNAQ-related reads into $(pwd)..."

########## SRR5602384 ##########
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/004/SRR5602384/SRR5602384_1.fastq.gz'
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/004/SRR5602384/SRR5602384_2.fastq.gz'

########## SRR5602389 ##########
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/009/SRR5602389/SRR5602389_1.fastq.gz'
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/009/SRR5602389/SRR5602389_2.fastq.gz'

########## SRR5602393 ##########
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/003/SRR5602393/SRR5602393_1.fastq.gz'
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/003/SRR5602393/SRR5602393_2.fastq.gz'

########## SRR5602414 ##########
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/004/SRR5602414/SRR5602414_1.fastq.gz'
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/004/SRR5602414/SRR5602414_2.fastq.gz'

########## SRR5602419 ##########
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/009/SRR5602419/SRR5602419_1.fastq.gz'
wget -c 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR560/009/SRR5602419/SRR5602419_2.fastq.gz'

echo "Download finished, now unzipping..."

gunzip -v *.fastq.gz

echo "Done. Uncompressed FASTQs are in: $(pwd)"
