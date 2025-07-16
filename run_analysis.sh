#!/bin/bash

# Create analysis directory if it doesn't exist
mkdir -p analysis

samtools view -b mapped_reads.sam > analysis/mapped_reads.bam

samtools sort -o analysis/sorted_mapped_reads.bam analysis/mapped_reads.bam 

samtools index analysis/sorted_mapped_reads.bam

freebayes -f ../../ref.fa -b analysis/sorted_mapped_reads.bam --min-alternate-count 5 > analysis/bowtie2.vcf

freebayes -f ../../ref.fa -b bam/PKD1.sorted.bam --min-alternate-count 5 > analysis/mapper1.vcf  

freebayes -f ../../ref.fa -b bam/PKD1P1.sorted.bam --min-alternate-count 5 > analysis/mapperp1.vcf

freebayes -f ../../ref.fa -b bam/PKD1P2.sorted.bam --min-alternate-count 5 > analysis/mapperp2.vcf    

freebayes -f ../../ref.fa -b bam/PKD1P3.sorted.bam --min-alternate-count 5 > analysis/mapperp3.vcf    

freebayes -f ../../ref.fa -b bam/PKD1P4.sorted.bam --min-alternate-count 5 > analysis/mapperp4.vcf    

freebayes -f ../../ref.fa -b bam/PKD1P5.sorted.bam --min-alternate-count 5 > analysis/mapperp5.vcf    

freebayes -f ../../ref.fa -b bam/PKD1P6.sorted.bam --min-alternate-count 5 > analysis/mapperp6.vcf

# Determine error rate from current directory
ERROR_RATE=$(basename $(pwd) | sed 's/err_//')

python ../../analysis/read_mapping_analysis.py -t unique_mappings.tsv -f ../../sim_reads/simulated_r1_err_${ERROR_RATE}.fq -s mapped_reads.sam -o analysis -p read_mapping_${ERROR_RATE}

python ../../analysis/snp_calling_analysis.py ../../sim_reads/all_genes_mutations.vcf analysis/mapper1.vcf analysis/mapperp1.vcf analysis/mapperp2.vcf analysis/mapperp3.vcf analysis/mapperp4.vcf analysis/mapperp5.vcf analysis/mapperp6.vcf --method2 analysis/bowtie2.vcf --tsv unique_mappings.tsv --sam mapped_reads.sam --fastq ../../sim_reads/simulated_r1_err_${ERROR_RATE}.fq --output-dir analysis/variant_calling_analysis

