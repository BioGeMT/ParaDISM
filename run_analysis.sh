#!/bin/bash

samtools view -b output/mapped_reads.sam > analysis/mapped_reads.bam

samtools sort -o analysis/sorted_mapped_reads.bam analysis/mapped_reads.bam 

samtools index analysis/sorted_mapped_reads.bam

freebayes -f ref.fa -b  analysis/sorted_mapped_reads.bam > analysis/bowtie2.vcf

freebayes -f ref.fa -b output/bam/PKD1.sorted.bam > analysis/mapper1.vcf  

freebayes -f ref.fa -b output/bam/PKD1P1.sorted.bam > analysis/mapperp1.vcf

freebayes -f ref.fa -b output/bam/PKD1P2.sorted.bam > analysis/mapperp2.vcf    

freebayes -f ref.fa -b output/bam/PKD1P3.sorted.bam > analysis/mapperp3.vcf    

freebayes -f ref.fa -b output/bam/PKD1P4.sorted.bam > analysis/mapperp4.vcf    

freebayes -f ref.fa -b output/bam/PKD1P5.sorted.bam > analysis/mapperp5.vcf    

freebayes -f ref.fa -b output/bam/PKD1P6.sorted.bam > analysis/mapperp6.vcf

python analysis/read_mapping_analysis.py -t output/unique_mappings.tsv -f simulated_r1.fq -s output/mapped_reads.sam -o analysis/read_mapping_analysis

python analysis/snp_calling_analysis.py all_genes_mutations.vcf analysis/mapper1.vcf analysis/mapperp1.vcf analysis/mapperp2.vcf analysis/mapperp3.vcf analysis/mapperp4.vcf analysis/mapperp5.vcf analysis/mapperp6.vcf --method2 analysis/bowtie2.vcf --output-dir analysis/variant_calling_analysis

