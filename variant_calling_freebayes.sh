#!/bin/bash

# FreeBayes variant calling on mapper outputs
OUTPUT_DIR="${1:-./output}"
REF_FILE="${2:-ref.fa}"

cd "$OUTPUT_DIR"

# Create analysis directory and convert SAM to BAM
mkdir -p analysis
samtools view -b mapped_reads.sam > analysis/mapped_reads.bam
samtools sort -o analysis/sorted_mapped_reads.bam analysis/mapped_reads.bam
samtools index analysis/sorted_mapped_reads.bam

# FreeBayes variant calling with enhanced QC + normalization
# Note: Some FreeBayes builds do not support --no-indels/--no-mnps/--no-complex.
# We keep SNP-only downstream via 'bcftools view -v snps'.
# Enforce thresholds in FreeBayes itself and restrict to SNPs only (use short flags for compatibility)
FB_PARAMS="-q 13 -m 20 -C 5 --min-coverage 8 -F 0.2 -i -X -u"

# Ensure the reference FASTA has a .fai index (works for consensus per-iteration)
samtools faidx "$REF_FILE"

#### FreeBayes + bcftools normalization and filtering (SNP-only)
# Bowtie2 aggregate BAM
freebayes -f "$REF_FILE" -b analysis/sorted_mapped_reads.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/bowtie2.freebayes.norm.vcf.gz
# Output normalized VCF (SNP-only already enforced by FreeBayes flags)
gunzip -c analysis/bowtie2.freebayes.norm.vcf.gz > analysis/bowtie2_freebayes_norm.vcf

## FreeBayes + normalization + SNP filtering + thresholding for per-gene BAMs
freebayes -f "$REF_FILE" -b bam/PKD1.sorted.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/mapper1.freebayes.norm.vcf.gz
gunzip -c analysis/mapper1.freebayes.norm.vcf.gz > analysis/mapper1_freebayes_norm.vcf

freebayes -f "$REF_FILE" -b bam/PKD1P1.sorted.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/mapperp1.freebayes.norm.vcf.gz
gunzip -c analysis/mapperp1.freebayes.norm.vcf.gz > analysis/mapperp1_freebayes_norm.vcf

freebayes -f "$REF_FILE" -b bam/PKD1P2.sorted.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/mapperp2.freebayes.norm.vcf.gz
gunzip -c analysis/mapperp2.freebayes.norm.vcf.gz > analysis/mapperp2_freebayes_norm.vcf

freebayes -f "$REF_FILE" -b bam/PKD1P3.sorted.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/mapperp3.freebayes.norm.vcf.gz
gunzip -c analysis/mapperp3.freebayes.norm.vcf.gz > analysis/mapperp3_freebayes_norm.vcf

freebayes -f "$REF_FILE" -b bam/PKD1P4.sorted.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/mapperp4.freebayes.norm.vcf.gz
gunzip -c analysis/mapperp4.freebayes.norm.vcf.gz > analysis/mapperp4_freebayes_norm.vcf

freebayes -f "$REF_FILE" -b bam/PKD1P5.sorted.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/mapperp5.freebayes.norm.vcf.gz
gunzip -c analysis/mapperp5.freebayes.norm.vcf.gz > analysis/mapperp5_freebayes_norm.vcf

freebayes -f "$REF_FILE" -b bam/PKD1P6.sorted.bam $FB_PARAMS | \
bcftools norm -f "$REF_FILE" -Oz -o analysis/mapperp6.freebayes.norm.vcf.gz
gunzip -c analysis/mapperp6.freebayes.norm.vcf.gz > analysis/mapperp6_freebayes_norm.vcf

# Copy final VCF files to variant_calling directory for consensus building
mkdir -p variant_calling
cp analysis/bowtie2_freebayes_norm.vcf variant_calling/
cp analysis/mapper1_freebayes_norm.vcf variant_calling/
cp analysis/mapperp1_freebayes_norm.vcf variant_calling/
cp analysis/mapperp2_freebayes_norm.vcf variant_calling/
cp analysis/mapperp3_freebayes_norm.vcf variant_calling/
cp analysis/mapperp4_freebayes_norm.vcf variant_calling/
cp analysis/mapperp5_freebayes_norm.vcf variant_calling/
cp analysis/mapperp6_freebayes_norm.vcf variant_calling/
