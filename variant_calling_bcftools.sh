#!/bin/bash

# bcftools variant calling on mapper outputs
OUTPUT_DIR="${1:-./output}"
REF_FILE="${2:-ref.fa}"

cd "$OUTPUT_DIR"

# Create analysis directory and convert SAM to BAM
mkdir -p analysis
samtools view -b mapped_reads.sam > analysis/mapped_reads.bam
samtools sort -o analysis/sorted_mapped_reads.bam analysis/mapped_reads.bam
samtools index analysis/sorted_mapped_reads.bam

# Ensure the reference FASTA has a .fai index (works for consensus per-iteration)
samtools faidx "$REF_FILE"

# bcftools variant calling with aligned parameters (SNP-only)
# Annotate is kept minimal; downstream filters use INFO/DP4 which is widely available
bcftools mpileup -Ou -f "$REF_FILE" --min-BQ 13 --min-MQ 20 analysis/sorted_mapped_reads.bam | \
bcftools call -mv -Oz -o analysis/bowtie2.vcf.gz

bcftools norm -f "$REF_FILE" analysis/bowtie2.vcf.gz -Oz -o analysis/bowtie2.norm.vcf.gz
# Keep SNPs only
bcftools view -v snps analysis/bowtie2.norm.vcf.gz -Oz -o analysis/bowtie2.norm.snps.vcf.gz
# Apply thresholds to mirror FreeBayes (no QUAL cutoff):
#  - depth >=8: sum(DP4) >= 8
#  - alt count >=5: DP4[2]+DP4[3] >= 5
#  - VAF >=0.2: (DP4[2]+DP4[3]) / sum(DP4) >= 0.2
bcftools filter -e '((INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])<8) || ((INFO/DP4[2]+INFO/DP4[3])<5) || ((INFO/DP4[2]+INFO/DP4[3])/(INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])<0.2)' analysis/bowtie2.norm.snps.vcf.gz -Oz -o analysis/bowtie2.norm.snps.flt.vcf.gz
gunzip -c analysis/bowtie2.norm.snps.flt.vcf.gz > analysis/bowtie2_bcftools_norm.vcf

# bcftools variant calling - Mapper results  
GENES=("PKD1" "PKD1P1" "PKD1P2" "PKD1P3" "PKD1P4" "PKD1P5" "PKD1P6")
MAPPER_NAMES=("mapper1" "mapperp1" "mapperp2" "mapperp3" "mapperp4" "mapperp5" "mapperp6")

for i in "${!GENES[@]}"; do
    gene="${GENES[$i]}"
    mapper_name="${MAPPER_NAMES[$i]}"
    
    # Call variants; normalize; keep SNPs; apply equivalent filters using INFO/DP4
    bcftools mpileup -Ou -f "$REF_FILE" --min-BQ 13 --min-MQ 20 bam/${gene}.sorted.bam | \
    bcftools call -mv -Oz -o analysis/${mapper_name}.vcf.gz

    bcftools norm -f "$REF_FILE" analysis/${mapper_name}.vcf.gz -Oz -o analysis/${mapper_name}.norm.vcf.gz
    bcftools view -v snps analysis/${mapper_name}.norm.vcf.gz -Oz -o analysis/${mapper_name}.norm.snps.vcf.gz
    bcftools filter -e '((INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])<8) || ((INFO/DP4[2]+INFO/DP4[3])<5) || ((INFO/DP4[2]+INFO/DP4[3])/(INFO/DP4[0]+INFO/DP4[1]+INFO/DP4[2]+INFO/DP4[3])<0.2)' analysis/${mapper_name}.norm.snps.vcf.gz -Oz -o analysis/${mapper_name}.norm.snps.flt.vcf.gz
    gunzip -c analysis/${mapper_name}.norm.snps.flt.vcf.gz > analysis/${mapper_name}_bcftools_norm.vcf
done

# Copy final VCF files to variant_calling directory for consensus building
mkdir -p variant_calling
cp analysis/bowtie2_bcftools_norm.vcf variant_calling/
for i in "${!MAPPER_NAMES[@]}"; do
    mapper_name="${MAPPER_NAMES[$i]}"
    cp analysis/${mapper_name}_bcftools_norm.vcf variant_calling/
done
