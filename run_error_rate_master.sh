#!/bin/bash

# Master script to run mapper and analysis pipeline on all error rates
echo "Starting error rate analysis master script..."

# Define error rates
ERROR_RATES=(0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.010)

# Create results directory
mkdir -p error_rate_results

# Pre-compute MSA and reference mapping (only once)
echo "=========================================="
echo "Pre-computing MSA and reference mapping..."
echo "=========================================="

echo "Processing MSA..."
mafft --auto ref.fa > error_rate_results/ref_seq_msa.aln

echo "Mapping reference sequences to MSA..."
python scripts/ref_2_msa.py --reference_fasta ref.fa --msa_file error_rate_results/ref_seq_msa.aln --output error_rate_results/ref_seq_msa.tsv

echo "✓ MSA processing complete - will reuse for all error rates"

# Function to run complete pipeline for one error rate
run_error_rate() {
    local rate=$1
    local rate_underscore=${rate//./_}
    
    echo "=========================================="
    echo "Processing error rate: $rate"
    echo "=========================================="
    
    # Input files
    R1_FILE="simulated_r1_err_${rate_underscore}.fq"
    R2_FILE="simulated_r2_err_${rate_underscore}.fq"
    
    # Check if input files exist
    if [[ ! -f "$R1_FILE" || ! -f "$R2_FILE" ]]; then
        echo "ERROR: Input files not found for error rate $rate"
        return 1
    fi
    
    # Create output directory for this error rate
    OUTPUT_DIR="error_rate_results/err_${rate_underscore}"
    ANALYSIS_DIR="$OUTPUT_DIR/analysis"
    mkdir -p "$OUTPUT_DIR" "$ANALYSIS_DIR"
    
    echo "Running mapper pipeline for error rate $rate..."
    
    # Step 1: Run mapper.sh with modified output directory
    SCRIPTS_DIR="./scripts"
    
    # Create output directory if it doesn't exist
    mkdir -p "$OUTPUT_DIR"
    
    # Copy pre-computed MSA files (no need to regenerate)
    cp error_rate_results/ref_seq_msa.aln "$OUTPUT_DIR/"
    cp error_rate_results/ref_seq_msa.tsv "$OUTPUT_DIR/"
    
    echo "Building Bowtie2 index..."
    bowtie2-build ref.fa "$OUTPUT_DIR/PKD1_index"
    
    echo "Aligning reads..."
    bowtie2 -x "$OUTPUT_DIR/PKD1_index" -1 "$R1_FILE" -2 "$R2_FILE" -S "$OUTPUT_DIR/mapped_reads.sam"
    
    echo "Mapping reads to reference sequences..."
    python "$SCRIPTS_DIR/read_2_gene.py" --sam "$OUTPUT_DIR/mapped_reads.sam" --output "$OUTPUT_DIR/mapped_reads.tsv"
    
    echo "Mapping reads..."
    python "$SCRIPTS_DIR/mapper_algo.py" --read_map "$OUTPUT_DIR/mapped_reads.tsv" --msa "$OUTPUT_DIR/ref_seq_msa.tsv" --output_file "$OUTPUT_DIR/unique_mappings.tsv"
    
    echo "Writing files..."
    python "$SCRIPTS_DIR/output.py" --tsv "$OUTPUT_DIR/unique_mappings.tsv" --r1 "$R1_FILE" --r2 "$R2_FILE" --ref ref.fa --fastq-dir "$OUTPUT_DIR/fastq" --bam-dir "$OUTPUT_DIR/bam"
    
    # Step 2: Run analysis pipeline
    echo "Running analysis for error rate $rate..."
    
    # Convert SAM to BAM and sort
    samtools view -b "$OUTPUT_DIR/mapped_reads.sam" > "$ANALYSIS_DIR/mapped_reads.bam"
    samtools sort -o "$ANALYSIS_DIR/sorted_mapped_reads.bam" "$ANALYSIS_DIR/mapped_reads.bam"
    samtools index "$ANALYSIS_DIR/sorted_mapped_reads.bam"
    
    # Run FreeBayes on bowtie2 results
    freebayes -f ref.fa -b "$ANALYSIS_DIR/sorted_mapped_reads.bam" > "$ANALYSIS_DIR/bowtie2.vcf"
    
    # Run FreeBayes on gene-specific BAM files
    if [[ -d "$OUTPUT_DIR/bam" ]]; then
        for bam_file in "$OUTPUT_DIR/bam"/*.sorted.bam; do
            if [[ -f "$bam_file" ]]; then
                gene_name=$(basename "$bam_file" .sorted.bam)
                freebayes -f ref.fa -b "$bam_file" > "$ANALYSIS_DIR/mapper_${gene_name}.vcf"
            fi
        done
    fi
    
    # Run read mapping analysis
    python analysis/read_mapping_analysis.py \
        -t "$OUTPUT_DIR/unique_mappings.tsv" \
        -f "$R1_FILE" \
        -s "$OUTPUT_DIR/mapped_reads.sam" \
        -o "$ANALYSIS_DIR" \
        -p "read_mapping_err_${rate_underscore}"
    
    # Run SNP calling analysis
    mapper_vcfs=()
    for vcf in "$ANALYSIS_DIR"/mapper_*.vcf; do
        if [[ -f "$vcf" ]]; then
            mapper_vcfs+=("$vcf")
        fi
    done
    
    if [[ ${#mapper_vcfs[@]} -gt 0 ]]; then
        python analysis/snp_calling_analysis.py \
            all_genes_mutations.vcf \
            "${mapper_vcfs[@]}" \
            --method2 "$ANALYSIS_DIR/bowtie2.vcf" \
            --tsv "$OUTPUT_DIR/unique_mappings.tsv" \
            --sam "$OUTPUT_DIR/mapped_reads.sam" \
            --fastq "$R1_FILE" \
            --output-dir "$ANALYSIS_DIR/variant_calling_analysis"
    fi
    
    echo "✓ Completed processing error rate $rate"
    echo ""
}

# Run pipeline for each error rate
for rate in "${ERROR_RATES[@]}"; do
    run_error_rate "$rate"
done

echo "=========================================="
echo "Error rate analysis completed!"
echo "Results saved in: error_rate_results/"
echo "=========================================="

# Create summary
echo "Creating summary of results..."
echo "Error_Rate,Mapper_Precision,Mapper_Recall,Mapper_F1,Bowtie2_Precision,Bowtie2_Recall,Bowtie2_F1" > error_rate_results/summary.csv

for rate in "${ERROR_RATES[@]}"; do
    rate_underscore=${rate//./_}
    metrics_file="error_rate_results/err_${rate_underscore}/analysis/read_mapping_err_${rate_underscore}_comprehensive_metrics.tsv"
    
    if [[ -f "$metrics_file" ]]; then
        # Extract aggregated metrics from the comprehensive file
        tail -n 2 "$metrics_file" | head -n 1 | awk -F'\t' -v rate="$rate" '{
            printf "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", 
            rate, $5, $6, $7, $12, $13, $14
        }' >> error_rate_results/summary.csv
    fi
done

echo "Summary created: error_rate_results/summary.csv"

# Generate plots
echo "Generating performance plots..."
python plot_error_rate_analysis.py --results-dir error_rate_results

echo "✓ Analysis complete! Check error_rate_results/ for all results and plots."