#!/bin/bash

# Plot iterative results using no-normalization variant confusion matrices
METHOD="${1:-freebayes}"  # freebayes or bcftools
ITERATIONS="${2:-3}"
OUTPUT_DIR="${3:-./output}"
R1_FILE="${4:-simulated_r1_err_0_001.fq}"
GROUND_TRUTH_VCF="${5:-all_genes_mutations.vcf}"

echo "=== Plotting Iterative Results (no-norm): $METHOD ==="
echo "Iterations: $ITERATIONS"
echo "Output directory: $OUTPUT_DIR"
echo "FASTQ (R1): $R1_FILE"
echo "Ground truth VCF: $GROUND_TRUTH_VCF"

PLOT_ROOT="plots_iterative_no_norm_${METHOD}"
mkdir -p "$PLOT_ROOT"

for ((i=1; i<=ITERATIONS; i++)); do
    echo ""
    echo "Processing iteration $i..."
    
    ITER_DIR="$OUTPUT_DIR/iter_$i"
    PLOT_DIR="$PLOT_ROOT/iter_$i"
    mkdir -p "$PLOT_DIR"
    
    if [[ -d "$ITER_DIR" ]]; then
        cd "$ITER_DIR"
        
        # Ensure SAM exists (recreate from BAM if needed)
        if [[ ! -f "mapped_reads.sam" && -f "analysis/sorted_mapped_reads.bam" ]]; then
            echo "Reconstructing mapped_reads.sam from BAM for iteration $i..."
            samtools view -h analysis/sorted_mapped_reads.bam > mapped_reads.sam
        fi

        # Run read mapping analysis if not already done (save to plot dir)
        if [[ ! -f "../../$PLOT_DIR/read_mapping_${METHOD}_iter_${i}_mapper_confusion_matrix.tsv" ]]; then
            echo "Running read mapping analysis for iteration $i..."
        python3 ../../analysis/read_mapping_analysis.py \
            -t unique_mappings.tsv \
            -f "../../$R1_FILE" \
            -s mapped_reads.sam \
            -o "../../$PLOT_DIR" \
            -p "read_mapping_${METHOD}_iter_${i}"
        fi
        
        # No variant-calling analysis here (no-norm mode focuses on read mapping only)
        
        # Generate plots (reads from plot dir, writes to plot dir)
        echo "Generating plots for iteration $i..."
        python3 ../../analysis/plot_iterative_analysis.py \
            --iteration $i \
            --method $METHOD \
            --iter-dir "../../$PLOT_DIR" \
            --output-dir "../../$PLOT_DIR"
            
        cd ../..
    else
        echo "Warning: Directory $ITER_DIR not found, skipping iteration $i"
    fi
done

echo ""
echo "Generating summary plots across iterations (no-norm)..."
python3 analysis/plot_iterative_summary.py \
    --method $METHOD \
    --iterations $ITERATIONS \
    --output-dir "$PLOT_ROOT" \
    --plot-dir "$PLOT_ROOT"

echo ""
echo "=== Plotting Complete (no-norm) ==="
echo "Results in: ${PLOT_ROOT}/"
echo "Summary plots:"
echo "  - read_mapping_summary_${METHOD}.png"
echo "  - variant_calling_summary_${METHOD}.png"
