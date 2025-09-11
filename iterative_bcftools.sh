#!/bin/bash
set -euo pipefail

# Iterative bcftools variant calling with consensus reference improvement
OUTPUT_DIR="${1:-./output}"
REF_FILE="${2:-ref.fa}"
ITERATIONS="${3:-2}"

echo "=== Iterative bcftools Pipeline ==="
echo "Iterations: $ITERATIONS"
echo "Output directory: $OUTPUT_DIR"
echo "Initial reference: $REF_FILE"

mkdir -p "$OUTPUT_DIR"
# MAP_REF evolves for mapping only; CALL_REF stays as the original for calling
MAP_REF="$REF_FILE"
CALL_REF="$REF_FILE"

for ((i=1; i<=ITERATIONS; i++)); do
    echo ""
    echo "=== Iteration $i/$ITERATIONS ==="
    
    ITER_DIR="$OUTPUT_DIR/iter_$i"
    mkdir -p "$ITER_DIR/variant_calling"
    cd "$ITER_DIR"
    
    echo "Running mapper.sh with mapping ref: $MAP_REF"
    SCRIPTS_DIR="../../scripts" OUTPUT_DIR="." ../../mapper.sh -r1 "../../${4}" -r2 "../../${5}" -ref "../../$MAP_REF"

    # Rebuild per-gene BAMs against the original calling reference
    echo "Re-aligning per-gene FASTQs to calling ref: $CALL_REF"
    python3 ../../scripts/output.py \
        --tsv unique_mappings.tsv \
        --r1 "../../${4}" \
        --r2 "../../${5}" \
        --ref "../../$CALL_REF" \
        --fastq-dir "fastq" \
        --bam-dir "bam"
    
    echo "Running bcftools variant calling against original ref..."
    ../../variant_calling_bcftools.sh "." "../../$CALL_REF"
    
    if [[ $i -lt $ITERATIONS ]]; then
        echo "Creating consensus reference for next iteration..."
        
        # Extract individual gene references from the original calling reference
        python3 - << EOF
from Bio import SeqIO
import sys
sys.path.append('../../')

# Read the calling reference relative to the iteration dir
sequences = list(SeqIO.parse("../../$CALL_REF", "fasta"))
for seq in sequences:
    gene_name = seq.id
    output_file = f"variant_calling/{gene_name}_ref.fa"
    with open(output_file, "w") as f:
        SeqIO.write(seq, f, "fasta")
EOF

        # Create mapping consensus sequences from SNP-only normalized VCFs
        GENES=("PKD1" "PKD1P1" "PKD1P2" "PKD1P3" "PKD1P4" "PKD1P5" "PKD1P6")
        MAPPER_NAMES=("mapper1" "mapperp1" "mapperp2" "mapperp3" "mapperp4" "mapperp5" "mapperp6")
        
        for j in "${!GENES[@]}"; do
            gene="${GENES[$j]}"
            mapper_name="${MAPPER_NAMES[$j]}"
            
            if [[ -f "variant_calling/${mapper_name}_bcftools_norm.vcf" ]]; then
                # Keep SNPs only and compress for bcftools consensus
                bcftools view -v snps "variant_calling/${mapper_name}_bcftools_norm.vcf" | \
                    bgzip -c > "variant_calling/${mapper_name}.snps.vcf.gz"
                tabix -p vcf "variant_calling/${mapper_name}.snps.vcf.gz"
                
                # Create consensus
                bcftools consensus -f "variant_calling/${gene}_ref.fa" "variant_calling/${mapper_name}.snps.vcf.gz" > "variant_calling/${gene}_consensus.fa"
                sed -i "s/^>.*/>$gene/" "variant_calling/${gene}_consensus.fa"
            fi
        done
        
        # Combine consensus sequences with strict guard (no fallback)
        shopt -s nullglob
        CONS_FILES=(variant_calling/*_consensus.fa)
        if (( ${#CONS_FILES[@]} == 0 )); then
            echo "Error: no per-gene consensus FASTAs were produced; cannot update mapping reference." >&2
            exit 1
        fi
        cat "${CONS_FILES[@]}" > "variant_calling/consensus_reference.fa"
        if [[ ! -s "variant_calling/consensus_reference.fa" ]]; then
            echo "Error: consensus_reference.fa is empty; cannot update mapping reference." >&2
            exit 1
        fi
        shopt -u nullglob
        MAP_REF="$OUTPUT_DIR/iter_$i/variant_calling/consensus_reference.fa"
        echo "Consensus reference created (mapping only): $MAP_REF"
        
        # Clean up intermediate files
        echo "Cleaning up intermediate files..."
        rm -f variant_calling/*.snps.vcf.gz variant_calling/*.snps.vcf.gz.tbi
        rm -f variant_calling/*_ref.fa variant_calling/*_consensus.fa
    fi
    
    # Clean up mapper and bcftools intermediate files
    echo "Cleaning up intermediate files..."
    rm -f PKD1_index.*.bt2
    rm -f analysis/mapped_reads.bam
    rm -f analysis/*.vcf.gz analysis/*.norm.vcf.gz analysis/*.norm.flt.vcf.gz
    
    cd ../..
done

echo ""
echo "=== Iterative bcftools Pipeline Complete ==="
echo "Results in: $OUTPUT_DIR/"
