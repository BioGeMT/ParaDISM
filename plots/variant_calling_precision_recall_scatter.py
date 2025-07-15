#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

COLORS = ['#FF8C00', '#008080']  

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 24,
    'axes.titlesize': 28,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 14
})

def load_variant_metrics(error_rate):
    rate_str = f"{error_rate:.3f}".replace('.', '_')
    file_path = f"error_rate_results/err_{rate_str}/analysis/variant_calling_analysis/comprehensive_variant_metrics.tsv"
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    return pd.read_csv(file_path, sep='\t')

def create_precision_recall_scatter(error_rate):
    df = load_variant_metrics(error_rate)
    
    if df is None:
        print(f"Could not load variant metrics for error rate {error_rate}")
        return
    
    print(f"\n=== Processing Variant Calling Precision vs Recall Scatter - Error Rate {error_rate} ===")
    
    gene_rows = df[df['Gene'] != 'AGGREGATED'].copy()
    
    genes = gene_rows['Gene'].tolist()
    mapper_precision = gene_rows['Method1_Precision'].tolist()
    mapper_recall = gene_rows['Method1_Recall'].tolist()
    bowtie2_precision = gene_rows['Method2_Precision'].tolist()
    bowtie2_recall = gene_rows['Method2_Recall'].tolist()
    
    plt.figure(figsize=(10, 8), dpi=300)
    
    plt.scatter(mapper_recall, mapper_precision, c=COLORS[0], s=100, alpha=0.8, 
                label='Mapper', marker='o')
    plt.scatter(bowtie2_recall, bowtie2_precision, c=COLORS[1], s=100, alpha=0.8, 
                label='Bowtie2', marker='s')
    
    for i, gene in enumerate(genes):
        plt.annotate(gene, (mapper_recall[i], mapper_precision[i]), 
                    xytext=(8, 0), textcoords='offset points', fontsize=10, 
                    color='black', fontweight='bold', ha='left', va='center')
        plt.annotate(gene, (bowtie2_recall[i], bowtie2_precision[i]), 
                    xytext=(8, 0), textcoords='offset points', fontsize=10, 
                    color='black', fontweight='bold', ha='left', va='center')
    
    plt.xlabel('Recall', fontweight='bold')
    plt.ylabel('Precision', fontweight='bold')
    plt.title(f'Variant Calling Precision vs Recall - Error Rate {error_rate*100:.1f}%', fontweight='bold')
    
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    
    plt.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0))
    
    plt.tight_layout()
    
    filename = f"variant_calling_precision_recall_scatter_{f'{error_rate:.3f}'.replace('.', '_')}.png"
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"Created scatter plot: {filename}")
    
    print(f"\n=== PRECISION vs RECALL DATA FOR ERROR RATE {error_rate} ===")
    print("Gene         | Mapper P/R        | Bowtie2 P/R")
    print("-" * 50)
    for i, gene in enumerate(genes):
        print(f"{gene:<12} | {mapper_precision[i]:.3f}/{mapper_recall[i]:.3f} | {bowtie2_precision[i]:.3f}/{bowtie2_recall[i]:.3f}")

def main():
    error_rates = [0.001, 0.005, 0.010]
    for rate in error_rates:
        create_precision_recall_scatter(rate)

if __name__ == "__main__":
    main()