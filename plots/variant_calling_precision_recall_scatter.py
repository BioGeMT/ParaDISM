#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.metrics import precision_score, recall_score

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
    file_path = f"error_rate_results/err_{rate_str}/analysis/variant_calling_analysis/variant_calling_confusion_matrix.tsv"
    
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
    
    # Create y_true and y_pred arrays for sklearn
    y_true_method1 = []
    y_pred_method1 = []
    y_true_method2 = []
    y_pred_method2 = []
    
    for _, row in df.iterrows():
        gene = row['Gene']
        
        # Method 1 - Add TP as correct predictions
        y_true_method1.extend([gene] * int(row['Method1_TP']))
        y_pred_method1.extend([gene] * int(row['Method1_TP']))
        
        # Method 1 - Add FP as incorrect predictions
        y_true_method1.extend(['OTHER'] * int(row['Method1_FP']))
        y_pred_method1.extend([gene] * int(row['Method1_FP']))
        
        # Method 1 - Add FN as missed predictions
        y_true_method1.extend([gene] * int(row['Method1_FN']))
        y_pred_method1.extend(['OTHER'] * int(row['Method1_FN']))
        
        # Method 2 - Same logic
        y_true_method2.extend([gene] * int(row['Method2_TP']))
        y_pred_method2.extend([gene] * int(row['Method2_TP']))
        
        y_true_method2.extend(['OTHER'] * int(row['Method2_FP']))
        y_pred_method2.extend([gene] * int(row['Method2_FP']))
        
        y_true_method2.extend([gene] * int(row['Method2_FN']))
        y_pred_method2.extend(['OTHER'] * int(row['Method2_FN']))
    
    genes = df['Gene'].unique().tolist()
    
    mapper_precision = []
    mapper_recall = []
    bowtie2_precision = []
    bowtie2_recall = []
    
    # Calculate per-gene metrics
    for gene in genes:
        # Create binary classification for this gene
        y_true_gene_m1 = [1 if x == gene else 0 for x in y_true_method1]
        y_pred_gene_m1 = [1 if x == gene else 0 for x in y_pred_method1]
        
        y_true_gene_m2 = [1 if x == gene else 0 for x in y_true_method2]
        y_pred_gene_m2 = [1 if x == gene else 0 for x in y_pred_method2]
        
        precision1 = precision_score(y_true_gene_m1, y_pred_gene_m1, zero_division=0)
        recall1 = recall_score(y_true_gene_m1, y_pred_gene_m1, zero_division=0)
        
        precision2 = precision_score(y_true_gene_m2, y_pred_gene_m2, zero_division=0)
        recall2 = recall_score(y_true_gene_m2, y_pred_gene_m2, zero_division=0)
        
        mapper_precision.append(precision1)
        mapper_recall.append(recall1)
        bowtie2_precision.append(precision2)
        bowtie2_recall.append(recall2)
    
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
    
    filename = f"plots/variant_calling_precision_recall_scatter_{f'{error_rate:.3f}'.replace('.', '_')}.png"
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"Created scatter plot: {filename}")
    
    print(f"\n=== PRECISION vs RECALL DATA FOR ERROR RATE {error_rate} ===")
    print("Gene         | Mapper P/R        | Bowtie2 P/R")
    print("-" * 50)
    for i, gene in enumerate(genes):
        print(f"{gene:<12} | {mapper_precision[i]:.3f}/{mapper_recall[i]:.3f} | {bowtie2_precision[i]:.3f}/{bowtie2_recall[i]:.3f}")

def main():
    error_rates = [0.001, 0.005]
    for rate in error_rates:
        create_precision_recall_scatter(rate)

if __name__ == "__main__":
    main()