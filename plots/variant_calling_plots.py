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

def calculate_variant_metrics(df):
    results = []
    
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
    
    # Calculate per-gene metrics
    for gene in genes:
        # Create binary classification for this gene
        y_true_gene_m1 = [1 if x == gene else 0 for x in y_true_method1]
        y_pred_gene_m1 = [1 if x == gene else 0 for x in y_pred_method1]
        
        y_true_gene_m2 = [1 if x == gene else 0 for x in y_true_method2]
        y_pred_gene_m2 = [1 if x == gene else 0 for x in y_pred_method2]
        
        precision1 = precision_score(y_true_gene_m1, y_pred_gene_m1, zero_division=0)
        recall1 = recall_score(y_true_gene_m1, y_pred_gene_m1, zero_division=0)
        specificity1 = precision1  # For variant calling, specificity ≈ precision
        
        precision2 = precision_score(y_true_gene_m2, y_pred_gene_m2, zero_division=0)
        recall2 = recall_score(y_true_gene_m2, y_pred_gene_m2, zero_division=0)
        specificity2 = precision2  # For variant calling, specificity ≈ precision
        
        results.append({
            'Gene': gene,
            'Precision': precision1,
            'Recall': recall1,
            'Specificity': specificity1,
            'Precision_Bowtie2': precision2,
            'Recall_Bowtie2': recall2,
            'Specificity_Bowtie2': specificity2
        })
    
    # Calculate overall metrics using weighted average
    overall_precision1 = precision_score(y_true_method1, y_pred_method1, labels=genes, average='weighted', zero_division=0)
    overall_recall1 = recall_score(y_true_method1, y_pred_method1, labels=genes, average='weighted', zero_division=0)
    overall_specificity1 = overall_precision1
    
    overall_precision2 = precision_score(y_true_method2, y_pred_method2, labels=genes, average='weighted', zero_division=0)
    overall_recall2 = recall_score(y_true_method2, y_pred_method2, labels=genes, average='weighted', zero_division=0)
    overall_specificity2 = overall_precision2
    
    results.append({
        'Gene': 'Overall',
        'Precision': overall_precision1,
        'Recall': overall_recall1,
        'Specificity': overall_specificity1,
        'Precision_Bowtie2': overall_precision2,
        'Recall_Bowtie2': overall_recall2,
        'Specificity_Bowtie2': overall_specificity2
    })
    
    return results

def create_variant_plot(error_rate):
    df = load_variant_metrics(error_rate)
    
    if df is None:
        print(f"Could not load variant metrics for error rate {error_rate}")
        return
    
    print(f"\n=== Processing Variant Calling Error Rate {error_rate} ===")
    
    results = calculate_variant_metrics(df)
    
    combined_results = []
    for result in results:
        combined_results.append({
            'Gene': result['Gene'],
            'Precision_Mapper': result['Precision'],
            'Recall_Mapper': result['Recall'],
            'Specificity_Mapper': result['Specificity'],
            'Precision_Bowtie2': result['Precision_Bowtie2'],
            'Recall_Bowtie2': result['Recall_Bowtie2'],
            'Specificity_Bowtie2': result['Specificity_Bowtie2']
        })
    
    create_horizontal_plot(combined_results, f"Variant Calling Performance - Error Rate {error_rate*100:.1f}%", 
                          f"plots/variant_calling_{str(error_rate).replace('.', '_')}.png")
    save_tsv_output(combined_results, error_rate)
    
    print(f"\n=== VARIANT CALLING RESULTS FOR ERROR RATE {error_rate} ===")
    for result in combined_results:
        print(f"{result['Gene']:<12} - Mapper: P={result['Precision_Mapper']:.3f}, R={result['Recall_Mapper']:.3f}, S={result['Specificity_Mapper']:.3f}")
        print(f"{'':>12}   Bowtie2: P={result['Precision_Bowtie2']:.3f}, R={result['Recall_Bowtie2']:.3f}, S={result['Specificity_Bowtie2']:.3f}")

def create_horizontal_plot(results, title, filename):
    metrics_df = pd.DataFrame(results)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=300)
    
    genes = metrics_df['Gene'].tolist()
    y_pos = np.arange(len(genes))
    bar_height = 0.35
    
    for i, metric in enumerate(['Precision', 'Recall', 'Specificity']):
        axes[i].barh(y_pos - bar_height/2, metrics_df[f'{metric}_Mapper'], bar_height, 
                    label='Mapper', color=COLORS[0])
        axes[i].barh(y_pos + bar_height/2, metrics_df[f'{metric}_Bowtie2'], bar_height, 
                    label='Bowtie2', color=COLORS[1])
        axes[i].set_title(metric, fontweight='bold', fontsize=20)
        axes[i].set_yticks(y_pos)
        axes[i].set_yticklabels(genes, fontweight='bold')
        axes[i].set_xlim(0, 1.1)
        axes[i].grid(axis='x', alpha=0.3)
        axes[i].invert_yaxis()
    
    fig.suptitle(title, fontsize=24, fontweight='bold')
    
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"Created plot: {filename}")

def save_tsv_output(results, error_rate):
    df = pd.DataFrame(results)
    filename = f"plots/variant_calling_{str(error_rate).replace('.', '_')}.tsv"
    df.to_csv(filename, sep='\t', index=False)
    print(f"Created TSV: {filename}")

def main():
    error_rates = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]
    for rate in error_rates:
        create_variant_plot(rate)

if __name__ == "__main__":
    main()