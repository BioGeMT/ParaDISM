#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report
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

def load_confusion_matrix(error_rate, method):
    rate_str = str(error_rate).replace('.', '_')
    file_path = f"error_rate_results/err_{rate_str}/analysis/read_mapping_err_{rate_str}_{method}_confusion_matrix.tsv"
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    return pd.read_csv(file_path, sep='\t', index_col=0)

def calculate_multiclass_metrics(confusion_matrix):
    # convert confusion matrix to sklearn format
    y_true = []
    y_pred = []
    
    for true_gene in confusion_matrix.index:
        for pred_gene in confusion_matrix.columns:
            count = confusion_matrix.loc[true_gene, pred_gene]
            y_true.extend([true_gene] * int(count))
            y_pred.extend([pred_gene] * int(count))
    
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    
    gene_classes = [col for col in confusion_matrix.columns if col != 'NONE']
    
    # use sklearn for precision and recall
    report = classification_report(y_true, y_pred, labels=gene_classes, output_dict=True, zero_division=0)
    
    print("\nSKLEARN CLASSIFICATION REPORT:")
    print(classification_report(y_true, y_pred, labels=gene_classes, target_names=gene_classes, zero_division=0))
    
    # calculate specificity manually (not in sklearn)
    specificity_per_class = []
    for gene in gene_classes:
        tp = confusion_matrix.loc[gene, gene] if gene in confusion_matrix.index else 0
        fp = confusion_matrix[gene].sum() - tp
        fn = confusion_matrix.loc[gene].sum() - tp if gene in confusion_matrix.index else 0
        tn = confusion_matrix.values.sum() - (tp + fp + fn)
        
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 1.0
        specificity_per_class.append(specificity)
    
    results = []
    
    for i, gene in enumerate(gene_classes):
        results.append({
            'Gene': gene,
            'Precision': report[gene]['precision'],
            'Recall': report[gene]['recall'],
            'Specificity': specificity_per_class[i]
        })
    
    # overall results using macro average
    results.append({
        'Gene': 'Overall',
        'Precision': report['macro avg']['precision'],
        'Recall': report['macro avg']['recall'],
        'Specificity': np.mean(specificity_per_class)
    })
    
    return results

def create_multiclass_plot(error_rate):
    mapper_cm = load_confusion_matrix(error_rate, 'mapper')
    bowtie2_cm = load_confusion_matrix(error_rate, 'bowtie2')
    
    if mapper_cm is None or bowtie2_cm is None:
        print(f"Could not load confusion matrices for error rate {error_rate}")
        return
    
    print(f"\n=== Processing Error Rate {error_rate} ===")
    
    mapper_results = calculate_multiclass_metrics(mapper_cm)
    bowtie2_results = calculate_multiclass_metrics(bowtie2_cm)
    
    combined_results = []
    for mapper_result, bowtie2_result in zip(mapper_results, bowtie2_results):
        combined_results.append({
            'Gene': mapper_result['Gene'],
            'Precision_Mapper': mapper_result['Precision'],
            'Recall_Mapper': mapper_result['Recall'],
            'Specificity_Mapper': mapper_result['Specificity'],
            'Precision_Bowtie2': bowtie2_result['Precision'],
            'Recall_Bowtie2': bowtie2_result['Recall'],
            'Specificity_Bowtie2': bowtie2_result['Specificity']
        })
    
    create_horizontal_plot(combined_results, "Read Mapping Performance", 
                          f"multiclass_read_mapping_{str(error_rate).replace('.', '_')}.png")
    save_tsv_output(combined_results, error_rate)
    
    print(f"\n=== MULTICLASS RESULTS FOR ERROR RATE {error_rate} ===")
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
    filename = f"multiclass_read_mapping_{str(error_rate).replace('.', '_')}.tsv"
    df.to_csv(filename, sep='\t', index=False)
    print(f"Created TSV: {filename}")

def main():
    error_rates = [0.001, 0.005]
    for rate in error_rates:
        create_multiclass_plot(rate)

if __name__ == "__main__":
    main()