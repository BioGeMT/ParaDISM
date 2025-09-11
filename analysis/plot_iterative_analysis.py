#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

COLORS = ['#FF8C00', '#008080', '#4CAF50']
METHOD_COLORS = {'mapper': '#FF8C00', 'bowtie2': '#008080', 'method1': '#FF8C00', 'method2': '#008080'}

plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 12
})

def load_confusion_matrix(filepath):
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None
    return pd.read_csv(filepath, sep='\t', index_col=0)

def calculate_multiclass_metrics(confusion_matrix):
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
    results = []
    for gene in gene_classes:
        tp = confusion_matrix.loc[gene, gene] if gene in confusion_matrix.index else 0
        fp = confusion_matrix[gene].sum() - tp
        fn = confusion_matrix.loc[gene].sum() - tp if gene in confusion_matrix.index else 0
        tn = confusion_matrix.values.sum() - (tp + fp + fn)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 1.0
        results.append({'Gene': gene, 'Precision': precision, 'Recall': recall, 'Specificity': specificity})
    return pd.DataFrame(results)

def save_confusion_matrix_csv(confusion_matrix, output_file):
    confusion_matrix.to_csv(output_file)

def create_combined_read_mapping_plot(df_mapper, df_bowtie2, title, output_file):
    genes = df_mapper['Gene'].tolist()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 8))
    y_pos = np.arange(len(genes))
    ax1.barh(y_pos - 0.2, df_mapper['Precision'], 0.4, label='Mapper', color=METHOD_COLORS['mapper'], alpha=0.8)
    ax1.barh(y_pos + 0.2, df_bowtie2['Precision'], 0.4, label='Bowtie2', color=METHOD_COLORS['bowtie2'], alpha=0.8)
    ax1.set_xlabel('Precision'); ax1.set_ylabel('Gene'); ax1.set_title('Precision')
    ax1.set_yticks(y_pos); ax1.set_yticklabels(genes); ax1.grid(axis='x', alpha=0.3); ax1.set_xlim(0, 1)
    ax2.barh(y_pos - 0.2, df_mapper['Recall'], 0.4, label='Mapper', color=METHOD_COLORS['mapper'], alpha=0.8)
    ax2.barh(y_pos + 0.2, df_bowtie2['Recall'], 0.4, label='Bowtie2', color=METHOD_COLORS['bowtie2'], alpha=0.8)
    ax2.set_xlabel('Recall'); ax2.set_title('Recall')
    ax2.set_yticks(y_pos); ax2.set_yticklabels(genes); ax2.grid(axis='x', alpha=0.3); ax2.set_xlim(0, 1)
    ax3.barh(y_pos - 0.2, df_mapper['Specificity'], 0.4, label='Mapper', color=METHOD_COLORS['mapper'], alpha=0.8)
    ax3.barh(y_pos + 0.2, df_bowtie2['Specificity'], 0.4, label='Bowtie2', color=METHOD_COLORS['bowtie2'], alpha=0.8)
    ax3.set_xlabel('Specificity'); ax3.set_title('Specificity')
    ax3.set_yticks(y_pos); ax3.set_yticklabels(genes); ax3.grid(axis='x', alpha=0.3); ax3.set_xlim(0, 1)
    handles, labels = ax1.get_legend_handles_labels(); plt.suptitle(title, fontsize=22)
    fig.legend(handles, labels, loc='upper right'); plt.tight_layout(); plt.savefig(output_file, dpi=300, bbox_inches='tight'); plt.close()

def main():
    parser = argparse.ArgumentParser(description='Plot iterative analysis results')
    parser.add_argument('--iteration', type=int, required=True)
    parser.add_argument('--method', required=True, choices=['freebayes', 'bcftools'])
    parser.add_argument('--iter-dir', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    method_name = args.method.capitalize()
    print(f"Processing {method_name} iteration {args.iteration}...")
    mapper_cm_file = f"{args.iter_dir}/read_mapping_{args.method}_iter_{args.iteration}_mapper_confusion_matrix.tsv"
    bowtie2_cm_file = f"{args.iter_dir}/read_mapping_{args.method}_iter_{args.iteration}_bowtie2_confusion_matrix.tsv"
    mapper_cm = load_confusion_matrix(mapper_cm_file)
    bowtie2_cm = load_confusion_matrix(bowtie2_cm_file)
    if mapper_cm is not None and bowtie2_cm is not None:
        save_confusion_matrix_csv(mapper_cm, f'{args.output_dir}/read_mapping_mapper_confusion_matrix.csv')
        save_confusion_matrix_csv(bowtie2_cm, f'{args.output_dir}/read_mapping_bowtie2_confusion_matrix.csv')
        mapper_metrics = calculate_multiclass_metrics(mapper_cm)
        bowtie2_metrics = calculate_multiclass_metrics(bowtie2_cm)
        create_combined_read_mapping_plot(mapper_metrics, bowtie2_metrics,
                                         f'{method_name} Iteration {args.iteration} - Read Mapping Metrics',
                                         f'{args.output_dir}/read_mapping_metrics.png')
        print("✓ Read mapping plots generated and confusion matrices saved as CSV")
    else:
        print("Warning: Could not load read mapping confusion matrices")
    print(f"✓ All plots generated in: {args.output_dir}")

if __name__ == '__main__':
    main()

