#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from pathlib import Path

COLORS = ['#FF8C00', '#008080', '#4CAF50']
METHOD_COLORS = {'mapper': '#FF8C00', 'bowtie2': '#008080', 'freebayes': '#E91E63', 'bcftools': '#2196F3'}

plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14
})

def load_confusion_matrix(filepath):
    if not os.path.exists(filepath):
        return None
    return pd.read_csv(filepath, sep='\t', index_col=0)

def calculate_overall_metrics(confusion_matrix):
    if confusion_matrix is None:
        return None
    gene_classes = [col for col in confusion_matrix.columns if col != 'NONE']
    precisions, recalls, specificities = [], [], []
    for gene in gene_classes:
        tp = confusion_matrix.loc[gene, gene] if gene in confusion_matrix.index else 0
        fp = confusion_matrix[gene].sum() - tp
        fn = confusion_matrix.loc[gene].sum() - tp if gene in confusion_matrix.index else 0
        tn = confusion_matrix.values.sum() - (tp + fp + fn)
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 1.0
        precisions.append(precision); recalls.append(recall); specificities.append(specificity)
    return {'precision': np.mean(precisions), 'recall': np.mean(recalls), 'specificity': np.mean(specificities)}

def collect_metrics_across_iterations(method, iterations, output_dir):
    read_mapping_data = {
        'iterations': [],
        'mapper_precision': [], 'mapper_recall': [], 'mapper_specificity': [],
        'bowtie2_precision': [], 'bowtie2_recall': [], 'bowtie2_specificity': []
    }
    for i in range(1, iterations + 1):
        iter_dir = f"{output_dir}/iter_{i}"
        mapper_cm_file = f"{iter_dir}/read_mapping_{method}_iter_{i}_mapper_confusion_matrix.tsv"
        bowtie2_cm_file = f"{iter_dir}/read_mapping_{method}_iter_{i}_bowtie2_confusion_matrix.tsv"
        mapper_cm = load_confusion_matrix(mapper_cm_file)
        bowtie2_cm = load_confusion_matrix(bowtie2_cm_file)
        mapper_metrics = calculate_overall_metrics(mapper_cm)
        bowtie2_metrics = calculate_overall_metrics(bowtie2_cm)
        if mapper_metrics and bowtie2_metrics:
            read_mapping_data['iterations'].append(i)
            read_mapping_data['mapper_precision'].append(mapper_metrics['precision'])
            read_mapping_data['mapper_recall'].append(mapper_metrics['recall'])
            read_mapping_data['mapper_specificity'].append(mapper_metrics['specificity'])
            read_mapping_data['bowtie2_precision'].append(bowtie2_metrics['precision'])
            read_mapping_data['bowtie2_recall'].append(bowtie2_metrics['recall'])
            read_mapping_data['bowtie2_specificity'].append(bowtie2_metrics['specificity'])
    return read_mapping_data

def plot_read_mapping_summary(data, method, output_file):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    iterations = data['iterations']
    ax1.plot(iterations, data['mapper_precision'], 'o-', color=METHOD_COLORS['mapper'], linewidth=3, markersize=8, label='Mapper')
    ax1.plot(iterations, data['bowtie2_precision'], 's-', color=METHOD_COLORS['bowtie2'], linewidth=3, markersize=8, label='Bowtie2')
    ax1.set_xlabel('Iteration'); ax1.set_ylabel('Precision'); ax1.set_title('Read Mapping Precision'); ax1.grid(True, alpha=0.3); ax1.set_ylim(0, 1.2); ax1.set_xticks(iterations); ax1.set_yticks(np.linspace(0, 1, 6))
    ax2.plot(iterations, data['mapper_recall'], 'o-', color=METHOD_COLORS['mapper'], linewidth=3, markersize=8, label='Mapper')
    ax2.plot(iterations, data['bowtie2_recall'], 's-', color=METHOD_COLORS['bowtie2'], linewidth=3, markersize=8, label='Bowtie2')
    ax2.set_xlabel('Iteration'); ax2.set_ylabel('Recall'); ax2.set_title('Read Mapping Recall'); ax2.grid(True, alpha=0.3); ax2.set_ylim(0, 1.2); ax2.set_xticks(iterations); ax2.set_yticks(np.linspace(0, 1, 6))
    ax3.plot(iterations, data['mapper_specificity'], 'o-', color=METHOD_COLORS['mapper'], linewidth=3, markersize=8, label='Mapper')
    ax3.plot(iterations, data['bowtie2_specificity'], 's-', color=METHOD_COLORS['bowtie2'], linewidth=3, markersize=8, label='Bowtie2')
    ax3.set_xlabel('Iteration'); ax3.set_ylabel('Specificity'); ax3.set_title('Read Mapping Specificity'); ax3.grid(True, alpha=0.3); ax3.set_ylim(0, 1.2); ax3.set_xticks(iterations); ax3.set_yticks(np.linspace(0, 1, 6))
    handles, labels = ax1.get_legend_handles_labels(); fig.legend(handles, labels, loc='upper right')
    plt.suptitle(f'{method.capitalize()} - Read Mapping Performance Across Iterations', fontsize=22)
    plt.tight_layout(); plt.savefig(output_file, dpi=300, bbox_inches='tight'); plt.close()

def main():
    parser = argparse.ArgumentParser(description='Create summary plots across iterations')
    parser.add_argument('--method', required=True, choices=['freebayes', 'bcftools'])
    parser.add_argument('--iterations', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--plot-dir', required=True)
    args = parser.parse_args()
    os.makedirs(args.plot_dir, exist_ok=True)
    print(f"Generating summary plots for {args.method} across {args.iterations} iterations...")
    read_mapping_data = collect_metrics_across_iterations(args.method, args.iterations, args.output_dir)
    if read_mapping_data['iterations']:
        plot_read_mapping_summary(read_mapping_data, args.method, f'{args.plot_dir}/read_mapping_summary_{args.method}.png')
        print('✓ Read mapping summary plot generated')
    else:
        print('Warning: No read mapping data found')
    print(f"✓ Summary plots saved in: {args.plot_dir}")

if __name__ == '__main__':
    main()
