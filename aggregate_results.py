#!/usr/bin/env python3
"""
Aggregate results across all seeds and create summary plots with error bars.
Creates per-gene plots (like per-seed plots) but aggregated across multiple seeds with error bars.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix


def parse_read_mapping_summary_per_gene(csv_path):
    """
    Parse read mapping summary CSV and extract per-gene metrics.
    Returns dict with per-gene metrics for mapper and direct.
    """
    df = pd.read_csv(csv_path)

    # Get column names
    ground_truth = df['Ground_Truth']
    mapper_pred = df['Mapper_Prediction']

    # Find direct prediction column
    direct_pred_col = None
    for col in df.columns:
        if col.endswith('_Prediction') and col != 'Mapper_Prediction':
            direct_pred_col = col

    if direct_pred_col is None:
        direct_pred = df.iloc[:, 3]
    else:
        direct_pred = df[direct_pred_col]

    # Get all labels (genes)
    all_labels = sorted(set(ground_truth) | set(mapper_pred) | set(direct_pred))

    # Calculate per-gene metrics
    def calc_per_gene_metrics(y_true, y_pred, labels):
        """Calculate precision, recall, specificity per gene."""
        # Calculate confusion matrix
        cm = confusion_matrix(y_true, y_pred, labels=labels)

        metrics = {}
        for i, label in enumerate(labels):
            if label == 'NONE':
                continue

            tp = cm[i, i]
            fp = np.sum(cm[:, i]) - tp
            fn = np.sum(cm[i, :]) - tp
            tn = np.sum(cm) - tp - fp - fn

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

            metrics[label] = {
                'precision': precision,
                'recall': recall,
                'specificity': specificity
            }

        # Calculate overall weighted metrics
        prec, rec, _, support = precision_recall_fscore_support(
            y_true, y_pred, labels=labels, average='weighted', zero_division=0
        )

        # Overall specificity
        specificities = []
        supports = []
        for i, label in enumerate(labels):
            tn = np.sum(cm) - np.sum(cm[i, :]) - np.sum(cm[:, i]) + cm[i, i]
            fp = np.sum(cm[:, i]) - cm[i, i]
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
            specificities.append(spec)
            supports.append(np.sum(cm[i, :]))
        total_support = np.sum(supports)
        overall_spec = np.sum(np.array(specificities) * np.array(supports)) / total_support if total_support > 0 else 0.0

        metrics['Overall'] = {
            'precision': prec,
            'recall': rec,
            'specificity': overall_spec
        }

        return metrics

    mapper_metrics = calc_per_gene_metrics(ground_truth, mapper_pred, all_labels)
    direct_metrics = calc_per_gene_metrics(ground_truth, direct_pred, all_labels)

    return {
        'mapper': mapper_metrics,
        'direct': direct_metrics
    }


def parse_variant_calling_summary_per_gene(csv_path):
    """
    Parse variant calling summary CSV and extract per-gene metrics.
    Returns dict with per-gene metrics for combined and direct.
    """
    if not Path(csv_path).exists():
        return None

    df = pd.read_csv(csv_path)

    combined_metrics = {}
    direct_metrics = {}

    for _, row in df.iterrows():
        gene = row['Gene']
        combined_metrics[gene] = {
            'precision': row['Combined_Precision'],
            'recall': row['Combined_Recall'],
            'specificity': row['Combined_Specificity']
        }
        direct_metrics[gene] = {
            'precision': row['Direct_Precision'],
            'recall': row['Direct_Recall'],
            'specificity': row['Direct_Specificity']
        }

    return {
        'combined': combined_metrics,
        'direct': direct_metrics
    }


def aggregate_read_mapping_results(sim_output_base, aligners, seed_start, seed_end):
    """Aggregate read mapping results per-gene across all seeds."""
    # Structure: {aligner: {approach: {gene: {metric: [values across seeds]}}}}
    results = {aligner: {'mapper': defaultdict(lambda: defaultdict(list)),
                         'direct': defaultdict(lambda: defaultdict(list))}
               for aligner in aligners}

    for seed in range(seed_start, seed_end + 1):
        seed_dir = Path(sim_output_base) / f"seed_{seed}"

        for aligner in aligners:
            csv_path = seed_dir / aligner / "read_mapping_analysis" / f"seed_{seed}_{aligner}_summary.csv"

            if csv_path.exists():
                try:
                    metrics = parse_read_mapping_summary_per_gene(csv_path)

                    # Collect mapper metrics
                    for gene, gene_metrics in metrics['mapper'].items():
                        for metric_name, value in gene_metrics.items():
                            results[aligner]['mapper'][gene][metric_name].append(value)

                    # Collect direct metrics
                    for gene, gene_metrics in metrics['direct'].items():
                        for metric_name, value in gene_metrics.items():
                            results[aligner]['direct'][gene][metric_name].append(value)

                except Exception as e:
                    print(f"Warning: Failed to parse {csv_path}: {e}", file=sys.stderr)
                    continue

    # Calculate means and stds
    aggregated = {}
    for aligner in aligners:
        aggregated[aligner] = {'mapper': {}, 'direct': {}}

        for approach in ['mapper', 'direct']:
            for gene, metrics in results[aligner][approach].items():
                aggregated[aligner][approach][gene] = {
                    'precision_mean': np.mean(metrics['precision']),
                    'precision_std': np.std(metrics['precision']),
                    'recall_mean': np.mean(metrics['recall']),
                    'recall_std': np.std(metrics['recall']),
                    'specificity_mean': np.mean(metrics['specificity']),
                    'specificity_std': np.std(metrics['specificity']),
                    'n_seeds': len(metrics['precision'])
                }

    return aggregated


def aggregate_variant_calling_results(sim_output_base, aligners, seed_start, seed_end):
    """Aggregate variant calling results per-gene across all seeds."""
    # Structure: {aligner: {approach: {gene: {metric: [values across seeds]}}}}
    results = {aligner: {'combined': defaultdict(lambda: defaultdict(list)),
                         'direct': defaultdict(lambda: defaultdict(list))}
               for aligner in aligners}

    for seed in range(seed_start, seed_end + 1):
        seed_dir = Path(sim_output_base) / f"seed_{seed}"

        for aligner in aligners:
            csv_path = seed_dir / aligner / "variant_calling_analysis" / f"paradism_seed_{seed}_{aligner}_summary.csv"

            if csv_path.exists():
                try:
                    metrics = parse_variant_calling_summary_per_gene(csv_path)
                    if metrics is None:
                        continue

                    # Collect combined metrics
                    for gene, gene_metrics in metrics['combined'].items():
                        for metric_name, value in gene_metrics.items():
                            results[aligner]['combined'][gene][metric_name].append(value)

                    # Collect direct metrics
                    for gene, gene_metrics in metrics['direct'].items():
                        for metric_name, value in gene_metrics.items():
                            results[aligner]['direct'][gene][metric_name].append(value)

                except Exception as e:
                    print(f"Warning: Failed to parse {csv_path}: {e}", file=sys.stderr)
                    continue

    # Calculate means and stds
    aggregated = {}
    for aligner in aligners:
        aggregated[aligner] = {'combined': {}, 'direct': {}}

        for approach in ['combined', 'direct']:
            for gene, metrics in results[aligner][approach].items():
                aggregated[aligner][approach][gene] = {
                    'precision_mean': np.mean(metrics['precision']),
                    'precision_std': np.std(metrics['precision']),
                    'recall_mean': np.mean(metrics['recall']),
                    'recall_std': np.std(metrics['recall']),
                    'specificity_mean': np.mean(metrics['specificity']),
                    'specificity_std': np.std(metrics['specificity']),
                    'n_seeds': len(metrics['precision'])
                }

    return aggregated


def create_read_mapping_plot(aggregated, aligner, output_path):
    """Create read mapping comparison plot for a single aligner with per-gene breakdown and error bars."""
    if aligner not in aggregated or 'mapper' not in aggregated[aligner]:
        print(f"Warning: No data found for {aligner}", file=sys.stderr)
        return

    # Get all genes
    genes = [g for g in aggregated[aligner]['mapper'].keys() if g != 'Overall']
    genes = sorted(genes)
    display_labels = genes + ['Overall']

    # Set font sizes
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'font.weight': 'bold',
        'axes.titleweight': 'bold',
        'axes.labelweight': 'bold'
    })

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    y_pos = np.arange(len(display_labels))
    bar_height = 0.35

    colors = {'mapper': 'red', 'direct': '#4682B4'}

    # Extract data for mapper
    mapper_prec = []
    mapper_prec_err = []
    mapper_rec = []
    mapper_rec_err = []
    mapper_spec = []
    mapper_spec_err = []

    for gene in display_labels:
        if gene in aggregated[aligner]['mapper']:
            mapper_prec.append(aggregated[aligner]['mapper'][gene]['precision_mean'])
            mapper_prec_err.append(aggregated[aligner]['mapper'][gene]['precision_std'])
            mapper_rec.append(aggregated[aligner]['mapper'][gene]['recall_mean'])
            mapper_rec_err.append(aggregated[aligner]['mapper'][gene]['recall_std'])
            mapper_spec.append(aggregated[aligner]['mapper'][gene]['specificity_mean'])
            mapper_spec_err.append(aggregated[aligner]['mapper'][gene]['specificity_std'])
        else:
            mapper_prec.append(0)
            mapper_prec_err.append(0)
            mapper_rec.append(0)
            mapper_rec_err.append(0)
            mapper_spec.append(0)
            mapper_spec_err.append(0)

    # Extract data for direct
    direct_prec = []
    direct_prec_err = []
    direct_rec = []
    direct_rec_err = []
    direct_spec = []
    direct_spec_err = []

    for gene in display_labels:
        if gene in aggregated[aligner]['direct']:
            direct_prec.append(aggregated[aligner]['direct'][gene]['precision_mean'])
            direct_prec_err.append(aggregated[aligner]['direct'][gene]['precision_std'])
            direct_rec.append(aggregated[aligner]['direct'][gene]['recall_mean'])
            direct_rec_err.append(aggregated[aligner]['direct'][gene]['recall_std'])
            direct_spec.append(aggregated[aligner]['direct'][gene]['specificity_mean'])
            direct_spec_err.append(aggregated[aligner]['direct'][gene]['specificity_std'])
        else:
            direct_prec.append(0)
            direct_prec_err.append(0)
            direct_rec.append(0)
            direct_rec_err.append(0)
            direct_spec.append(0)
            direct_spec_err.append(0)

    # Panel 1: Precision
    axes[0].barh(y_pos - bar_height/2, mapper_prec, bar_height,
                 xerr=mapper_prec_err, label='paradism', color=colors['mapper'], capsize=5)
    axes[0].barh(y_pos + bar_height/2, direct_prec, bar_height,
                 xerr=direct_prec_err, label=aligner, color=colors['direct'], capsize=5)

    # Panel 2: Recall
    axes[1].barh(y_pos - bar_height/2, mapper_rec, bar_height,
                 xerr=mapper_rec_err, label='paradism', color=colors['mapper'], capsize=5)
    axes[1].barh(y_pos + bar_height/2, direct_rec, bar_height,
                 xerr=direct_rec_err, label=aligner, color=colors['direct'], capsize=5)

    # Panel 3: Specificity
    axes[2].barh(y_pos - bar_height/2, mapper_spec, bar_height,
                 xerr=mapper_spec_err, label='paradism', color=colors['mapper'], capsize=5)
    axes[2].barh(y_pos + bar_height/2, direct_spec, bar_height,
                 xerr=direct_spec_err, label=aligner, color=colors['direct'], capsize=5)

    # Format axes
    for idx, (ax, title) in enumerate(zip(axes, ['Precision', 'Recall', 'Specificity'])):
        ax.set_yticks(y_pos)
        ax.set_yticklabels(display_labels, weight='bold')
        ax.set_xlabel(title, weight='bold')
        ax.set_title(title, weight='bold')
        ax.set_xlim([0, 1.05])
        ax.grid(axis='x', alpha=0.3)

    # Add legend
    handles, labels_legend = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_legend, loc='upper right', bbox_to_anchor=(1.02, 1.02),
               fontsize=16, frameon=True, fancybox=True, shadow=True)

    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_variant_calling_plot(aggregated, aligner, output_path):
    """Create variant calling comparison plot for a single aligner with per-gene breakdown and error bars."""
    if aligner not in aggregated or 'combined' not in aggregated[aligner]:
        print(f"Warning: No data found for {aligner}", file=sys.stderr)
        return

    # Get all genes
    genes = [g for g in aggregated[aligner]['combined'].keys() if g != 'Overall']
    genes = sorted(genes)
    display_labels = genes + ['Overall']

    # Set font sizes
    plt.rcParams.update({
        'font.size': 14,
        'axes.titlesize': 18,
        'axes.labelsize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'font.weight': 'bold',
        'axes.titleweight': 'bold',
        'axes.labelweight': 'bold'
    })

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    y_pos = np.arange(len(display_labels))
    bar_height = 0.35

    colors = {'combined': 'red', 'direct': '#4682B4'}

    # Extract data for combined
    combined_prec = []
    combined_prec_err = []
    combined_rec = []
    combined_rec_err = []
    combined_spec = []
    combined_spec_err = []

    for gene in display_labels:
        if gene in aggregated[aligner]['combined']:
            combined_prec.append(aggregated[aligner]['combined'][gene]['precision_mean'])
            combined_prec_err.append(aggregated[aligner]['combined'][gene]['precision_std'])
            combined_rec.append(aggregated[aligner]['combined'][gene]['recall_mean'])
            combined_rec_err.append(aggregated[aligner]['combined'][gene]['recall_std'])
            combined_spec.append(aggregated[aligner]['combined'][gene]['specificity_mean'])
            combined_spec_err.append(aggregated[aligner]['combined'][gene]['specificity_std'])
        else:
            combined_prec.append(0)
            combined_prec_err.append(0)
            combined_rec.append(0)
            combined_rec_err.append(0)
            combined_spec.append(0)
            combined_spec_err.append(0)

    # Extract data for direct
    direct_prec = []
    direct_prec_err = []
    direct_rec = []
    direct_rec_err = []
    direct_spec = []
    direct_spec_err = []

    for gene in display_labels:
        if gene in aggregated[aligner]['direct']:
            direct_prec.append(aggregated[aligner]['direct'][gene]['precision_mean'])
            direct_prec_err.append(aggregated[aligner]['direct'][gene]['precision_std'])
            direct_rec.append(aggregated[aligner]['direct'][gene]['recall_mean'])
            direct_rec_err.append(aggregated[aligner]['direct'][gene]['recall_std'])
            direct_spec.append(aggregated[aligner]['direct'][gene]['specificity_mean'])
            direct_spec_err.append(aggregated[aligner]['direct'][gene]['specificity_std'])
        else:
            direct_prec.append(0)
            direct_prec_err.append(0)
            direct_rec.append(0)
            direct_rec_err.append(0)
            direct_spec.append(0)
            direct_spec_err.append(0)

    # Panel 1: Precision
    axes[0].barh(y_pos - bar_height/2, combined_prec, bar_height,
                 xerr=combined_prec_err, label='paradism', color=colors['combined'], capsize=5)
    axes[0].barh(y_pos + bar_height/2, direct_prec, bar_height,
                 xerr=direct_prec_err, label=aligner, color=colors['direct'], capsize=5)

    # Panel 2: Recall
    axes[1].barh(y_pos - bar_height/2, combined_rec, bar_height,
                 xerr=combined_rec_err, label='paradism', color=colors['combined'], capsize=5)
    axes[1].barh(y_pos + bar_height/2, direct_rec, bar_height,
                 xerr=direct_rec_err, label=aligner, color=colors['direct'], capsize=5)

    # Panel 3: Specificity
    axes[2].barh(y_pos - bar_height/2, combined_spec, bar_height,
                 xerr=combined_spec_err, label='paradism', color=colors['combined'], capsize=5)
    axes[2].barh(y_pos + bar_height/2, direct_spec, bar_height,
                 xerr=direct_spec_err, label=aligner, color=colors['direct'], capsize=5)

    # Format axes
    for idx, (ax, title) in enumerate(zip(axes, ['Precision', 'Recall', 'Specificity'])):
        ax.set_yticks(y_pos)
        ax.set_yticklabels(display_labels, weight='bold')
        ax.set_xlabel(title, weight='bold')
        ax.set_title(title, weight='bold')
        ax.set_xlim([0, 1.05])
        ax.grid(axis='x', alpha=0.3)

    # Add legend
    handles, labels_legend = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_legend, loc='upper right', bbox_to_anchor=(1.02, 1.02),
               fontsize=16, frameon=True, fancybox=True, shadow=True)

    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_read_mapping_overall_comparison(aggregated, aligners, output_path):
    """Create overall read mapping comparison across all aligners with vertical bars."""
    # Use BWA-MEM2 as the ParaDISM baseline
    paradism_aligner = 'bwa-mem2'

    if paradism_aligner not in aggregated:
        print(f"Warning: {paradism_aligner} not found for ParaDISM baseline", file=sys.stderr)
        return

    # Check if 'Overall' exists
    if 'Overall' not in aggregated[paradism_aligner]['mapper']:
        print("Warning: No Overall metrics found", file=sys.stderr)
        return

    # Set font sizes
    plt.rcParams.update({
        'font.size': 11,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 11,
        'font.weight': 'bold',
        'axes.titleweight': 'bold',
        'axes.labelweight': 'bold'
    })

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Prepare data: ParaDISM + 3 direct aligners
    labels = ['paradism', 'bwa-mem2', 'bowtie2', 'minimap2']
    colors = ['red', '#4682B4', '#4682B4', '#4682B4']

    # Extract ParaDISM metrics (mapper approach)
    paradism_prec = aggregated[paradism_aligner]['mapper']['Overall']['precision_mean']
    paradism_prec_err = aggregated[paradism_aligner]['mapper']['Overall']['precision_std']
    paradism_rec = aggregated[paradism_aligner]['mapper']['Overall']['recall_mean']
    paradism_rec_err = aggregated[paradism_aligner]['mapper']['Overall']['recall_std']
    paradism_spec = aggregated[paradism_aligner]['mapper']['Overall']['specificity_mean']
    paradism_spec_err = aggregated[paradism_aligner]['mapper']['Overall']['specificity_std']

    # Extract direct aligner metrics
    precision_values = [paradism_prec]
    precision_errors = [paradism_prec_err]
    recall_values = [paradism_rec]
    recall_errors = [paradism_rec_err]
    specificity_values = [paradism_spec]
    specificity_errors = [paradism_spec_err]

    for aligner in aligners:
        if aligner in aggregated and 'Overall' in aggregated[aligner]['direct']:
            precision_values.append(aggregated[aligner]['direct']['Overall']['precision_mean'])
            precision_errors.append(aggregated[aligner]['direct']['Overall']['precision_std'])
            recall_values.append(aggregated[aligner]['direct']['Overall']['recall_mean'])
            recall_errors.append(aggregated[aligner]['direct']['Overall']['recall_std'])
            specificity_values.append(aggregated[aligner]['direct']['Overall']['specificity_mean'])
            specificity_errors.append(aggregated[aligner]['direct']['Overall']['specificity_std'])
        else:
            precision_values.append(0)
            precision_errors.append(0)
            recall_values.append(0)
            recall_errors.append(0)
            specificity_values.append(0)
            specificity_errors.append(0)

    x_pos = np.arange(len(labels))
    bar_width = 0.6

    # Panel 1: Precision
    axes[0].bar(x_pos, precision_values, bar_width, yerr=precision_errors,
                color=colors, capsize=15)
    axes[0].set_ylabel('Precision', weight='bold')
    axes[0].set_title('Precision', weight='bold')
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(labels, weight='bold')
    axes[0].set_ylim([0, 1.05])
    axes[0].grid(axis='y', alpha=0.3)

    # Panel 2: Recall
    axes[1].bar(x_pos, recall_values, bar_width, yerr=recall_errors,
                color=colors, capsize=15)
    axes[1].set_ylabel('Recall', weight='bold')
    axes[1].set_title('Recall', weight='bold')
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(labels, weight='bold')
    axes[1].set_ylim([0, 1.05])
    axes[1].grid(axis='y', alpha=0.3)

    # Panel 3: Specificity
    axes[2].bar(x_pos, specificity_values, bar_width, yerr=specificity_errors,
                color=colors, capsize=15)
    axes[2].set_ylabel('Specificity', weight='bold')
    axes[2].set_title('Specificity', weight='bold')
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels(labels, weight='bold')
    axes[2].set_ylim([0, 1.05])
    axes[2].grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_variant_calling_overall_comparison(aggregated, aligners, output_path):
    """Create overall variant calling comparison across all aligners with vertical bars."""
    # Use BWA-MEM2 as the ParaDISM baseline
    paradism_aligner = 'bwa-mem2'

    if paradism_aligner not in aggregated:
        print(f"Warning: {paradism_aligner} not found for ParaDISM baseline", file=sys.stderr)
        return

    # Check if 'Overall' exists
    if 'Overall' not in aggregated[paradism_aligner]['combined']:
        print("Warning: No Overall metrics found", file=sys.stderr)
        return

    # Set font sizes
    plt.rcParams.update({
        'font.size': 11,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 11,
        'font.weight': 'bold',
        'axes.titleweight': 'bold',
        'axes.labelweight': 'bold'
    })

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Prepare data: ParaDISM + 3 direct aligners
    labels = ['paradism', 'bwa-mem2', 'bowtie2', 'minimap2']
    colors = ['red', '#4682B4', '#4682B4', '#4682B4']

    # Extract ParaDISM metrics (combined approach)
    paradism_prec = aggregated[paradism_aligner]['combined']['Overall']['precision_mean']
    paradism_prec_err = aggregated[paradism_aligner]['combined']['Overall']['precision_std']
    paradism_rec = aggregated[paradism_aligner]['combined']['Overall']['recall_mean']
    paradism_rec_err = aggregated[paradism_aligner]['combined']['Overall']['recall_std']
    paradism_spec = aggregated[paradism_aligner]['combined']['Overall']['specificity_mean']
    paradism_spec_err = aggregated[paradism_aligner]['combined']['Overall']['specificity_std']

    # Extract direct aligner metrics
    precision_values = [paradism_prec]
    precision_errors = [paradism_prec_err]
    recall_values = [paradism_rec]
    recall_errors = [paradism_rec_err]
    specificity_values = [paradism_spec]
    specificity_errors = [paradism_spec_err]

    for aligner in aligners:
        if aligner in aggregated and 'Overall' in aggregated[aligner]['direct']:
            precision_values.append(aggregated[aligner]['direct']['Overall']['precision_mean'])
            precision_errors.append(aggregated[aligner]['direct']['Overall']['precision_std'])
            recall_values.append(aggregated[aligner]['direct']['Overall']['recall_mean'])
            recall_errors.append(aggregated[aligner]['direct']['Overall']['recall_std'])
            specificity_values.append(aggregated[aligner]['direct']['Overall']['specificity_mean'])
            specificity_errors.append(aggregated[aligner]['direct']['Overall']['specificity_std'])
        else:
            precision_values.append(0)
            precision_errors.append(0)
            recall_values.append(0)
            recall_errors.append(0)
            specificity_values.append(0)
            specificity_errors.append(0)

    x_pos = np.arange(len(labels))
    bar_width = 0.6

    # Panel 1: Precision
    axes[0].bar(x_pos, precision_values, bar_width, yerr=precision_errors,
                color=colors, capsize=15)
    axes[0].set_ylabel('Precision', weight='bold')
    axes[0].set_title('Precision', weight='bold')
    axes[0].set_xticks(x_pos)
    axes[0].set_xticklabels(labels, weight='bold')
    axes[0].set_ylim([0, 1.05])
    axes[0].grid(axis='y', alpha=0.3)

    # Panel 2: Recall
    axes[1].bar(x_pos, recall_values, bar_width, yerr=recall_errors,
                color=colors, capsize=15)
    axes[1].set_ylabel('Recall', weight='bold')
    axes[1].set_title('Recall', weight='bold')
    axes[1].set_xticks(x_pos)
    axes[1].set_xticklabels(labels, weight='bold')
    axes[1].set_ylim([0, 1.05])
    axes[1].grid(axis='y', alpha=0.3)

    # Panel 3: Specificity
    axes[2].bar(x_pos, specificity_values, bar_width, yerr=specificity_errors,
                color=colors, capsize=15)
    axes[2].set_ylabel('Specificity', weight='bold')
    axes[2].set_title('Specificity', weight='bold')
    axes[2].set_xticks(x_pos)
    axes[2].set_xticklabels(labels, weight='bold')
    axes[2].set_ylim([0, 1.05])
    axes[2].grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate results across all seeds and create summary plots with per-gene breakdown'
    )
    parser.add_argument(
        '--sim-output-base',
        default='sim_output',
        help='Base directory for simulation outputs (default: sim_output)'
    )
    parser.add_argument(
        '--seed-start',
        type=int,
        default=1,
        help='Starting seed number (default: 1)'
    )
    parser.add_argument(
        '--seed-end',
        type=int,
        default=1000,
        help='Ending seed number (default: 1000)'
    )
    parser.add_argument(
        '--aligners',
        nargs='+',
        default=['bowtie2', 'bwa-mem2', 'minimap2'],
        help='List of aligners to aggregate (default: bowtie2 bwa-mem2 minimap2)'
    )
    parser.add_argument(
        '--output-dir',
        default='aggregated_results',
        help='Output directory for aggregated results (default: aggregated_results)'
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Aggregating read mapping results...")
    read_mapping_agg = aggregate_read_mapping_results(
        args.sim_output_base, args.aligners, args.seed_start, args.seed_end
    )

    print("Aggregating variant calling results...")
    variant_calling_agg = aggregate_variant_calling_results(
        args.sim_output_base, args.aligners, args.seed_start, args.seed_end
    )

    # Create CSV files
    print("Creating CSV files...")

    # Read mapping CSV
    rm_data = []
    for aligner in args.aligners:
        if aligner in read_mapping_agg:
            for approach in ['mapper', 'direct']:
                for gene, metrics in read_mapping_agg[aligner][approach].items():
                    rm_data.append({
                        'Aligner': aligner.upper(),
                        'Approach': approach.capitalize(),
                        'Gene': gene,
                        'Precision_Mean': metrics['precision_mean'],
                        'Precision_Std': metrics['precision_std'],
                        'Recall_Mean': metrics['recall_mean'],
                        'Recall_Std': metrics['recall_std'],
                        'Specificity_Mean': metrics['specificity_mean'],
                        'Specificity_Std': metrics['specificity_std'],
                        'N_Seeds': metrics['n_seeds']
                    })
    if rm_data:
        rm_df = pd.DataFrame(rm_data)
        rm_csv = output_dir / 'read_mapping_aggregated_summary.csv'
        rm_df.to_csv(rm_csv, index=False)
        print(f"  Saved: {rm_csv}")

    # Variant calling CSV
    vc_data = []
    for aligner in args.aligners:
        if aligner in variant_calling_agg:
            for approach in ['combined', 'direct']:
                for gene, metrics in variant_calling_agg[aligner][approach].items():
                    vc_data.append({
                        'Aligner': aligner.upper(),
                        'Approach': approach.capitalize(),
                        'Gene': gene,
                        'Precision_Mean': metrics['precision_mean'],
                        'Precision_Std': metrics['precision_std'],
                        'Recall_Mean': metrics['recall_mean'],
                        'Recall_Std': metrics['recall_std'],
                        'Specificity_Mean': metrics['specificity_mean'],
                        'Specificity_Std': metrics['specificity_std'],
                        'N_Seeds': metrics['n_seeds']
                    })
    if vc_data:
        vc_df = pd.DataFrame(vc_data)
        vc_csv = output_dir / 'variant_calling_aggregated_summary.csv'
        vc_df.to_csv(vc_csv, index=False)
        print(f"  Saved: {vc_csv}")

    # Create plots (one per aligner)
    print("Creating plots...")

    for aligner in args.aligners:
        if aligner in read_mapping_agg:
            rm_plot = output_dir / f'read_mapping_aggregated_{aligner}_comparison.png'
            create_read_mapping_plot(read_mapping_agg, aligner, rm_plot)
            print(f"  Saved: {rm_plot}")

        if aligner in variant_calling_agg:
            vc_plot = output_dir / f'variant_calling_aggregated_{aligner}_comparison.png'
            create_variant_calling_plot(variant_calling_agg, aligner, vc_plot)
            print(f"  Saved: {vc_plot}")

    # Create overall comparison plots (cross-aligner)
    print("Creating overall comparison plots...")

    if read_mapping_agg:
        rm_overall_plot = output_dir / 'read_mapping_overall_comparison.png'
        create_read_mapping_overall_comparison(read_mapping_agg, args.aligners, rm_overall_plot)
        print(f"  Saved: {rm_overall_plot}")

    if variant_calling_agg:
        vc_overall_plot = output_dir / 'variant_calling_overall_comparison.png'
        create_variant_calling_overall_comparison(variant_calling_agg, args.aligners, vc_overall_plot)
        print(f"  Saved: {vc_overall_plot}")

    print(f"\nAggregation complete! Results saved to: {output_dir}/")
    print(f"  Read mapping: {len([a for a in args.aligners if a in read_mapping_agg])} aligners")
    print(f"  Variant calling: {len([a for a in args.aligners if a in variant_calling_agg])} aligners")


if __name__ == "__main__":
    main()
