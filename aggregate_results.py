#!/usr/bin/env python3
"""
Aggregate results across all seeds and create summary plots with error bars.
Creates one read mapping plot and one variant calling plot with aggregated metrics.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix


def parse_read_mapping_summary(csv_path):
    """
    Parse read mapping summary CSV and extract overall metrics.
    Returns dict with precision, recall, specificity for mapper and direct.
    """
    df = pd.read_csv(csv_path)
    
    # Get column names
    # Expected: Read_ID, Ground_Truth, Mapper_Prediction, [ALIGNER]_Prediction, Mapper_Correct, [ALIGNER]_Correct
    ground_truth = df['Ground_Truth']
    mapper_pred = df['Mapper_Prediction']
    
    # Find direct prediction column (should be column index 3, but check by name pattern)
    direct_pred_col = None
    direct_correct_col = None
    for col in df.columns:
        if col.endswith('_Prediction') and col != 'Mapper_Prediction':
            direct_pred_col = col
        elif col.endswith('_Correct') and col != 'Mapper_Correct':
            direct_correct_col = col
    
    if direct_pred_col is None:
        # Fallback: use column index 3
        direct_pred = df.iloc[:, 3]
    else:
        direct_pred = df[direct_pred_col]
    
    # Calculate confusion matrix metrics
    all_labels = sorted(set(ground_truth) | set(mapper_pred) | set(direct_pred))
    
    # Mapper metrics
    mapper_prec, mapper_rec, _, mapper_support = precision_recall_fscore_support(
        ground_truth, mapper_pred, labels=all_labels, average='weighted', zero_division=0
    )
    
    # Direct metrics
    direct_prec, direct_rec, _, direct_support = precision_recall_fscore_support(
        ground_truth, direct_pred, labels=all_labels, average='weighted', zero_division=0
    )
    
    # Calculate specificity
    def calc_specificity(y_true, y_pred, labels):
        cm = confusion_matrix(y_true, y_pred, labels=labels)
        # Weighted specificity
        specificities = []
        supports = []
        for i, label in enumerate(labels):
            tn = np.sum(cm) - np.sum(cm[i, :]) - np.sum(cm[:, i]) + cm[i, i]
            fp = np.sum(cm[:, i]) - cm[i, i]
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
            specificities.append(spec)
            supports.append(np.sum(cm[i, :]))
        # Weighted average
        total_support = np.sum(supports)
        return np.sum(np.array(specificities) * np.array(supports)) / total_support if total_support > 0 else 0.0
    
    mapper_spec = calc_specificity(ground_truth, mapper_pred, all_labels)
    direct_spec = calc_specificity(ground_truth, direct_pred, all_labels)
    
    return {
        'mapper_precision': mapper_prec,
        'mapper_recall': mapper_rec,
        'mapper_specificity': mapper_spec,
        'direct_precision': direct_prec,
        'direct_recall': direct_rec,
        'direct_specificity': direct_spec
    }


def parse_variant_calling_summary(csv_path):
    """
    Parse variant calling summary CSV and extract overall metrics.
    Returns dict with precision, recall, specificity for combined and direct.
    """
    df = pd.read_csv(csv_path)
    
    # Get overall row
    overall = df[df['Gene'] == 'Overall'].iloc[0]
    
    return {
        'combined_precision': overall['Combined_Precision'],
        'combined_recall': overall['Combined_Recall'],
        'combined_specificity': overall['Combined_Specificity'],
        'direct_precision': overall['Direct_Precision'],
        'direct_recall': overall['Direct_Recall'],
        'direct_specificity': overall['Direct_Specificity']
    }


def aggregate_read_mapping_results(sim_output_base, aligners, seed_start, seed_end):
    """Aggregate read mapping results across all seeds."""
    results = {aligner: defaultdict(list) for aligner in aligners}
    
    for seed in range(seed_start, seed_end + 1):
        seed_dir = Path(sim_output_base) / f"seed_{seed}"
        
        for aligner in aligners:
            csv_path = seed_dir / aligner / "read_mapping_analysis" / f"seed_{seed}_{aligner}_summary.csv"
            
            if csv_path.exists():
                try:
                    metrics = parse_read_mapping_summary(csv_path)
                    results[aligner]['mapper_precision'].append(metrics['mapper_precision'])
                    results[aligner]['mapper_recall'].append(metrics['mapper_recall'])
                    results[aligner]['mapper_specificity'].append(metrics['mapper_specificity'])
                    results[aligner]['direct_precision'].append(metrics['direct_precision'])
                    results[aligner]['direct_recall'].append(metrics['direct_recall'])
                    results[aligner]['direct_specificity'].append(metrics['direct_specificity'])
                except Exception as e:
                    print(f"Warning: Failed to parse {csv_path}: {e}", file=sys.stderr)
                    continue
    
    # Calculate means and stds
    aggregated = {}
    for aligner in aligners:
        aggregated[aligner] = {
            'mapper_precision_mean': np.mean(results[aligner]['mapper_precision']),
            'mapper_precision_std': np.std(results[aligner]['mapper_precision']),
            'mapper_recall_mean': np.mean(results[aligner]['mapper_recall']),
            'mapper_recall_std': np.std(results[aligner]['mapper_recall']),
            'mapper_specificity_mean': np.mean(results[aligner]['mapper_specificity']),
            'mapper_specificity_std': np.std(results[aligner]['mapper_specificity']),
            'direct_precision_mean': np.mean(results[aligner]['direct_precision']),
            'direct_precision_std': np.std(results[aligner]['direct_precision']),
            'direct_recall_mean': np.mean(results[aligner]['direct_recall']),
            'direct_recall_std': np.std(results[aligner]['direct_recall']),
            'direct_specificity_mean': np.mean(results[aligner]['direct_specificity']),
            'direct_specificity_std': np.std(results[aligner]['direct_specificity']),
            'n_seeds': len(results[aligner]['mapper_precision'])
        }
    
    return aggregated


def aggregate_variant_calling_results(sim_output_base, aligners, seed_start, seed_end):
    """Aggregate variant calling results across all seeds."""
    results = {aligner: defaultdict(list) for aligner in aligners}
    
    for seed in range(seed_start, seed_end + 1):
        seed_dir = Path(sim_output_base) / f"seed_{seed}"
        
        for aligner in aligners:
            csv_path = seed_dir / aligner / "variant_calling_analysis" / f"paradism_seed_{seed}_{aligner}_summary.csv"
            
            if csv_path.exists():
                try:
                    metrics = parse_variant_calling_summary(csv_path)
                    results[aligner]['combined_precision'].append(metrics['combined_precision'])
                    results[aligner]['combined_recall'].append(metrics['combined_recall'])
                    results[aligner]['combined_specificity'].append(metrics['combined_specificity'])
                    results[aligner]['direct_precision'].append(metrics['direct_precision'])
                    results[aligner]['direct_recall'].append(metrics['direct_recall'])
                    results[aligner]['direct_specificity'].append(metrics['direct_specificity'])
                except Exception as e:
                    print(f"Warning: Failed to parse {csv_path}: {e}", file=sys.stderr)
                    continue
    
    # Calculate means and stds
    aggregated = {}
    for aligner in aligners:
        aggregated[aligner] = {
            'combined_precision_mean': np.mean(results[aligner]['combined_precision']),
            'combined_precision_std': np.std(results[aligner]['combined_precision']),
            'combined_recall_mean': np.mean(results[aligner]['combined_recall']),
            'combined_recall_std': np.std(results[aligner]['combined_recall']),
            'combined_specificity_mean': np.mean(results[aligner]['combined_specificity']),
            'combined_specificity_std': np.std(results[aligner]['combined_specificity']),
            'direct_precision_mean': np.mean(results[aligner]['direct_precision']),
            'direct_precision_std': np.std(results[aligner]['direct_precision']),
            'direct_recall_mean': np.mean(results[aligner]['direct_recall']),
            'direct_recall_std': np.std(results[aligner]['direct_recall']),
            'direct_specificity_mean': np.mean(results[aligner]['direct_specificity']),
            'direct_specificity_std': np.std(results[aligner]['direct_specificity']),
            'n_seeds': len(results[aligner]['combined_precision'])
        }
    
    return aggregated


def create_read_mapping_plot(aggregated, output_path):
    """Create read mapping comparison plot with error bars."""
    aligners = list(aggregated.keys())
    display_labels = aligners
    
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
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    y_pos = np.arange(len(display_labels))
    bar_height = 0.35
    
    # Extract data
    mapper_prec = [aggregated[a]['mapper_precision_mean'] for a in aligners]
    mapper_prec_err = [aggregated[a]['mapper_precision_std'] for a in aligners]
    direct_prec = [aggregated[a]['direct_precision_mean'] for a in aligners]
    direct_prec_err = [aggregated[a]['direct_precision_std'] for a in aligners]
    
    mapper_rec = [aggregated[a]['mapper_recall_mean'] for a in aligners]
    mapper_rec_err = [aggregated[a]['mapper_recall_std'] for a in aligners]
    direct_rec = [aggregated[a]['direct_recall_mean'] for a in aligners]
    direct_rec_err = [aggregated[a]['direct_recall_std'] for a in aligners]
    
    mapper_spec = [aggregated[a]['mapper_specificity_mean'] for a in aligners]
    mapper_spec_err = [aggregated[a]['mapper_specificity_std'] for a in aligners]
    direct_spec = [aggregated[a]['direct_specificity_mean'] for a in aligners]
    direct_spec_err = [aggregated[a]['direct_specificity_std'] for a in aligners]
    
    # Panel 1: Precision
    axes[0].barh(y_pos - bar_height/2, mapper_prec, bar_height, 
                 xerr=mapper_prec_err, label='Mapper', color='red', capsize=5)
    axes[0].barh(y_pos + bar_height/2, direct_prec, bar_height,
                 xerr=direct_prec_err, label='Direct', color='#4682B4', capsize=5)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels([a.upper() for a in display_labels], weight='bold')
    axes[0].set_xlabel('Precision', weight='bold')
    axes[0].set_title('Precision', weight='bold')
    axes[0].set_xlim([0, 1.05])
    axes[0].grid(axis='x', alpha=0.3)
    
    # Panel 2: Recall
    axes[1].barh(y_pos - bar_height/2, mapper_rec, bar_height,
                 xerr=mapper_rec_err, label='Mapper', color='red', capsize=5)
    axes[1].barh(y_pos + bar_height/2, direct_rec, bar_height,
                 xerr=direct_rec_err, label='Direct', color='#4682B4', capsize=5)
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels([a.upper() for a in display_labels], weight='bold')
    axes[1].set_xlabel('Recall', weight='bold')
    axes[1].set_title('Recall', weight='bold')
    axes[1].set_xlim([0, 1.05])
    axes[1].grid(axis='x', alpha=0.3)
    
    # Panel 3: Specificity
    axes[2].barh(y_pos - bar_height/2, mapper_spec, bar_height,
                 xerr=mapper_spec_err, label='Mapper', color='red', capsize=5)
    axes[2].barh(y_pos + bar_height/2, direct_spec, bar_height,
                 xerr=direct_spec_err, label='Direct', color='#4682B4', capsize=5)
    axes[2].set_yticks(y_pos)
    axes[2].set_yticklabels([a.upper() for a in display_labels], weight='bold')
    axes[2].set_xlabel('Specificity', weight='bold')
    axes[2].set_title('Specificity', weight='bold')
    axes[2].set_xlim([0, 1.05])
    axes[2].grid(axis='x', alpha=0.3)
    
    # Add legend
    handles, labels_legend = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_legend, loc='upper right', bbox_to_anchor=(1.02, 1.02),
               fontsize=16, frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_variant_calling_plot(aggregated, output_path):
    """Create variant calling comparison plot with error bars."""
    aligners = list(aggregated.keys())
    display_labels = aligners
    
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
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    y_pos = np.arange(len(display_labels))
    bar_height = 0.35
    
    # Extract data
    combined_prec = [aggregated[a]['combined_precision_mean'] for a in aligners]
    combined_prec_err = [aggregated[a]['combined_precision_std'] for a in aligners]
    direct_prec = [aggregated[a]['direct_precision_mean'] for a in aligners]
    direct_prec_err = [aggregated[a]['direct_precision_std'] for a in aligners]
    
    combined_rec = [aggregated[a]['combined_recall_mean'] for a in aligners]
    combined_rec_err = [aggregated[a]['combined_recall_std'] for a in aligners]
    direct_rec = [aggregated[a]['direct_recall_mean'] for a in aligners]
    direct_rec_err = [aggregated[a]['direct_recall_std'] for a in aligners]
    
    combined_spec = [aggregated[a]['combined_specificity_mean'] for a in aligners]
    combined_spec_err = [aggregated[a]['combined_specificity_std'] for a in aligners]
    direct_spec = [aggregated[a]['direct_specificity_mean'] for a in aligners]
    direct_spec_err = [aggregated[a]['direct_specificity_std'] for a in aligners]
    
    # Panel 1: Precision
    axes[0].barh(y_pos - bar_height/2, combined_prec, bar_height,
                 xerr=combined_prec_err, label='Combined', color='red', capsize=5)
    axes[0].barh(y_pos + bar_height/2, direct_prec, bar_height,
                 xerr=direct_prec_err, label='Direct', color='#4682B4', capsize=5)
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels([a.upper() for a in display_labels], weight='bold')
    axes[0].set_xlabel('Precision', weight='bold')
    axes[0].set_title('Precision', weight='bold')
    axes[0].set_xlim([0, 1.05])
    axes[0].grid(axis='x', alpha=0.3)
    
    # Panel 2: Recall
    axes[1].barh(y_pos - bar_height/2, combined_rec, bar_height,
                 xerr=combined_rec_err, label='Combined', color='red', capsize=5)
    axes[1].barh(y_pos + bar_height/2, direct_rec, bar_height,
                 xerr=direct_rec_err, label='Direct', color='#4682B4', capsize=5)
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels([a.upper() for a in display_labels], weight='bold')
    axes[1].set_xlabel('Recall', weight='bold')
    axes[1].set_title('Recall', weight='bold')
    axes[1].set_xlim([0, 1.05])
    axes[1].grid(axis='x', alpha=0.3)
    
    # Panel 3: Specificity
    axes[2].barh(y_pos - bar_height/2, combined_spec, bar_height,
                 xerr=combined_spec_err, label='Combined', color='red', capsize=5)
    axes[2].barh(y_pos + bar_height/2, direct_spec, bar_height,
                 xerr=direct_spec_err, label='Direct', color='#4682B4', capsize=5)
    axes[2].set_yticks(y_pos)
    axes[2].set_yticklabels([a.upper() for a in display_labels], weight='bold')
    axes[2].set_xlabel('Specificity', weight='bold')
    axes[2].set_title('Specificity', weight='bold')
    axes[2].set_xlim([0, 1.05])
    axes[2].grid(axis='x', alpha=0.3)
    
    # Add legend
    handles, labels_legend = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_legend, loc='upper right', bbox_to_anchor=(1.02, 1.02),
               fontsize=16, frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate results across all seeds and create summary plots'
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
            rm_data.append({
                'Aligner': aligner.upper(),
                'Mapper_Precision_Mean': read_mapping_agg[aligner]['mapper_precision_mean'],
                'Mapper_Precision_Std': read_mapping_agg[aligner]['mapper_precision_std'],
                'Mapper_Recall_Mean': read_mapping_agg[aligner]['mapper_recall_mean'],
                'Mapper_Recall_Std': read_mapping_agg[aligner]['mapper_recall_std'],
                'Mapper_Specificity_Mean': read_mapping_agg[aligner]['mapper_specificity_mean'],
                'Mapper_Specificity_Std': read_mapping_agg[aligner]['mapper_specificity_std'],
                'Direct_Precision_Mean': read_mapping_agg[aligner]['direct_precision_mean'],
                'Direct_Precision_Std': read_mapping_agg[aligner]['direct_precision_std'],
                'Direct_Recall_Mean': read_mapping_agg[aligner]['direct_recall_mean'],
                'Direct_Recall_Std': read_mapping_agg[aligner]['direct_recall_std'],
                'Direct_Specificity_Mean': read_mapping_agg[aligner]['direct_specificity_mean'],
                'Direct_Specificity_Std': read_mapping_agg[aligner]['direct_specificity_std'],
                'N_Seeds': read_mapping_agg[aligner]['n_seeds']
            })
    rm_df = pd.DataFrame(rm_data)
    rm_csv = output_dir / 'read_mapping_aggregated_summary.csv'
    rm_df.to_csv(rm_csv, index=False)
    print(f"  Saved: {rm_csv}")
    
    # Variant calling CSV
    vc_data = []
    for aligner in args.aligners:
        if aligner in variant_calling_agg:
            vc_data.append({
                'Aligner': aligner.upper(),
                'Combined_Precision_Mean': variant_calling_agg[aligner]['combined_precision_mean'],
                'Combined_Precision_Std': variant_calling_agg[aligner]['combined_precision_std'],
                'Combined_Recall_Mean': variant_calling_agg[aligner]['combined_recall_mean'],
                'Combined_Recall_Std': variant_calling_agg[aligner]['combined_recall_std'],
                'Combined_Specificity_Mean': variant_calling_agg[aligner]['combined_specificity_mean'],
                'Combined_Specificity_Std': variant_calling_agg[aligner]['combined_specificity_std'],
                'Direct_Precision_Mean': variant_calling_agg[aligner]['direct_precision_mean'],
                'Direct_Precision_Std': variant_calling_agg[aligner]['direct_precision_std'],
                'Direct_Recall_Mean': variant_calling_agg[aligner]['direct_recall_mean'],
                'Direct_Recall_Std': variant_calling_agg[aligner]['direct_recall_std'],
                'Direct_Specificity_Mean': variant_calling_agg[aligner]['direct_specificity_mean'],
                'Direct_Specificity_Std': variant_calling_agg[aligner]['direct_specificity_std'],
                'N_Seeds': variant_calling_agg[aligner]['n_seeds']
            })
    vc_df = pd.DataFrame(vc_data)
    vc_csv = output_dir / 'variant_calling_aggregated_summary.csv'
    vc_df.to_csv(vc_csv, index=False)
    print(f"  Saved: {vc_csv}")
    
    # Create plots
    print("Creating plots...")
    
    if read_mapping_agg:
        rm_plot = output_dir / 'read_mapping_aggregated_comparison.png'
        create_read_mapping_plot(read_mapping_agg, rm_plot)
        print(f"  Saved: {rm_plot}")
    
    if variant_calling_agg:
        vc_plot = output_dir / 'variant_calling_aggregated_comparison.png'
        create_variant_calling_plot(variant_calling_agg, vc_plot)
        print(f"  Saved: {vc_plot}")
    
    print(f"\nAggregation complete! Results saved to: {output_dir}/")
    print(f"  Read mapping: {len([a for a in args.aligners if a in read_mapping_agg])} aligners")
    print(f"  Variant calling: {len([a for a in args.aligners if a in variant_calling_agg])} aligners")


if __name__ == "__main__":
    main()

