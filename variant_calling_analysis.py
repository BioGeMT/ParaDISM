#!/usr/bin/env python3
"""
Analyze variant calling results with confusion matrices.
Compares combined VCF and direct VCF against ground truth.
"""

from typing import Any


import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support, accuracy_score
import matplotlib.pyplot as plt


def parse_vcf(vcf_path, filter_snps_only=False):
    """
    Parse VCF file and return set of variants.
    Variants are stored as tuples: (CHROM, POS, REF, ALT)
    """
    variants = set[Any]()
    
    if not Path(vcf_path).exists():
        return variants
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            # Filter SNPs only if requested
            if filter_snps_only:
                # Check if it's a SNP (single nucleotide substitution)
                if len(ref) != 1 or len(alt) != 1:
                    continue
                # Also check INFO field if TYPE is present
                if len(parts) > 7:
                    info = parts[7]
                    if 'TYPE=' in info:
                        if 'TYPE=SNP' not in info and 'TYPE=snp' not in info:
                            continue
            
            variants.add((chrom, pos, ref, alt))
    
    return variants


def parse_all_positions_from_tsv(tsv_path):
    """
    Parse TSV file to get all positions (GENE, POS) and which are variants.
    Returns: (all_positions_set, variant_positions_set)
    """
    all_positions = set()
    variant_positions = set()
    
    if not Path(tsv_path).exists():
        return all_positions, variant_positions
    
    with open(tsv_path, 'r') as f:
        header = f.readline()  # Skip header
        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:  # Skip empty lines
                continue
            
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            
            gene = parts[0].strip()
            pos_str = parts[1].strip()
            mut_type = parts[4].strip()
            
            # Skip if gene or position is empty
            if not gene or not pos_str:
                continue
            
            try:
                pos = int(pos_str)
            except ValueError:
                # Skip invalid position lines
                continue
            
            position = (gene, pos)
            all_positions.add(position)
            
            # If it's not MATCH, it's a variant position
            if mut_type != 'MATCH':
                variant_positions.add(position)
    
    return all_positions, variant_positions


def parse_vcf_with_info(vcf_path):
    """
    Parse VCF file and return dict with variant info.
    Key: (CHROM, POS, REF, ALT)
    Value: dict with INFO fields
    """
    variants = {}
    
    if not Path(vcf_path).exists():
        return variants
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            info = {}
            if len(parts) > 7:
                info_str = parts[7]
                for item in info_str.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info[key] = value
            
            variants[(chrom, pos, ref, alt)] = info
    
    return variants


def get_all_positions(variants):
    """Get all unique (CHROM, POS) positions from variants."""
    positions = set()
    for chrom, pos, ref, alt in variants:
        positions.add((chrom, pos))
    return positions


def compare_variants(ground_truth, predicted):
    """
    Compare predicted variants against ground truth.
    Returns TP, FP, FN sets and metrics.
    """
    tp = ground_truth & predicted  # True Positives: in both
    fp = predicted - ground_truth  # False Positives: predicted but not in truth
    fn = ground_truth - predicted  # False Negatives: in truth but not predicted
    
    return tp, fp, fn


def compare_variants_with_positions(ground_truth_variants, predicted_variants, all_positions):
    """
    Compare predicted variants against ground truth with complete position information.
    Returns TP, FP, FN, TN sets.
    """
    # Get all positions from predicted variants
    predicted_positions = {(chrom, pos) for chrom, pos, ref, alt in predicted_variants}
    ground_truth_positions = {(chrom, pos) for chrom, pos, ref, alt in ground_truth_variants}
    
    # Variant-level comparison
    tp = ground_truth_variants & predicted_variants  # True Positives: exact match
    fp = predicted_variants - ground_truth_variants  # False Positives: predicted but not in truth
    fn = ground_truth_variants - predicted_variants  # False Negatives: in truth but not predicted
    
    # Calculate True Negatives: positions with no variant in truth AND no variant predicted
    # These are positions in all_positions that are NOT in ground_truth_positions AND NOT in predicted_positions
    tn_positions = all_positions - ground_truth_positions - predicted_positions
    
    # Count TN: positions correctly identified as non-variant
    tn_count = len(tn_positions)
    
    # Also count FP positions: positions called as variant but shouldn't be
    fp_positions = predicted_positions - ground_truth_positions
    fp_position_count = len(fp_positions)
    
    return tp, fp, fn, tn_count, fp_position_count


def compare_variants(ground_truth, predicted):
    """
    Compare predicted variants against ground truth (legacy method for backward compatibility).
    Returns TP, FP, FN sets.
    """
    tp = ground_truth & predicted  # True Positives: in both
    fp = predicted - ground_truth  # False Positives: predicted but not in truth
    fn = ground_truth - predicted  # False Negatives: in truth but not predicted
    
    return tp, fp, fn


def calculate_metrics(tp, fp, fn, tn=0, fp_positions=0):
    """Calculate precision, recall, F1 score, and specificity."""
    tp_count = len(tp)
    fp_count = len(fp)
    fn_count = len(fn)
    
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    
    # Specificity = TN / (TN + FP)
    # If we have true TN count (from position information), use it
    # Otherwise use approximation
    if tn > 0:
        # True specificity using position-level information
        # FP positions = positions incorrectly called as variant
        specificity = tn / (tn + fp_positions) if (tn + fp_positions) > 0 else 0.0
    else:
        # Fallback approximation when we don't have all positions
        # TN approximation: variants we correctly didn't call as false positives
        tn_estimate = tp_count + fn_count
        specificity = tn_estimate / (tn_estimate + fp_count) if (tn_estimate + fp_count) > 0 else 0.0
    
    return {
        'tp': tp_count,
        'fp': fp_count,
        'fn': fn_count,
        'tn': tn,
        'precision': precision,
        'recall': recall,
        'specificity': specificity,
        'f1': f1
    }


def calculate_per_gene_metrics_legacy(ground_truth, combined, direct):
    """Calculate metrics per gene without position information (legacy method)."""
    # Get all genes
    genes = set()
    for chrom, pos, ref, alt in ground_truth | combined | direct:
        genes.add(chrom)
    
    gene_metrics = {}
    for gene in sorted(genes):
        # Filter variants for this gene
        gt_gene = {v for v in ground_truth if v[0] == gene}
        combined_gene = {v for v in combined if v[0] == gene}
        direct_gene = {v for v in direct if v[0] == gene}
        
        tp_combined, fp_combined, fn_combined = compare_variants(gt_gene, combined_gene)
        tp_direct, fp_direct, fn_direct = compare_variants(gt_gene, direct_gene)
        
        metrics_combined = calculate_metrics(tp_combined, fp_combined, fn_combined)
        metrics_direct = calculate_metrics(tp_direct, fp_direct, fn_direct)
        
        gene_metrics[gene] = {
            'combined': metrics_combined,
            'direct': metrics_direct
        }
    
    return gene_metrics


def calculate_per_gene_metrics(ground_truth, combined, direct, all_positions):
    """Calculate metrics per gene using position information."""
    # Get all genes
    genes = set()
    for chrom, pos, ref, alt in ground_truth | combined | direct:
        genes.add(chrom)
    
    gene_metrics = {}
    for gene in sorted(genes):
        # Filter variants for this gene
        gt_gene = {v for v in ground_truth if v[0] == gene}
        combined_gene = {v for v in combined if v[0] == gene}
        direct_gene = {v for v in direct if v[0] == gene}
        
        # Filter positions for this gene
        gene_positions = {(g, p) for g, p in all_positions if g == gene}
        
        # Use position-aware comparison
        tp_combined, fp_combined, fn_combined, tn_combined, fp_pos_combined = compare_variants_with_positions(
            gt_gene, combined_gene, gene_positions
        )
        tp_direct, fp_direct, fn_direct, tn_direct, fp_pos_direct = compare_variants_with_positions(
            gt_gene, direct_gene, gene_positions
        )
        
        metrics_combined = calculate_metrics(tp_combined, fp_combined, fn_combined, tn_combined, fp_pos_combined)
        metrics_direct = calculate_metrics(tp_direct, fp_direct, fn_direct, tn_direct, fp_pos_direct)
        
        gene_metrics[gene] = {
            'combined': metrics_combined,
            'direct': metrics_direct
        }
    
    return gene_metrics


def create_grouped_bar_chart(gene_metrics, aligner_name, output_path):
    """
    Create a 3-panel horizontal grouped bar chart showing precision, recall, and specificity.
    Each panel has genes + overall, each with 2 bars (combined vs direct).
    """
    genes = sorted([g for g in gene_metrics.keys()])
    display_labels = genes + ['Overall']
    
    # Calculate overall metrics
    all_combined = {'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0, 'fp_positions': 0}
    all_direct = {'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0, 'fp_positions': 0}
    
    for gene, metrics in gene_metrics.items():
        all_combined['tp'] += metrics['combined']['tp']
        all_combined['fp'] += metrics['combined']['fp']
        all_combined['fn'] += metrics['combined']['fn']
        # Get TN and FP positions from metrics if available
        if 'tn' in metrics['combined']:
            all_combined['tn'] += metrics['combined']['tn']
        # For FP positions, we need to track them separately
        # Since we're aggregating, we'll use FP count as proxy for FP positions
        # (each FP variant represents a position incorrectly called)
        all_direct['tp'] += metrics['direct']['tp']
        all_direct['fp'] += metrics['direct']['fp']
        all_direct['fn'] += metrics['direct']['fn']
        if 'tn' in metrics['direct']:
            all_direct['tn'] += metrics['direct']['tn']
    
    # Calculate overall metrics using the actual counts
    if all_combined['tn'] > 0:
        # Use true TN values from position information
        # For FP positions in specificity calculation, use FP count as proxy
        # (each false positive variant = one position incorrectly called)
        overall_combined = calculate_metrics(
            set(), set(), set(), 
            tn=all_combined['tn'], 
            fp_positions=all_combined['fp']  # FP variant count = FP positions
        )
        overall_combined['tp'] = all_combined['tp']
        overall_combined['fp'] = all_combined['fp']
        overall_combined['fn'] = all_combined['fn']
        # Recalculate precision and recall with actual counts
        overall_combined['precision'] = all_combined['tp'] / (all_combined['tp'] + all_combined['fp']) if (all_combined['tp'] + all_combined['fp']) > 0 else 0.0
        overall_combined['recall'] = all_combined['tp'] / (all_combined['tp'] + all_combined['fn']) if (all_combined['tp'] + all_combined['fn']) > 0 else 0.0
        
        overall_direct = calculate_metrics(
            set(), set(), set(),
            tn=all_direct['tn'],
            fp_positions=all_direct['fp']  # FP variant count = FP positions
        )
        overall_direct['tp'] = all_direct['tp']
        overall_direct['fp'] = all_direct['fp']
        overall_direct['fn'] = all_direct['fn']
        # Recalculate precision and recall with actual counts
        overall_direct['precision'] = all_direct['tp'] / (all_direct['tp'] + all_direct['fp']) if (all_direct['tp'] + all_direct['fp']) > 0 else 0.0
        overall_direct['recall'] = all_direct['tp'] / (all_direct['tp'] + all_direct['fn']) if (all_direct['tp'] + all_direct['fn']) > 0 else 0.0
    else:
        # Fallback to approximation
        overall_combined = calculate_metrics(set(), set(), set())
        overall_combined['tp'] = all_combined['tp']
        overall_combined['fp'] = all_combined['fp']
        overall_combined['fn'] = all_combined['fn']
        overall_combined['precision'] = all_combined['tp'] / (all_combined['tp'] + all_combined['fp']) if (all_combined['tp'] + all_combined['fp']) > 0 else 0.0
        overall_combined['recall'] = all_combined['tp'] / (all_combined['tp'] + all_combined['fn']) if (all_combined['tp'] + all_combined['fn']) > 0 else 0.0
        tn_combined = all_combined['tp'] + all_combined['fn']
        overall_combined['specificity'] = tn_combined / (tn_combined + all_combined['fp']) if (tn_combined + all_combined['fp']) > 0 else 0.0
        
        overall_direct = calculate_metrics(set(), set(), set())
        overall_direct['tp'] = all_direct['tp']
        overall_direct['fp'] = all_direct['fp']
        overall_direct['fn'] = all_direct['fn']
        overall_direct['precision'] = all_direct['tp'] / (all_direct['tp'] + all_direct['fp']) if (all_direct['tp'] + all_direct['fp']) > 0 else 0.0
        overall_direct['recall'] = all_direct['tp'] / (all_direct['tp'] + all_direct['fn']) if (all_direct['tp'] + all_direct['fn']) > 0 else 0.0
        tn_direct = all_direct['tp'] + all_direct['fn']
        overall_direct['specificity'] = tn_direct / (tn_direct + all_direct['fp']) if (tn_direct + all_direct['fp']) > 0 else 0.0
    
    # Extract metrics for genes
    combined_precision = [gene_metrics[g]['combined']['precision'] for g in genes]
    combined_recall = [gene_metrics[g]['combined']['recall'] for g in genes]
    combined_specificity = [gene_metrics[g]['combined']['specificity'] for g in genes]
    
    direct_precision = [gene_metrics[g]['direct']['precision'] for g in genes]
    direct_recall = [gene_metrics[g]['direct']['recall'] for g in genes]
    direct_specificity = [gene_metrics[g]['direct']['specificity'] for g in genes]
    
    # Add overall metrics
    combined_precision.append(overall_combined['precision'])
    combined_recall.append(overall_combined['recall'])
    combined_specificity.append(overall_combined['specificity'])
    
    direct_precision.append(overall_direct['precision'])
    direct_recall.append(overall_direct['recall'])
    direct_specificity.append(overall_direct['specificity'])
    
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
    
    # Panel 1: Precision
    axes[0].barh(y_pos - bar_height/2, combined_precision, bar_height, label='paradism', color='red')
    axes[0].barh(y_pos + bar_height/2, direct_precision, bar_height, label=aligner_name, color='#4682B4')
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(display_labels, weight='bold')
    axes[0].set_xlabel('Precision', weight='bold')
    axes[0].set_title('Precision', weight='bold')
    axes[0].set_xlim([0, 1.05])
    axes[0].grid(axis='x', alpha=0.3)
    
    # Panel 2: Recall
    axes[1].barh(y_pos - bar_height/2, combined_recall, bar_height, label='paradism', color='red')
    axes[1].barh(y_pos + bar_height/2, direct_recall, bar_height, label=aligner_name, color='#4682B4')
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels(display_labels, weight='bold')
    axes[1].set_xlabel('Recall', weight='bold')
    axes[1].set_title('Recall', weight='bold')
    axes[1].set_xlim([0, 1.05])
    axes[1].grid(axis='x', alpha=0.3)
    
    # Panel 3: Specificity
    axes[2].barh(y_pos - bar_height/2, combined_specificity, bar_height, label='paradism', color='red')
    axes[2].barh(y_pos + bar_height/2, direct_specificity, bar_height, label=aligner_name, color='#4682B4')
    axes[2].set_yticks(y_pos)
    axes[2].set_yticklabels(display_labels, weight='bold')
    axes[2].set_xlabel('Specificity', weight='bold')
    axes[2].set_title('Specificity', weight='bold')
    axes[2].set_xlim([0, 1.05])
    axes[2].grid(axis='x', alpha=0.3)
    
    # Add single legend at top right (moved further away)
    handles, labels_legend = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels_legend, loc='upper right', bbox_to_anchor=(1.02, 1.02),
               fontsize=16, frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def print_summary(metrics_combined, metrics_direct, aligner_name):
    """Print concise summary."""
    print(f"\nParaDISM:         Prec={metrics_combined['precision']:.4f}  Rec={metrics_combined['recall']:.4f}  Spec={metrics_combined['specificity']:.4f}")
    print(f"{aligner_name.upper():<14}Prec={metrics_direct['precision']:.4f}  Rec={metrics_direct['recall']:.4f}  Spec={metrics_direct['specificity']:.4f}")
    print(f"\nTrue Positives:   ParaDISM={metrics_combined['tp']}, Direct={metrics_direct['tp']}")
    print(f"False Positives:  ParaDISM={metrics_combined['fp']}, Direct={metrics_direct['fp']}")
    print(f"False Negatives:  ParaDISM={metrics_combined['fn']}, Direct={metrics_direct['fn']}")
    if metrics_combined.get('tn', 0) > 0:
        print(f"True Negatives:   ParaDISM={metrics_combined['tn']}, Direct={metrics_direct.get('tn', 0)}")
        print(f"  (Using position-level specificity calculation)")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze variant calling results with confusion matrices'
    )
    parser.add_argument(
        '--aligner',
        default='bwa-mem2',
        choices=['bowtie2', 'bwa-mem2', 'minimap2'],
        help='Aligner to analyze (default: bwa-mem2)'
    )
    parser.add_argument(
        '--ground-truth',
        default='all_genes_mutations.vcf',
        help='Ground truth VCF file (default: all_genes_mutations.vcf)'
    )
    parser.add_argument(
        '--combined-vcf',
        default=None,
        help='Combined VCF file path (default: sim_freebayes_out/{aligner}/{aligner}_combined.vcf)'
    )
    parser.add_argument(
        '--direct-vcf',
        default=None,
        help='Direct VCF file path (default: sim_freebayes_out/{aligner}/{aligner}_direct.vcf)'
    )
    parser.add_argument(
        '--output-prefix',
        default=None,
        help='Output prefix for plots and reports (default: metrics_{aligner})'
    )
    parser.add_argument(
        '--analysis-dir',
        default=None,
        help='Directory to store analysis outputs (default: sim_variant_calling_analysis/{aligner})'
    )
    parser.add_argument(
        '--positions-tsv',
        default='all_genes_mutations.tsv',
        help='TSV file with all positions (GENE, POS, TYPE) for true specificity calculation (default: all_genes_mutations.tsv)'
    )
    
    args = parser.parse_args()
    
    # Set default paths
    if args.combined_vcf is None:
        args.combined_vcf = f'sim_freebayes_out/{args.aligner}/{args.aligner}_combined.vcf'
    if args.direct_vcf is None:
        args.direct_vcf = f'sim_freebayes_out/{args.aligner}/{args.aligner}_direct.vcf'
    
    # Set default output prefix
    if args.output_prefix is None:
        args.output_prefix = f'metrics_{args.aligner}'
    
    # Create output directory
    if args.analysis_dir is None:
        analysis_dir = Path('sim_variant_calling_analysis') / args.aligner
    else:
        analysis_dir = Path(args.analysis_dir)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # Update output prefix to include directory
    args.output_prefix = str(analysis_dir / args.output_prefix)
    
    # Check required files exist
    for filepath in [args.ground_truth, args.combined_vcf, args.direct_vcf]:
        if not Path(filepath).exists():
            print(f"Error: File not found: {filepath}")
            sys.exit(1)
    
    # Parse all positions from TSV for true specificity calculation
    print("Parsing all positions from TSV file...")
    all_positions, variant_positions = parse_all_positions_from_tsv(args.positions_tsv)
    if len(all_positions) > 0:
        print(f"  Found {len(all_positions)} total positions")
        print(f"  Found {len(variant_positions)} variant positions")
        use_position_info = True
    else:
        print("  Warning: TSV file not found or empty, using approximation for specificity")
        use_position_info = False
        all_positions = set()
    
    print("Parsing ground truth VCF (SNPs only)...")
    ground_truth = parse_vcf(args.ground_truth, filter_snps_only=True)
    print(f"  Found {len(ground_truth)} SNPs in ground truth")
    
    print("Parsing ParaDISM combined VCF...")
    combined = parse_vcf(args.combined_vcf)
    print(f"  Found {len(combined)} variants in ParaDISM combined VCF")
    
    print("Parsing direct VCF...")
    direct = parse_vcf(args.direct_vcf)
    print(f"  Found {len(direct)} variants in direct VCF")
    
    # Calculate overall metrics
    if use_position_info:
        tp_combined, fp_combined, fn_combined, tn_combined, fp_pos_combined = compare_variants_with_positions(
            ground_truth, combined, all_positions
        )
        tp_direct, fp_direct, fn_direct, tn_direct, fp_pos_direct = compare_variants_with_positions(
            ground_truth, direct, all_positions
        )
        metrics_combined = calculate_metrics(tp_combined, fp_combined, fn_combined, tn_combined, fp_pos_combined)
        metrics_direct = calculate_metrics(tp_direct, fp_direct, fn_direct, tn_direct, fp_pos_direct)
    else:
        tp_combined, fp_combined, fn_combined = compare_variants(ground_truth, combined)
        tp_direct, fp_direct, fn_direct = compare_variants(ground_truth, direct)
        metrics_combined = calculate_metrics(tp_combined, fp_combined, fn_combined)
        metrics_direct = calculate_metrics(tp_direct, fp_direct, fn_direct)
    
    # Calculate per-gene metrics
    print("\nCalculating per-gene metrics...")
    if use_position_info:
        gene_metrics = calculate_per_gene_metrics(ground_truth, combined, direct, all_positions)
    else:
        # Fallback: use old method without position info
        gene_metrics = calculate_per_gene_metrics_legacy(ground_truth, combined, direct)
    
    # Create grouped bar chart
    print("Creating comparison plot...")
    create_grouped_bar_chart(
        gene_metrics, args.aligner,
        f"{args.output_prefix}_comparison.png"
    )
    
    # Save summary to CSV
    summary_data = []
    for gene in sorted(gene_metrics.keys()):
        metrics = gene_metrics[gene]
        summary_data.append({
            'Gene': gene,
            'Combined_TP': metrics['combined']['tp'],
            'Combined_FP': metrics['combined']['fp'],
            'Combined_FN': metrics['combined']['fn'],
            'Combined_Precision': metrics['combined']['precision'],
            'Combined_Recall': metrics['combined']['recall'],
            'Combined_Specificity': metrics['combined']['specificity'],
            'Direct_TP': metrics['direct']['tp'],
            'Direct_FP': metrics['direct']['fp'],
            'Direct_FN': metrics['direct']['fn'],
            'Direct_Precision': metrics['direct']['precision'],
            'Direct_Recall': metrics['direct']['recall'],
            'Direct_Specificity': metrics['direct']['specificity'],
        })
    
    # Add overall row
    summary_data.append({
        'Gene': 'Overall',
        'Combined_TP': metrics_combined['tp'],
        'Combined_FP': metrics_combined['fp'],
        'Combined_FN': metrics_combined['fn'],
        'Combined_Precision': metrics_combined['precision'],
        'Combined_Recall': metrics_combined['recall'],
        'Combined_Specificity': metrics_combined['specificity'],
        'Direct_TP': metrics_direct['tp'],
        'Direct_FP': metrics_direct['fp'],
        'Direct_FN': metrics_direct['fn'],
        'Direct_Precision': metrics_direct['precision'],
        'Direct_Recall': metrics_direct['recall'],
        'Direct_Specificity': metrics_direct['specificity'],
    })
    
    df = pd.DataFrame(summary_data)
    csv_path = f"{args.output_prefix}_summary.csv"
    df.to_csv(csv_path, index=False)
    
    # Print summary
    print_summary(metrics_combined, metrics_direct, args.aligner)
    print(f"\nOutputs saved to: sim_variant_calling_analysis/{args.aligner}/")
    print(f"  Plot: {Path(args.output_prefix).name}_comparison.png")
    print(f"  CSV:  {Path(csv_path).name}")


if __name__ == "__main__":
    main()
