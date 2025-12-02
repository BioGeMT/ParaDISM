#!/usr/bin/env python3
"""
Analyze simulated read mapping results with confusion matrices.
Compares mapper predictions and BWA direct alignment against ground truth.
"""

import argparse
import sys
from pathlib import Path
import pysam
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support, accuracy_score
import matplotlib.pyplot as plt


def parse_ground_truth_from_fastq(fastq_path):
    """Extract ground truth gene for each read from FASTQ headers.

    Supports two formats:
    1) Original simulator:
         @Read0 PKD1P3; 9509-9899
       -> read_id = "Read0", gene = "PKD1P3"

    2) DWGSIM-style (gene embedded in read name, no separate field):
         @PKD1_20058_20275_0_1_0_0_2:3:0_2:1:0_0/1
       -> read_id = "PKD1_20058_20275_0_1_0_0_2:3:0_2:1:0_0"
          gene    = "PKD1"
    """
    ground_truth = {}

    with open(fastq_path, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            _ = f.readline()  # seq
            _ = f.readline()  # plus
            _ = f.readline()  # qual

            parts = header.strip().split()
            if not parts:
                continue

            # Strip '@' and any /1 or /2 mate suffix from read ID
            read_id = parts[0][1:]
            if "/" in read_id:
                read_id = read_id.split("/", 1)[0]

            gene = None

            # Case 1: original simulator format with explicit gene field
            if len(parts) >= 2:
                gene = parts[1].rstrip(';')
            else:
                # Case 2: DWGSIM-style name; gene is the prefix before first '_'
                # Example: PKD1_20058_20275_0_1_0_0_2:... -> PKD1
                if "_" in read_id:
                    gene = read_id.split("_", 1)[0]

            if gene is not None:
                ground_truth[read_id] = gene

    return ground_truth


def parse_mapper_results(tsv_path):
    """Parse mapper unique_mappings.tsv output."""
    mapper_predictions = {}
    partial_reads = set()

    with open(tsv_path, 'r') as f:
        _ = f.readline()  # Skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue

            read_name = parts[0]
            gene = parts[1] if parts[1] != "NONE" else "NONE"
            
            # Check for unpaired read indicators (both old and new format)
            is_partial = False
            if gene.endswith("(1/2 read mates)"):
                partial_reads.add(read_name)
                is_partial = True
                # Strip the old format tag
                gene = gene.replace(" (1/2 read mates)", "")
            elif " (plus)" in gene or " (minus)" in gene:
                partial_reads.add(read_name)
                is_partial = True
                # Strip the new format tag (plus/minus)
                gene = gene.replace(" (plus)", "").replace(" (minus)", "")

            # Remove strand suffix (+/-) if present from read name
            read_name = read_name.rstrip('+-')
            # Store gene without strand tags for comparison
            mapper_predictions[read_name] = gene

    return mapper_predictions, partial_reads


def parse_bwa_sam(sam_path):
    """Parse BWA SAM file to extract primary alignments."""
    bwa_predictions = {}

    mode = "rb" if sam_path.endswith(('.bam', '.cram')) else "r"

    with pysam.AlignmentFile(sam_path, mode) as samfile:
        for read in samfile.fetch(until_eof=True):
            # Skip secondary and supplementary alignments
            if read.is_secondary or read.is_supplementary:
                continue

            read_name = read.query_name

            if read.is_unmapped:
                gene = "NONE"
            else:
                # Reference name is the gene
                gene = samfile.get_reference_name(read.reference_id)

            # For paired-end, we may see both mates - just take first
            if read_name not in bwa_predictions:
                bwa_predictions[read_name] = gene

    return bwa_predictions


def calculate_specificity(y_true, y_pred, labels):
    """
    Calculate specificity for each class.
    Specificity = TN / (TN + FP)
    """
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    specificities = []

    for i, label in enumerate(labels):
        # True Negatives: all samples not in this class that were correctly not predicted as this class
        # False Positives: all samples not in this class that were incorrectly predicted as this class

        # Sum of all cells except row i and column i
        tn = np.sum(cm) - np.sum(cm[i, :]) - np.sum(cm[:, i]) + cm[i, i]

        # Sum of column i except row i (false positives)
        fp = np.sum(cm[:, i]) - cm[i, i]

        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        specificities.append(specificity)

    return np.array(specificities)


def calculate_metrics(y_true, y_pred, labels):
    """Calculate precision, recall, and specificity for each class and overall."""
    # Per-class metrics
    precision, recall, _, support = precision_recall_fscore_support(
        y_true, y_pred, labels=labels, average=None, zero_division=0
    )
    specificity = calculate_specificity(y_true, y_pred, labels)

    # Overall weighted averages (weighted by support)
    precision_overall, recall_overall, _, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=labels, average='weighted', zero_division=0
    )

    # Calculate overall specificity as weighted average
    total_support = np.sum(support)
    specificity_overall = np.sum(specificity * support) / total_support if total_support > 0 else 0.0

    return {
        'precision': precision,
        'recall': recall,
        'specificity': specificity,
        'precision_overall': precision_overall,
        'recall_overall': recall_overall,
        'specificity_overall': specificity_overall,
        'support': support
    }


def create_grouped_bar_chart(metrics_mapper, metrics_direct, labels, aligner_name, output_path):
    """
    Create a 3-panel horizontal grouped bar chart showing precision, recall, and specificity.
    Each panel has 8 groups (7 genes + overall), each with 2 bars (mapper vs direct aligner).
    """
    # Prepare data
    genes = [label for label in labels if label != 'NONE']
    display_labels = genes + ['Overall']

    # Extract metrics for genes (excluding NONE)
    mapper_precision = [metrics_mapper['precision'][i] for i, label in enumerate(labels) if label != 'NONE']
    mapper_recall = [metrics_mapper['recall'][i] for i, label in enumerate(labels) if label != 'NONE']
    mapper_specificity = [metrics_mapper['specificity'][i] for i, label in enumerate(labels) if label != 'NONE']

    direct_precision = [metrics_direct['precision'][i] for i, label in enumerate(labels) if label != 'NONE']
    direct_recall = [metrics_direct['recall'][i] for i, label in enumerate(labels) if label != 'NONE']
    direct_specificity = [metrics_direct['specificity'][i] for i, label in enumerate(labels) if label != 'NONE']

    # Add overall metrics
    mapper_precision.append(metrics_mapper['precision_overall'])
    mapper_recall.append(metrics_mapper['recall_overall'])
    mapper_specificity.append(metrics_mapper['specificity_overall'])

    direct_precision.append(metrics_direct['precision_overall'])
    direct_recall.append(metrics_direct['recall_overall'])
    direct_specificity.append(metrics_direct['specificity_overall'])

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
    axes[0].barh(y_pos - bar_height/2, mapper_precision, bar_height, label='paradism', color='red')
    axes[0].barh(y_pos + bar_height/2, direct_precision, bar_height, label=aligner_name, color='#4682B4')
    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(display_labels, weight='bold')
    axes[0].set_xlabel('Precision', weight='bold')
    axes[0].set_title('Precision', weight='bold')
    axes[0].set_xlim([0, 1.05])
    axes[0].grid(axis='x', alpha=0.3)
    
    # Panel 2: Recall
    axes[1].barh(y_pos - bar_height/2, mapper_recall, bar_height, label='paradism', color='red')
    axes[1].barh(y_pos + bar_height/2, direct_recall, bar_height, label=aligner_name, color='#4682B4')
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels(display_labels, weight='bold')
    axes[1].set_xlabel('Recall', weight='bold')
    axes[1].set_title('Recall', weight='bold')
    axes[1].set_xlim([0, 1.05])
    axes[1].grid(axis='x', alpha=0.3)
    
    # Panel 3: Specificity
    axes[2].barh(y_pos - bar_height/2, mapper_specificity, bar_height, label='paradism', color='red')
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


def print_summary(metrics_mapper, metrics_direct, aligner_name):
    """Print concise summary."""
    print(f"\nMapper:           Prec={metrics_mapper['precision_overall']:.4f}  Rec={metrics_mapper['recall_overall']:.4f}  Spec={metrics_mapper['specificity_overall']:.4f}")
    print(f"{aligner_name.upper():<14}Prec={metrics_direct['precision_overall']:.4f}  Rec={metrics_direct['recall_overall']:.4f}  Spec={metrics_direct['specificity_overall']:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze simulated read mapping results with confusion matrices'
    )
    parser.add_argument(
        '--aligner',
        default='bwa-mem2',
        choices=['bowtie2', 'bwa-mem2', 'minimap2'],
        help='Aligner to analyze (default: bwa-mem2)'
    )
    parser.add_argument(
        '--fastq-r1',
        default='simulated_r1_err_0_010.fq',
        help='Simulated R1 FASTQ file (contains ground truth)'
    )
    parser.add_argument(
        '--mapper-tsv',
        default=None,
        help='Path to ParaDISM unique_mappings TSV (default based on aligner)'
    )
    parser.add_argument(
        '--direct-sam',
        default=None,
        help='Path to baseline SAM/BAM (default based on aligner)'
    )
    parser.add_argument(
        '--analysis-dir',
        default=None,
        help='Directory to store analysis outputs (default: sim_read_mapping_analysis/{aligner})'
    )
    parser.add_argument(
        '--output-prefix',
        default=None,
        help='Output prefix for plots and reports (default: metrics_{aligner})'
    )

    args = parser.parse_args()

    # Set default output prefix based on aligner
    if args.output_prefix is None:
        args.output_prefix = f'metrics_{args.aligner}'

    if args.mapper_tsv is None:
        args.mapper_tsv = f'sim_output/{args.aligner}/{args.aligner}_unique_mappings.tsv'
    if args.direct_sam is None:
        args.direct_sam = f'sim_output/{args.aligner}/direct/{args.aligner}_direct.sorted.bam'

    if args.analysis_dir is None:
        analysis_dir = Path('sim_read_mapping_analysis') / args.aligner
    else:
        analysis_dir = Path(args.analysis_dir)
    analysis_dir.mkdir(parents=True, exist_ok=True)

    args.output_prefix = str(analysis_dir / args.output_prefix)

    # Check files exist
    for filepath in [args.fastq_r1, args.mapper_tsv, args.direct_sam]:
        if not Path(filepath).exists():
            print(f"Error: File not found: {filepath}")
            sys.exit(1)

    ground_truth = parse_ground_truth_from_fastq(args.fastq_r1)
    mapper_predictions, partial_reads = parse_mapper_results(args.mapper_tsv)
    direct_predictions = parse_bwa_sam(args.direct_sam)

    # Get all possible genes (labels)
    all_genes = set(ground_truth.values())
    all_genes.add("NONE")
    labels = sorted(list(all_genes))

    # Align all predictions to ground truth
    # Only include reads that have ground truth
    read_ids = [rid for rid in sorted(ground_truth.keys()) if rid not in partial_reads]

    y_true = [ground_truth[rid] for rid in read_ids]
    y_mapper = [mapper_predictions.get(rid, "NONE") for rid in read_ids]
    y_direct = [direct_predictions.get(rid, "NONE") for rid in read_ids]

    # Calculate metrics
    metrics_mapper = calculate_metrics(y_true, y_mapper, labels)
    metrics_direct = calculate_metrics(y_true, y_direct, labels)

    # Create grouped bar chart
    create_grouped_bar_chart(
        metrics_mapper, metrics_direct, labels, args.aligner,
        f"{args.output_prefix}_comparison.png"
    )

    # Save summary to CSV
    summary_data = {
        'Read_ID': read_ids,
        'Ground_Truth': y_true,
        'Mapper_Prediction': y_mapper,
        f'{args.aligner.upper()}_Prediction': y_direct,
        'Mapper_Correct': [1 if y_true[i] == y_mapper[i] else 0 for i in range(len(read_ids))],
        f'{args.aligner.upper()}_Correct': [1 if y_true[i] == y_direct[i] else 0 for i in range(len(read_ids))]
    }
    df = pd.DataFrame(summary_data)
    csv_path = f"{args.output_prefix}_summary.csv"
    df.to_csv(csv_path, index=False)

    # Print summary
    print_summary(metrics_mapper, metrics_direct, args.aligner)
    print(f"\nOutputs saved to: sim_read_mapping_analysis/{args.aligner}/")
    print(f"  Plot: {Path(args.output_prefix).name}_comparison.png")
    print(f"  CSV:  {Path(csv_path).name}")


if __name__ == "__main__":
    main()
