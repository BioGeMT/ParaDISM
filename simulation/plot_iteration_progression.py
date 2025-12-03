#!/usr/bin/env python3
"""
Plot precision, recall, and specificity across iterations.
Shows ParaDISM improvement vs constant BWA-MEM2 baseline.
"""

import argparse
import sys
from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import confusion_matrix


def parse_metrics_from_csv(csv_path: Path) -> dict:
    """Parse overall metrics from read_mapping_analysis CSV output."""
    if not csv_path.exists():
        return None
    
    df = pd.read_csv(csv_path)
    
    # Find overall metrics row
    overall_row = df[df['Ground_Truth'] == 'Overall']
    if overall_row.empty:
        # Try to calculate from all rows
        # Get mapper and direct predictions
        mapper_pred = df['Mapper_Prediction']
        direct_pred = df.iloc[:, 3]  # Direct prediction column
        
        # Calculate confusion matrix metrics
        from sklearn.metrics import confusion_matrix, precision_recall_fscore_support
        
        # Get all unique labels
        all_labels = sorted(set(df['Ground_Truth']) | set(mapper_pred) | set(direct_pred))
        
        # Calculate metrics for mapper
        cm_mapper = confusion_matrix(df['Ground_Truth'], mapper_pred, labels=all_labels)
        prec_mapper, rec_mapper, _, _ = precision_recall_fscore_support(
            df['Ground_Truth'], mapper_pred, labels=all_labels, average='weighted', zero_division=0
        )
        
        # Calculate specificity for mapper
        tn_mapper = np.sum(cm_mapper) - np.sum(cm_mapper.diagonal()) - np.sum(np.sum(cm_mapper, axis=1) - cm_mapper.diagonal())
        fp_mapper = np.sum(np.sum(cm_mapper, axis=0) - cm_mapper.diagonal())
        spec_mapper = tn_mapper / (tn_mapper + fp_mapper) if (tn_mapper + fp_mapper) > 0 else 0.0
        
        # Calculate metrics for direct
        cm_direct = confusion_matrix(df['Ground_Truth'], direct_pred, labels=all_labels)
        prec_direct, rec_direct, _, _ = precision_recall_fscore_support(
            df['Ground_Truth'], direct_pred, labels=all_labels, average='weighted', zero_division=0
        )
        
        # Calculate specificity for direct
        tn_direct = np.sum(cm_direct) - np.sum(cm_direct.diagonal()) - np.sum(np.sum(cm_direct, axis=1) - cm_direct.diagonal())
        fp_direct = np.sum(np.sum(cm_direct, axis=0) - cm_direct.diagonal())
        spec_direct = tn_direct / (tn_direct + fp_direct) if (tn_direct + fp_direct) > 0 else 0.0
        
        return {
            'mapper': {'precision': prec_mapper, 'recall': rec_mapper, 'specificity': spec_mapper},
            'direct': {'precision': prec_direct, 'recall': rec_direct, 'specificity': spec_direct}
        }
    
    # If we have overall row, extract metrics
    # This is a simplified approach - actual CSV structure may vary
    return {
        'mapper': {'precision': 0.0, 'recall': 0.0, 'specificity': 0.0},
        'direct': {'precision': 0.0, 'recall': 0.0, 'specificity': 0.0}
    }


def calculate_specificity(y_true, y_pred, labels):
    """Calculate overall specificity."""
    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    
    # Overall specificity: sum of all TNs / (sum of all TNs + sum of all FPs)
    total_tn = 0
    total_fp = 0
    
    for i, label in enumerate(labels):
        # TN for this class: all samples not in this class that were correctly not predicted as this class
        tn = np.sum(cm) - np.sum(cm[i, :]) - np.sum(cm[:, i]) + cm[i, i]
        # FP for this class: all samples not in this class that were incorrectly predicted as this class
        fp = np.sum(cm[:, i]) - cm[i, i]
        total_tn += tn
        total_fp += fp
    
    return total_tn / (total_tn + total_fp) if (total_tn + total_fp) > 0 else 0.0


def extract_metrics_from_summary(summary_path: Path) -> dict:
    """Extract overall metrics from summary CSV."""
    if not summary_path.exists():
        return None
    
    # Read the summary CSV from read_mapping_analysis
    # Format: Read_ID, Ground_Truth, Mapper_Prediction, BWA-MEM2_Prediction, ...
    df = pd.read_csv(summary_path)
    
    from sklearn.metrics import precision_recall_fscore_support
    
    # Get predictions
    ground_truth = df['Ground_Truth']

    # Handle both 'Mapper_Prediction' and 'ParaDISM_Prediction' column names
    if 'Mapper_Prediction' in df.columns:
        mapper_pred = df['Mapper_Prediction']
    elif 'ParaDISM_Prediction' in df.columns:
        mapper_pred = df['ParaDISM_Prediction']
    else:
        raise KeyError("Could not find Mapper_Prediction or ParaDISM_Prediction column")

    # Find direct prediction column
    direct_pred_col = None
    for col in df.columns:
        if col.endswith('_Prediction') and col not in ['Mapper_Prediction', 'ParaDISM_Prediction']:
            direct_pred_col = col
            break
    
    if direct_pred_col is None:
        # Try BWA-MEM2_Prediction or similar
        if 'BWA-MEM2_Prediction' in df.columns:
            direct_pred_col = 'BWA-MEM2_Prediction'
        else:
            direct_pred = df.iloc[:, 3]  # Fallback to 4th column
    else:
        direct_pred = df[direct_pred_col]
    
    # Get all labels
    all_labels = sorted(set(ground_truth) | set(mapper_pred) | set(direct_pred))
    
    # Calculate metrics for mapper (ParaDISM)
    prec_mapper, rec_mapper, _, support_mapper = precision_recall_fscore_support(
        ground_truth, mapper_pred, labels=all_labels, average='weighted', zero_division=0
    )
    spec_mapper = calculate_specificity(ground_truth, mapper_pred, all_labels)
    
    # Calculate metrics for direct (BWA-MEM2)
    prec_direct, rec_direct, _, support_direct = precision_recall_fscore_support(
        ground_truth, direct_pred, labels=all_labels, average='weighted', zero_division=0
    )
    spec_direct = calculate_specificity(ground_truth, direct_pred, all_labels)
    
    return {
        'mapper': {
            'precision': prec_mapper,
            'recall': rec_mapper,
            'specificity': spec_mapper
        },
        'direct': {
            'precision': prec_direct,
            'recall': rec_direct,
            'specificity': spec_direct
        }
    }


def find_max_iteration(analysis_dir: Path, seed: int, aligner: str) -> int:
    """Find the highest iteration with a summary CSV for this seed/aligner."""
    pattern = re.compile(rf"seed_{seed}_{aligner}_iter(\d+)_summary\.csv$")
    max_iter = 0
    if analysis_dir.exists():
        for path in analysis_dir.glob(f"seed_{seed}_{aligner}_iter*_summary.csv"):
            match = pattern.search(path.name)
            if match:
                try:
                    num = int(match.group(1))
                    if num > max_iter:
                        max_iter = num
                except ValueError:
                    continue
    return max_iter


def collect_iteration_metrics_for_seed(analysis_dir: Path, seed: int, aligner: str) -> tuple[list, list, list, list, list, list]:
    """Collect metrics for each iteration found on disk (no padding beyond max)."""
    paradism_precision = []
    paradism_recall = []
    paradism_specificity = []
    bwa_precision = []
    bwa_recall = []
    bwa_specificity = []
    
    max_iter = find_max_iteration(analysis_dir, seed, aligner)
    if max_iter == 0:
        print(f"Warning: No iteration summaries found for seed {seed} ({aligner}) in {analysis_dir}", file=sys.stderr)
        return (paradism_precision, paradism_recall, paradism_specificity,
                bwa_precision, bwa_recall, bwa_specificity)
    
    bwa_prec_base = None
    bwa_rec_base = None
    bwa_spec_base = None
    
    for iter_num in range(1, max_iter + 1):
        iter_csv = analysis_dir / f"seed_{seed}_{aligner}_iter{iter_num}_summary.csv"
        if iter_csv.exists():
            metrics = extract_metrics_from_summary(iter_csv)
            if metrics:
                paradism_precision.append(metrics['mapper']['precision'])
                paradism_recall.append(metrics['mapper']['recall'])
                paradism_specificity.append(metrics['mapper']['specificity'])
                
                if iter_num == 1:
                    bwa_prec_base = metrics['direct']['precision']
                    bwa_rec_base = metrics['direct']['recall']
                    bwa_spec_base = metrics['direct']['specificity']
                
                # BWA baseline stays constant across iterations for this seed
                bwa_precision.append(bwa_prec_base)
                bwa_recall.append(bwa_rec_base)
                bwa_specificity.append(bwa_spec_base)
            else:
                print(f"Warning: Could not parse metrics for {iter_csv}", file=sys.stderr)
        else:
            # Missing iteration file; stop at the last available
            break
    
    return (paradism_precision, paradism_recall, paradism_specificity,
            bwa_precision, bwa_recall, bwa_specificity)


def collect_metrics_across_seeds(
    sim_output_base: Path,
    seed_start: int,
    seed_end: int,
    aligner: str
) -> tuple[list, list, list, list, list, list]:
    """Collect metrics across all seeds and calculate mean/std for error bars."""
    import statistics

    # Determine maximum iteration observed across seeds
    max_iter_overall = 0
    per_seed_metrics = {}
    for seed in range(seed_start, seed_end + 1):
        analysis_dir = sim_output_base / f"seed_{seed}" / aligner / "read_mapping_analysis"
        (para_prec, para_rec, para_spec,
         bwa_prec, bwa_rec, bwa_spec) = collect_iteration_metrics_for_seed(
            analysis_dir, seed, aligner
        )

        per_seed_metrics[seed] = (para_prec, para_rec, para_spec, bwa_prec, bwa_rec, bwa_spec)
        max_iter_overall = max(max_iter_overall, len(para_prec))

    if max_iter_overall == 0:
        return ([], [], [], [], [], [], [], [], [], [], [], [], 0)

    # Collect metrics for each iteration index only from seeds that have that iteration
    all_para_prec = [[] for _ in range(max_iter_overall)]
    all_para_rec = [[] for _ in range(max_iter_overall)]
    all_para_spec = [[] for _ in range(max_iter_overall)]
    all_bwa_prec = [[] for _ in range(max_iter_overall)]
    all_bwa_rec = [[] for _ in range(max_iter_overall)]
    all_bwa_spec = [[] for _ in range(max_iter_overall)]

    for seed, (para_prec, para_rec, para_spec, bwa_prec, bwa_rec, bwa_spec) in per_seed_metrics.items():
        for i in range(len(para_prec)):
            all_para_prec[i].append(para_prec[i])
            all_para_rec[i].append(para_rec[i])
            all_para_spec[i].append(para_spec[i])
            all_bwa_prec[i].append(bwa_prec[i])
            all_bwa_rec[i].append(bwa_rec[i])
            all_bwa_spec[i].append(bwa_spec[i])

    # Calculate mean and std for each iteration (only over seeds that reached that iteration)
    para_prec_mean = []
    para_prec_std = []
    para_rec_mean = []
    para_rec_std = []
    para_spec_mean = []
    para_spec_std = []
    bwa_prec_mean = []
    bwa_prec_std = []
    bwa_rec_mean = []
    bwa_rec_std = []
    bwa_spec_mean = []
    bwa_spec_std = []
    
    for i in range(max_iter_overall):
        # ParaDISM metrics
        if all_para_prec[i]:
            para_prec_mean.append(statistics.mean(all_para_prec[i]))
            para_prec_std.append(statistics.stdev(all_para_prec[i]) if len(all_para_prec[i]) > 1 else 0.0)
            para_rec_mean.append(statistics.mean(all_para_rec[i]))
            para_rec_std.append(statistics.stdev(all_para_rec[i]) if len(all_para_rec[i]) > 1 else 0.0)
            para_spec_mean.append(statistics.mean(all_para_spec[i]))
            para_spec_std.append(statistics.stdev(all_para_spec[i]) if len(all_para_spec[i]) > 1 else 0.0)
        else:
            para_prec_mean.append(np.nan)
            para_prec_std.append(0.0)
            para_rec_mean.append(np.nan)
            para_rec_std.append(0.0)
            para_spec_mean.append(np.nan)
            para_spec_std.append(0.0)
        
        # BWA metrics
        if all_bwa_prec[i]:
            bwa_prec_mean.append(statistics.mean(all_bwa_prec[i]))
            bwa_prec_std.append(statistics.stdev(all_bwa_prec[i]) if len(all_bwa_prec[i]) > 1 else 0.0)
            bwa_rec_mean.append(statistics.mean(all_bwa_rec[i]))
            bwa_rec_std.append(statistics.stdev(all_bwa_rec[i]) if len(all_bwa_rec[i]) > 1 else 0.0)
            bwa_spec_mean.append(statistics.mean(all_bwa_spec[i]))
            bwa_spec_std.append(statistics.stdev(all_bwa_spec[i]) if len(all_bwa_spec[i]) > 1 else 0.0)
        else:
            bwa_prec_mean.append(np.nan)
            bwa_prec_std.append(0.0)
            bwa_rec_mean.append(np.nan)
            bwa_rec_std.append(0.0)
            bwa_spec_mean.append(np.nan)
            bwa_spec_std.append(0.0)
    
    return (
        para_prec_mean, para_prec_std,
        para_rec_mean, para_rec_std,
        para_spec_mean, para_spec_std,
        bwa_prec_mean, bwa_prec_std,
        bwa_rec_mean, bwa_rec_std,
        bwa_spec_mean, bwa_spec_std,
        max_iter_overall
    )


def create_iteration_plots(
    sim_output_base: Path,
    seed_start: int,
    seed_end: int,
    output_dir: Path,
    aligner: str,
    iterations: int | None = None,
) -> None:
    """Create plots showing metric progression across iterations with error bars."""
    
    # Collect metrics across all seeds
    (
        para_prec_mean, para_prec_std,
        para_rec_mean, para_rec_std,
        para_spec_mean, para_spec_std,
        bwa_prec_mean, bwa_prec_std,
        bwa_rec_mean, bwa_rec_std,
        bwa_spec_mean, bwa_spec_std,
        max_iter_overall
    ) = collect_metrics_across_seeds(
        sim_output_base, seed_start, seed_end, aligner
    )
    
    if not para_prec_mean:
        print(f"Error: No metrics found", file=sys.stderr)
        return

    # Determine how many iterations to plot: auto-detect unless explicitly capped
    target_iterations = iterations if iterations and iterations > 0 else max_iter_overall
    if target_iterations > max_iter_overall:
        print(f"Warning: requested iterations={target_iterations} exceeds observed max={max_iter_overall}; capping to {max_iter_overall}", file=sys.stderr)
        target_iterations = max_iter_overall

    # Truncate metrics to the chosen iteration count
    def trim(arr: list[float]) -> list[float]:
        return arr[:target_iterations]

    para_prec_mean = trim(para_prec_mean)
    para_prec_std = trim(para_prec_std)
    para_rec_mean = trim(para_rec_mean)
    para_rec_std = trim(para_rec_std)
    para_spec_mean = trim(para_spec_mean)
    para_spec_std = trim(para_spec_std)
    bwa_prec_mean = trim(bwa_prec_mean)
    bwa_prec_std = trim(bwa_prec_std)
    bwa_rec_mean = trim(bwa_rec_mean)
    bwa_rec_std = trim(bwa_rec_std)
    bwa_spec_mean = trim(bwa_spec_mean)
    bwa_spec_std = trim(bwa_spec_std)
    
    # Create iteration numbers (1 to target_iterations, 1-based indexing)
    iter_nums = list(range(1, target_iterations + 1))
    
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
    
    # Plot 1: Precision
    axes[0].errorbar(iter_nums, para_prec_mean, yerr=para_prec_std, 
                     fmt='o-', color='red', linewidth=2, markersize=8, 
                     capsize=5, capthick=2, label='ParaDISM', alpha=0.8)
    # BWA is constant across iterations - error bar only at iteration 1, then straight line
    axes[0].errorbar([iter_nums[0]], [bwa_prec_mean[0]], yerr=[bwa_prec_std[0]],
                     fmt='s', color='#4682B4', linewidth=2, markersize=8,
                     capsize=5, capthick=2, label='BWA-MEM2 (direct)', alpha=0.8)
    axes[0].plot(iter_nums, [bwa_prec_mean[0]] * len(iter_nums), ':', 
                 color='#4682B4', linewidth=2, alpha=0.8)
    axes[0].set_xlabel('Iteration', weight='bold')
    axes[0].set_ylabel('Precision', weight='bold')
    axes[0].set_title(f'Precision Across Iterations (n={seed_end - seed_start + 1} seeds)', weight='bold')
    axes[0].set_xlim([0.5, iterations + 0.5])
    axes[0].set_ylim([0, 1.05])
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    axes[0].set_xticks(iter_nums)
    
    # Plot 2: Recall
    axes[1].errorbar(iter_nums, para_rec_mean, yerr=para_rec_std,
                     fmt='o-', color='red', linewidth=2, markersize=8,
                     capsize=5, capthick=2, label='ParaDISM', alpha=0.8)
    # BWA is constant across iterations - error bar only at iteration 1, then straight line
    axes[1].errorbar([iter_nums[0]], [bwa_rec_mean[0]], yerr=[bwa_rec_std[0]],
                     fmt='s', color='#4682B4', linewidth=2, markersize=8,
                     capsize=5, capthick=2, label='BWA-MEM2 (direct)', alpha=0.8)
    axes[1].plot(iter_nums, [bwa_rec_mean[0]] * len(iter_nums), ':', 
                 color='#4682B4', linewidth=2, alpha=0.8)
    axes[1].set_xlabel('Iteration', weight='bold')
    axes[1].set_ylabel('Recall', weight='bold')
    axes[1].set_title(f'Recall Across Iterations (n={seed_end - seed_start + 1} seeds)', weight='bold')
    axes[1].set_xlim([0.5, iterations + 0.5])
    axes[1].set_ylim([0, 1.05])
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    axes[1].set_xticks(iter_nums)
    
    # Plot 3: Specificity
    axes[2].errorbar(iter_nums, para_spec_mean, yerr=para_spec_std,
                     fmt='o-', color='red', linewidth=2, markersize=8,
                     capsize=5, capthick=2, label='ParaDISM', alpha=0.8)
    # BWA is constant across iterations - error bar only at iteration 1, then straight line
    axes[2].errorbar([iter_nums[0]], [bwa_spec_mean[0]], yerr=[bwa_spec_std[0]],
                     fmt='s', color='#4682B4', linewidth=2, markersize=8,
                     capsize=5, capthick=2, label='BWA-MEM2 (direct)', alpha=0.8)
    axes[2].plot(iter_nums, [bwa_spec_mean[0]] * len(iter_nums), ':', 
                 color='#4682B4', linewidth=2, alpha=0.8)
    axes[2].set_xlabel('Iteration', weight='bold')
    axes[2].set_ylabel('Specificity', weight='bold')
    axes[2].set_title(f'Specificity Across Iterations (n={seed_end - seed_start + 1} seeds)', weight='bold')
    axes[2].set_xlim([0.5, iterations + 0.5])
    axes[2].set_ylim([0, 1.05])
    axes[2].grid(alpha=0.3)
    axes[2].legend()
    axes[2].set_xticks(iter_nums)
    
    plt.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f'iteration_progression_{aligner}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save CSV with iteration progression data
    csv_path = output_dir / f'iteration_progression_{aligner}.csv'
    import pandas as pd
    progression_data = {
        'Iteration': iter_nums,
        'ParaDISM_Precision_Mean': para_prec_mean,
        'ParaDISM_Precision_Std': para_prec_std,
        'ParaDISM_Recall_Mean': para_rec_mean,
        'ParaDISM_Recall_Std': para_rec_std,
        'ParaDISM_Specificity_Mean': para_spec_mean,
        'ParaDISM_Specificity_Std': para_spec_std,
        'BWA-MEM2_Precision_Mean': bwa_prec_mean,
        'BWA-MEM2_Precision_Std': bwa_prec_std,
        'BWA-MEM2_Recall_Mean': bwa_rec_mean,
        'BWA-MEM2_Recall_Std': bwa_rec_std,
        'BWA-MEM2_Specificity_Mean': bwa_spec_mean,
        'BWA-MEM2_Specificity_Std': bwa_spec_std,
    }
    df = pd.DataFrame(progression_data)
    df.to_csv(csv_path, index=False)
    
    print(f"Saved iteration progression plot: {output_path}")
    print(f"Saved iteration progression CSV: {csv_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot precision, recall, and specificity across iterations with error bars'
    )
    parser.add_argument(
        '--sim-output-base',
        type=Path,
        required=True,
        help='Base directory containing seed_* subdirectories'
    )
    parser.add_argument(
        '--seed-start',
        type=int,
        required=True,
        help='Starting seed number'
    )
    parser.add_argument(
        '--seed-end',
        type=int,
        required=True,
        help='Ending seed number'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        required=True,
        help='Output directory for plots'
    )
    parser.add_argument(
        '--aligner',
        default='bwa-mem2',
        help='Aligner name'
    )
    parser.add_argument(
        '--iterations',
        type=int,
        default=0,
        help='Number of iterations to plot (0 = auto-detect maximum per seed and cap to overall max)'
    )
    
    args = parser.parse_args()
    
    create_iteration_plots(
        args.sim_output_base,
        args.seed_start,
        args.seed_end,
        args.output_dir,
        args.aligner,
        args.iterations
    )


if __name__ == '__main__':
    main()
