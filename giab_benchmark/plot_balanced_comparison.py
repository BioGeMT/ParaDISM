#!/usr/bin/env python3
"""
Create 2-panel grouped bar plot comparing ParaDISM vs Base Aligner.
Groups: G60 threshold only
Panels: Precision, Recall

SNPs only, GT=1/1. FreeBayes: ploidy=2, min-alternate-count=5.
"""

import csv
import json
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_DIR = SCRIPT_DIR / "vcf_out"

TRUTH_VCF = SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact.vcf.gz"

# Coordinate mapping: gene -> (chr16_start, chr16_end)
GENE_COORDS = {
    "PKD1": (2088707, 2135898),
    "PKD1P1": (16310133, 16334190),
    "PKD1P2": (16356223, 16377507),
    "PKD1P3": (14911550, 14935708),
    "PKD1P4": (18334399, 18352476),
    "PKD1P5": (18374520, 18402014),
    "PKD1P6": (15125138, 15154873),
}

def load_truth_variants():
    """Load truth variants from VCF, converted to gene coordinates."""
    variants = set()
    if not TRUTH_VCF.exists():
        return variants

    cmd = ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(TRUTH_VCF)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
                # Convert to gene coordinates
                for gene, (start, end) in GENE_COORDS.items():
                    if start <= pos <= end:
                        gene_pos = pos - start + 1
                        variants.add((gene, gene_pos, ref, alt))
                        break

    return variants


def load_called_filtered(vcf_path):
    """Load called variants with no additional filtering."""
    variants = set()
    if not vcf_path.exists():
        return variants

    cmd = ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
                variants.add((chrom, pos, ref, alt))

    return variants


def calculate_metrics(truth_variants, called_variants):
    """Calculate precision, recall, and F1."""
    tp = len(truth_variants & called_variants)
    fp = len(called_variants - truth_variants)
    fn = len(truth_variants - called_variants)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0

    return precision, recall, f1, tp, fp, fn


def get_metrics(output_dir, truth, is_baseline=False):
    """Get precision/recall/F1 for a VCF."""
    if is_baseline:
        vcf_path = output_dir / "iteration_1/variant_calling_basealigner_balanced/variants.vcf.gz"
    else:
        vcf_path = output_dir / "final_outputs/variant_calling_balanced/variants.vcf.gz"

    if not vcf_path.exists():
        print(f"Warning: {vcf_path} not found")
        return None, None, 0, 0, 0

    called = load_called_filtered(vcf_path)
    precision, recall, f1, tp, fp, fn = calculate_metrics(truth, called)

    return precision, recall, f1, tp, fp, fn


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print("Loading ground truth variants...")
    truth = load_truth_variants()
    print(f"  {len(truth)} truth variants")

    # Get metrics for G60 only
    g60_dir = SCRIPT_DIR / "giab_hg002_output_bowtie2_G60_min5_qfilters"

    print("\nExtracting metrics...")

    # ParaDISM (final_outputs)
    paradism_g60_prec, paradism_g60_rec, paradism_g60_f1, p_g60_tp, p_g60_fp, _ = get_metrics(g60_dir, truth, is_baseline=False)

    # Base aligner (iteration_1)
    base_g60_prec, base_g60_rec, base_g60_f1, b_g60_tp, b_g60_fp, _ = get_metrics(g60_dir, truth, is_baseline=True)

    # Check for missing data
    missing = []
    if paradism_g60_prec is None: missing.append("ParaDISM G60")
    if base_g60_prec is None: missing.append("Base G60")

    if missing:
        print(f"Error: Missing data for: {', '.join(missing)}")
        print("Run variant calling first: bash giab_benchmark/call_variants_balanced_filters.sh")
        sys.exit(1)

    print(f"G60: ParaDISM TP={p_g60_tp} FP={p_g60_fp} Prec={paradism_g60_prec:.3f} Rec={paradism_g60_rec:.3f} F1={paradism_g60_f1:.3f}")
    print(f"     Base     TP={b_g60_tp} FP={b_g60_fp} Prec={base_g60_prec:.3f} Rec={base_g60_rec:.3f} F1={base_g60_f1:.3f}")

    # Calculate FN for confusion matrices
    p_g60_fn = len(truth) - p_g60_tp
    b_g60_fn = len(truth) - b_g60_tp

    # Save metrics to JSON
    metrics = {
        "filters": {
            "note": "No additional filtering applied in this script."
        },
        "truth_variants": len(truth),
        "G60": {
            "ParaDISM": {
                "TP": p_g60_tp,
                "FP": p_g60_fp,
                "FN": p_g60_fn,
                "precision": round(paradism_g60_prec, 3),
                "recall": round(paradism_g60_rec, 3),
                "f1": round(paradism_g60_f1, 3),
                "confusion_matrix": {
                    "TP": p_g60_tp,
                    "FP": p_g60_fp,
                    "FN": p_g60_fn
                }
            },
            "BaseAligner": {
                "TP": b_g60_tp,
                "FP": b_g60_fp,
                "FN": b_g60_fn,
                "precision": round(base_g60_prec, 3),
                "recall": round(base_g60_rec, 3),
                "f1": round(base_g60_f1, 3),
                "confusion_matrix": {
                    "TP": b_g60_tp,
                    "FP": b_g60_fp,
                    "FN": b_g60_fn
                }
            }
        }
    }

    metrics_path = OUTPUT_DIR / "balanced_comparison_metrics.json"
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"\nMetrics saved to {metrics_path}")

    # Save metrics to CSV
    csv_path = OUTPUT_DIR / "balanced_comparison_metrics.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["dataset", "method", "TP", "FP", "FN", "precision", "recall", "f1"])
        writer.writerow(["G60", "ParaDISM", p_g60_tp, p_g60_fp, p_g60_fn, round(paradism_g60_prec, 3), round(paradism_g60_rec, 3), round(paradism_g60_f1, 3)])
        writer.writerow(["G60", "BaseAligner", b_g60_tp, b_g60_fp, b_g60_fn, round(base_g60_prec, 3), round(base_g60_rec, 3), round(base_g60_f1, 3)])
    print(f"Metrics saved to {csv_path}")

    # Save overall confusion matrices to CSV
    overall_cm_path = OUTPUT_DIR / "confusion_matrices_overall.csv"
    with open(overall_cm_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["dataset", "method", "TP", "FP", "FN"])
        writer.writerow(["G60", "ParaDISM", p_g60_tp, p_g60_fp, p_g60_fn])
        writer.writerow(["G60", "BaseAligner", b_g60_tp, b_g60_fp, b_g60_fn])
    print(f"Confusion matrices saved to {overall_cm_path}")

    # Print confusion matrices
    print("\n" + "="*60)
    print("Confusion Matrices (TP/FP/FN)")
    print("="*60)
    print(f"\nG60 ParaDISM:    TP={p_g60_tp:2d}  FP={p_g60_fp:2d}  FN={p_g60_fn:2d}")
    print(f"G60 BaseAligner: TP={b_g60_tp:2d}  FP={b_g60_fp:2d}  FN={b_g60_fn:2d}")

    # Create confusion matrix plots
    conf_data = [
        ("G60", p_g60_tp, p_g60_fp, p_g60_fn, b_g60_tp, b_g60_fp, b_g60_fn),
    ]

    for name, p_tp, p_fp, p_fn, b_tp, b_fp, b_fn in conf_data:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        for idx, (method, tp, fp, fn) in enumerate([("ParaDISM", p_tp, p_fp, p_fn),
                                                      ("Base Aligner", b_tp, b_fp, b_fn)]):
            ax = axes[idx]

            # Create 2x2 matrix: [[TP, FP], [FN, TN]]
            # For variant calling, TN is not typically computed
            matrix = np.array([[tp, fp], [fn, 0]])
            labels = [["TP", "FP"], ["FN", "—"]]

            im = ax.imshow(matrix, cmap='Blues', vmin=0, vmax=max(tp+fp+fn, 1))

            # Add text annotations
            for i in range(2):
                for j in range(2):
                    if labels[i][j] != "—":
                        text = f"{labels[i][j]}\n{matrix[i, j]}"
                        color = "white" if matrix[i, j] > (tp+fp+fn)/2 else "black"
                    else:
                        text = "—"
                        color = "gray"
                    ax.text(j, i, text, ha="center", va="center", fontsize=14, color=color)

            ax.set_xticks([0, 1])
            ax.set_yticks([0, 1])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_title(f"{method}", fontsize=12, fontweight='bold')

        plt.suptitle(f"Confusion Matrices - {name}\n(SNPs only, GT=1/1; ploidy=2, min-alt=5)",
                     fontsize=12, fontweight='bold', y=1.02)
        plt.tight_layout()

        output_path = OUTPUT_DIR / f"confusion_matrix_{name}"
        plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight')
        plt.savefig(f"{output_path}.pdf", bbox_inches='tight')
        print(f"\nConfusion matrix saved to {output_path}.png and {output_path}.pdf")
        plt.close()

    # Create plot for G60
    datasets = [
        ("G60", paradism_g60_prec, paradism_g60_rec, paradism_g60_f1, base_g60_prec, base_g60_rec, base_g60_f1),
    ]

    for name, p_prec, p_rec, p_f1, b_prec, b_rec, b_f1 in datasets:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))

        x = np.arange(1)
        width = 0.35

        # Panel 1: Precision
        ax1 = axes[0]
        bars1 = ax1.bar(x - width/2, [p_prec], width, label='ParaDISM', color='#2E86AB')
        bars2 = ax1.bar(x + width/2, [b_prec], width, label='Base Aligner', color='#A23B72')

        ax1.set_ylabel('Precision', fontsize=12)
        ax1.set_title('Precision', fontsize=14, fontweight='bold')
        ax1.set_xticks([])
        ax1.set_ylim(0, 1)
        ax1.legend(loc='upper right')

        for bar in bars1:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)
        for bar in bars2:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)

        # Panel 2: Recall
        ax2 = axes[1]
        bars3 = ax2.bar(x - width/2, [p_rec], width, label='ParaDISM', color='#2E86AB')
        bars4 = ax2.bar(x + width/2, [b_rec], width, label='Base Aligner', color='#A23B72')

        ax2.set_ylabel('Recall', fontsize=12)
        ax2.set_title('Recall', fontsize=14, fontweight='bold')
        ax2.set_xticks([])
        ax2.set_ylim(0, 1)
        ax2.legend(loc='upper right')

        for bar in bars3:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)
        for bar in bars4:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)

        # Panel 3: F1
        ax3 = axes[2]
        bars5 = ax3.bar(x - width/2, [p_f1], width, label='ParaDISM', color='#2E86AB')
        bars6 = ax3.bar(x + width/2, [b_f1], width, label='Base Aligner', color='#A23B72')

        ax3.set_ylabel('F1', fontsize=12)
        ax3.set_title('F1', fontsize=14, fontweight='bold')
        ax3.set_xticks([])
        ax3.set_ylim(0, 1)
        ax3.legend(loc='upper right')

        for bar in bars5:
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)
        for bar in bars6:
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)

        plt.suptitle(f'ParaDISM vs Base Aligner - {name}\n(SNPs only, GT=1/1; ploidy=2, min-alt=5)',
                     fontsize=12, fontweight='bold', y=1.05)
        plt.tight_layout()

        output_path = OUTPUT_DIR / f"balanced_comparison_{name}"
        plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight')
        plt.savefig(f"{output_path}.pdf", bbox_inches='tight')
        print(f"\nPlot saved to {output_path}.png and {output_path}.pdf")
        plt.close()


def per_gene_analysis():
    """Generate per-gene metrics."""
    print("\n" + "="*60)
    print("Per-Gene Analysis")
    print("="*60)

    # Create output directory
    per_gene_dir = OUTPUT_DIR / "per_gene_metrics"
    per_gene_dir.mkdir(exist_ok=True)

    # Load truth per gene
    truth_by_gene = {gene: set() for gene in GENE_COORDS.keys()}
    cmd = ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(TRUTH_VCF)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 4:
                chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
                for gene, (start, end) in GENE_COORDS.items():
                    if start <= pos <= end:
                        gene_pos = pos - start + 1
                        truth_by_gene[gene].add((gene, gene_pos, ref, alt))
                        break

    # Process G60 only
    all_rows = []

    for dataset in ["G60"]:
        dir_name = SCRIPT_DIR / f"giab_hg002_output_bowtie2_{dataset}_min5_qfilters"

        for method, vcf_subpath in [("ParaDISM", "final_outputs/variant_calling_balanced/variants.vcf.gz"),
                                     ("BaseAligner", "iteration_1/variant_calling_basealigner_balanced/variants.vcf.gz")]:
            vcf_path = dir_name / vcf_subpath
            if not vcf_path.exists():
                continue

            # Load called variants by gene
            called_by_gene = {gene: set() for gene in GENE_COORDS.keys()}
            cmd = ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(vcf_path)]
            result = subprocess.run(cmd, capture_output=True, text=True)

            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
                        if chrom in called_by_gene:
                            called_by_gene[chrom].add((chrom, pos, ref, alt))

            # Calculate per-gene metrics
            for gene in GENE_COORDS.keys():
                truth = truth_by_gene[gene]
                called = called_by_gene[gene]

                tp = len(truth & called)
                fp = len(called - truth)
                fn = len(truth - called)

                precision = tp / (tp + fp) if (tp + fp) > 0 else 0
                recall = tp / (tp + fn) if (tp + fn) > 0 else 0
                f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0

                all_rows.append({
                    "dataset": dataset,
                    "method": method,
                    "gene": gene,
                    "truth_count": len(truth),
                    "TP": tp,
                    "FP": fp,
                    "FN": fn,
                    "precision": round(precision, 3),
                    "recall": round(recall, 3),
                    "f1": round(f1, 3)
                })

    # Save per-gene CSV
    csv_path = per_gene_dir / "per_gene_metrics.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["dataset", "method", "gene", "truth_count", "TP", "FP", "FN", "precision", "recall", "f1"],
        )
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nPer-gene metrics saved to {csv_path}")

    # Save per-gene JSON
    json_path = per_gene_dir / "per_gene_metrics.json"
    with open(json_path, 'w') as f:
        json.dump(all_rows, f, indent=2)
    print(f"Per-gene metrics saved to {json_path}")

    # Print summary
    print(f"\nPer-gene summary:")
    for gene in GENE_COORDS.keys():
        truth_count = len(truth_by_gene[gene])
        if truth_count > 0:
            print(f"  {gene}: {truth_count} truth variants")

    # Generate per-gene plots for G60 only
    print("\nGenerating per-gene plots (G60)...")

    # Filter to G60 and genes with truth variants
    g60_data = [r for r in all_rows if r["dataset"] == "G60" and r["truth_count"] > 0]
    genes_with_truth = sorted(set(r["gene"] for r in g60_data))

    for gene in genes_with_truth:
        gene_rows = [r for r in g60_data if r["gene"] == gene]
        paradism_row = next((r for r in gene_rows if r["method"] == "ParaDISM"), None)
        base_row = next((r for r in gene_rows if r["method"] == "BaseAligner"), None)

        if not paradism_row or not base_row:
            continue

        fig, axes = plt.subplots(1, 2, figsize=(8, 4))

        x = np.arange(1)
        width = 0.35

        # Panel 1: Precision
        ax1 = axes[0]
        bars1 = ax1.bar(x - width/2, [paradism_row["precision"]], width, label='ParaDISM', color='#2E86AB')
        bars2 = ax1.bar(x + width/2, [base_row["precision"]], width, label='Base Aligner', color='#A23B72')

        ax1.set_ylabel('Precision', fontsize=12)
        ax1.set_title('Precision', fontsize=14, fontweight='bold')
        ax1.set_xticks([])
        ax1.set_ylim(0, 1)
        ax1.legend(loc='upper right')

        for bar in bars1:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)
        for bar in bars2:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)

        # Panel 2: Recall
        ax2 = axes[1]
        bars3 = ax2.bar(x - width/2, [paradism_row["recall"]], width, label='ParaDISM', color='#2E86AB')
        bars4 = ax2.bar(x + width/2, [base_row["recall"]], width, label='Base Aligner', color='#A23B72')

        ax2.set_ylabel('Recall', fontsize=12)
        ax2.set_title('Recall', fontsize=14, fontweight='bold')
        ax2.set_xticks([])
        ax2.set_ylim(0, 1)
        ax2.legend(loc='upper right')

        for bar in bars3:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)
        for bar in bars4:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                     f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)

        plt.suptitle(f'{gene} - G60\n(n={paradism_row["truth_count"]} truth variants)',
                     fontsize=12, fontweight='bold', y=1.05)
        plt.tight_layout()

        output_path = per_gene_dir / f"G60_{gene}"
        plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight')
        plt.savefig(f"{output_path}.pdf", bbox_inches='tight')
        print(f"  Saved {output_path}.png")
        plt.close()

        # Per-gene confusion matrix plot
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        for idx, (method, row) in enumerate([("ParaDISM", paradism_row), ("Base Aligner", base_row)]):
            ax = axes[idx]
            tp = row["TP"]
            fp = row["FP"]
            fn = row["FN"]

            matrix = np.array([[tp, fp], [fn, 0]])
            labels = [["TP", "FP"], ["FN", "—"]]

            ax.imshow(matrix, cmap='Blues', vmin=0, vmax=max(tp + fp + fn, 1))

            for i in range(2):
                for j in range(2):
                    if labels[i][j] != "—":
                        text = f"{labels[i][j]}\n{matrix[i, j]}"
                        color = "white" if matrix[i, j] > (tp + fp + fn) / 2 else "black"
                    else:
                        text = "—"
                        color = "gray"
                    ax.text(j, i, text, ha="center", va="center", fontsize=14, color=color)

            ax.set_xticks([0, 1])
            ax.set_yticks([0, 1])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_title(f"{method}", fontsize=12, fontweight='bold')

        plt.suptitle(f"Confusion Matrices - {gene} (G60)",
                     fontsize=12, fontweight='bold', y=1.02)
        plt.tight_layout()

        output_cm = per_gene_dir / f"G60_{gene}_confusion"
        plt.savefig(f"{output_cm}.png", dpi=150, bbox_inches='tight')
        plt.savefig(f"{output_cm}.pdf", bbox_inches='tight')
        print(f"  Saved {output_cm}.png")
        plt.close()


if __name__ == "__main__":
    main()
    per_gene_analysis()
