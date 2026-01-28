#!/usr/bin/env python3
"""
Create 2-panel grouped bar plot comparing ParaDISM vs Base Aligner.
Groups: G40, G60 thresholds
Panels: Precision, Recall
"""

import matplotlib.pyplot as plt
import numpy as np
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent

TRUTH_VCF = SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_illumina_only.vcf.gz"


def load_variants(vcf_path):
    """Load variants from VCF as set of (chrom, pos, ref, alt)."""
    variants = set()
    if not vcf_path.exists():
        return variants

    cmd = ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.strip().split('\n'):
        if line:
            parts = line.split('\t')
            if len(parts) >= 4:
                variants.add((parts[0], parts[1], parts[2], parts[3]))

    return variants


def calculate_metrics(truth_variants, called_variants):
    """Calculate precision and recall."""
    tp = len(truth_variants & called_variants)
    fp = len(called_variants - truth_variants)
    fn = len(truth_variants - called_variants)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0

    return precision, recall, tp, fp, fn


def get_metrics(output_dir, is_baseline=False):
    """Get precision/recall for a VCF."""
    if is_baseline:
        vcf_path = output_dir / "iteration_1/variant_calling_balanced/variants.vcf.gz"
    else:
        vcf_path = output_dir / "final_outputs/variant_calling_balanced/variants.vcf.gz"

    if not vcf_path.exists():
        print(f"Warning: {vcf_path} not found")
        return None, None

    truth = load_variants(TRUTH_VCF)
    called = load_variants(vcf_path)

    precision, recall, tp, fp, fn = calculate_metrics(truth, called)

    return precision, recall


def main():
    print("Loading ground truth...")
    truth = load_variants(TRUTH_VCF)
    print(f"  {len(truth)} truth variants")

    # Get metrics for G40 and G60
    g40_dir = SCRIPT_DIR / "giab_hg002_output_bowtie2_G40_min5_qfilters"
    g60_dir = SCRIPT_DIR / "giab_hg002_output_bowtie2_G60_min5_qfilters"

    print("\nExtracting metrics...")

    # ParaDISM (final_outputs)
    paradism_g40_prec, paradism_g40_rec = get_metrics(g40_dir, is_baseline=False)
    paradism_g60_prec, paradism_g60_rec = get_metrics(g60_dir, is_baseline=False)

    # Base aligner (iteration_1)
    base_g40_prec, base_g40_rec = get_metrics(g40_dir, is_baseline=True)
    base_g60_prec, base_g60_rec = get_metrics(g60_dir, is_baseline=True)

    # Check for missing data
    missing = []
    if paradism_g40_prec is None: missing.append("ParaDISM G40")
    if paradism_g60_prec is None: missing.append("ParaDISM G60")
    if base_g40_prec is None: missing.append("Base G40")
    if base_g60_prec is None: missing.append("Base G60")

    if missing:
        print(f"Error: Missing data for: {', '.join(missing)}")
        print("Run variant calling first: bash giab_benchmark/call_variants_balanced_filters.sh")
        sys.exit(1)

    print(f"\nG40: ParaDISM Prec={paradism_g40_prec:.3f} Rec={paradism_g40_rec:.3f}, Base Prec={base_g40_prec:.3f} Rec={base_g40_rec:.3f}")
    print(f"G60: ParaDISM Prec={paradism_g60_prec:.3f} Rec={paradism_g60_rec:.3f}, Base Prec={base_g60_prec:.3f} Rec={base_g60_rec:.3f}")

    # Create grouped bar plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    x = np.arange(2)  # G40, G60
    width = 0.35

    # Panel 1: Precision
    ax1 = axes[0]
    paradism_prec = [paradism_g40_prec, paradism_g60_prec]
    base_prec = [base_g40_prec, base_g60_prec]

    bars1 = ax1.bar(x - width/2, paradism_prec, width, label='ParaDISM', color='#2E86AB')
    bars2 = ax1.bar(x + width/2, base_prec, width, label='Base Aligner', color='#A23B72')

    ax1.set_ylabel('Precision', fontsize=12)
    ax1.set_title('Precision', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(['G40', 'G60'], fontsize=11)
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
    paradism_rec = [paradism_g40_rec, paradism_g60_rec]
    base_rec = [base_g40_rec, base_g60_rec]

    bars3 = ax2.bar(x - width/2, paradism_rec, width, label='ParaDISM', color='#2E86AB')
    bars4 = ax2.bar(x + width/2, base_rec, width, label='Base Aligner', color='#A23B72')

    ax2.set_ylabel('Recall', fontsize=12)
    ax2.set_title('Recall', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(['G40', 'G60'], fontsize=11)
    ax2.set_ylim(0, 1)
    ax2.legend(loc='upper right')

    for bar in bars3:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)
    for bar in bars4:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=10)

    plt.suptitle('ParaDISM vs Base Aligner (Balanced Filters)', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    output_path = SCRIPT_DIR / "balanced_comparison_plot"
    plt.savefig(f"{output_path}.png", dpi=150, bbox_inches='tight')
    plt.savefig(f"{output_path}.pdf", bbox_inches='tight')
    print(f"\nPlot saved to {output_path}.png and {output_path}.pdf")


if __name__ == "__main__":
    main()
