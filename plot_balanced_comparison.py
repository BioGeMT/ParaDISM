#!/usr/bin/env python3
"""
Create 2-panel bar plot comparing ParaDISM vs Base Aligner for balanced filters.
Panels: Precision, Recall
"""

import matplotlib.pyplot as plt
import numpy as np

# Data from precision_recall_illumina_only.tsv (balanced filters only)
# Truth set: 42 variants from HG002_PKD1_genes_SNPs_illumina_only.vcf.gz
truth_count = 42

# ParaDISM balanced
paradism_precision = 0.6957
paradism_recall = 0.3810

# Base aligner balanced
basealigner_precision = 0.6842
basealigner_recall = 0.3095

print(f"ParaDISM balanced: Precision={paradism_precision:.4f}, Recall={paradism_recall:.4f}")
print(f"Base aligner balanced: Precision={basealigner_precision:.4f}, Recall={basealigner_recall:.4f}")

# Create figure with 2 panels
fig, axes = plt.subplots(1, 2, figsize=(8, 5))

methods = ['ParaDISM', 'Base Aligner']
colors = ['#2E86AB', '#A23B72']  # Blue for ParaDISM, magenta for base aligner
x = np.arange(len(methods))
width = 0.6

# Panel 1: Precision
ax1 = axes[0]
precision_vals = [paradism_precision, basealigner_precision]
bars1 = ax1.bar(x, precision_vals, width, color=colors)
ax1.set_ylabel('Precision', fontsize=12)
ax1.set_title('Precision', fontsize=14, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(methods, fontsize=11)
ax1.set_ylim(0, 1)
for bar, val in zip(bars1, precision_vals):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
             f'{val:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

# Panel 2: Recall
ax2 = axes[1]
recall_vals = [paradism_recall, basealigner_recall]
bars2 = ax2.bar(x, recall_vals, width, color=colors)
ax2.set_ylabel('Recall', fontsize=12)
ax2.set_title('Recall', fontsize=14, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(methods, fontsize=11)
ax2.set_ylim(0, 1)
for bar, val in zip(bars2, recall_vals):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
             f'{val:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

plt.suptitle('ParaDISM vs Base Aligner (Balanced Filters)', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('balanced_comparison_plot.png', dpi=150, bbox_inches='tight')
plt.savefig('balanced_comparison_plot.pdf', bbox_inches='tight')
print("\nPlot saved to balanced_comparison_plot.png and balanced_comparison_plot.pdf")
