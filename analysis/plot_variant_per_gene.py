#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 12
})

COLORS = {'mapper': '#FF8C00', 'bowtie2': '#008080'}

def main():
    ap = argparse.ArgumentParser(description='Plot per-gene SNP precision/recall (horizontal bars)')
    ap.add_argument('--tsv', required=True, help='per-gene confusion TSV')
    ap.add_argument('--output', required=True, help='output PNG path')
    ap.add_argument('--title', default='Variant Calling per Gene (SNP)')
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep='\t')
    # Compute precision/recall
    df['Mapper_Precision'] = df.apply(lambda r: (r['Mapper_TP']/(r['Mapper_TP']+r['Mapper_FP'])) if (r['Mapper_TP']+r['Mapper_FP'])>0 else 0, axis=1)
    df['Mapper_Recall'] = df.apply(lambda r: (r['Mapper_TP']/(r['Mapper_TP']+r['Mapper_FN'])) if (r['Mapper_TP']+r['Mapper_FN'])>0 else 0, axis=1)
    df['Bowtie2_Precision'] = df.apply(lambda r: (r['Bowtie2_TP']/(r['Bowtie2_TP']+r['Bowtie2_FP'])) if (r['Bowtie2_TP']+r['Bowtie2_FP'])>0 else 0, axis=1)
    df['Bowtie2_Recall'] = df.apply(lambda r: (r['Bowtie2_TP']/(r['Bowtie2_TP']+r['Bowtie2_FN'])) if (r['Bowtie2_TP']+r['Bowtie2_FN'])>0 else 0, axis=1)

    genes = df['Gene'].tolist()
    y_pos = np.arange(len(genes))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, max(6, 0.6*len(genes))))

    # Precision bars
    ax1.barh(y_pos - 0.2, df['Mapper_Precision'], 0.4, label='Mapper', color=COLORS['mapper'], alpha=0.85)
    ax1.barh(y_pos + 0.2, df['Bowtie2_Precision'], 0.4, label='Bowtie2', color=COLORS['bowtie2'], alpha=0.85)
    ax1.set_yticks(y_pos); ax1.set_yticklabels(genes)
    ax1.set_xlim(0, 1); ax1.set_xlabel('Precision'); ax1.set_title('SNP Precision'); ax1.grid(axis='x', alpha=0.3)

    # Recall bars
    ax2.barh(y_pos - 0.2, df['Mapper_Recall'], 0.4, label='Mapper', color=COLORS['mapper'], alpha=0.85)
    ax2.barh(y_pos + 0.2, df['Bowtie2_Recall'], 0.4, label='Bowtie2', color=COLORS['bowtie2'], alpha=0.85)
    ax2.set_yticks(y_pos); ax2.set_yticklabels(genes)
    ax2.set_xlim(0, 1); ax2.set_xlabel('Recall'); ax2.set_title('SNP Recall'); ax2.grid(axis='x', alpha=0.3)

    handles, labels = ax1.get_legend_handles_labels(); fig.legend(handles, labels, loc='upper right')
    fig.suptitle(args.title, fontsize=22)
    plt.tight_layout(); os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()

