#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 12
})

METHOD_COLORS = {
    'mapper': '#FF8C00',
    'bowtie2': '#008080'
}

def read_metrics(csv_path):
    df = pd.read_csv(csv_path)
    df = df.set_index('Type')
    out = {}
    # Only SNV metrics
    if 'SNV' in df.index:
        out['SNV'] = {
            'precision': float(df.loc['SNV', 'Precision']),
            'recall': float(df.loc['SNV', 'Recall'])
        }
    return out

def main():
    ap = argparse.ArgumentParser(description='Plot VCFCompare metrics per iteration (precision/recall)')
    ap.add_argument('--iteration', type=int, required=True)
    ap.add_argument('--method', required=True, choices=['freebayes', 'bcftools'])
    ap.add_argument('--mapper-csv', required=True)
    ap.add_argument('--bowtie-csv', required=True)
    ap.add_argument('--output-dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    mapper = read_metrics(args.mapper_csv)
    bowtie = read_metrics(args.bowtie_csv)

    labels = ['Mapper', 'Bowtie2']
    x = range(len(labels))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    y_pos = list(x)

    # Precision (SNP only) - horizontal bars
    m_prec = mapper.get('SNV', {}).get('precision', 0)
    b_prec = bowtie.get('SNV', {}).get('precision', 0)
    ax1.barh([p - 0.2 for p in y_pos], [m_prec], height=0.4, label='Mapper', color=METHOD_COLORS['mapper'], alpha=0.85)
    ax1.barh([p + 0.2 for p in y_pos], [b_prec], height=0.4, label='Bowtie2', color=METHOD_COLORS['bowtie2'], alpha=0.85)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels)
    ax1.set_xlim(0, 1)
    ax1.set_title('SNP Precision')
    ax1.grid(axis='x', alpha=0.3)

    # Recall (SNP only) - horizontal bars
    m_rec = mapper.get('SNV', {}).get('recall', 0)
    b_rec = bowtie.get('SNV', {}).get('recall', 0)
    ax2.barh([p - 0.2 for p in y_pos], [m_rec], height=0.4, label='Mapper', color=METHOD_COLORS['mapper'], alpha=0.85)
    ax2.barh([p + 0.2 for p in y_pos], [b_rec], height=0.4, label='Bowtie2', color=METHOD_COLORS['bowtie2'], alpha=0.85)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels)
    ax2.set_xlim(0, 1)
    ax2.set_title('SNP Recall')
    ax2.grid(axis='x', alpha=0.3)

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    fig.suptitle(f'{args.method.capitalize()} Iteration {args.iteration} - VCFCompare SNP Metrics', fontsize=22)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'vcfcompare_metrics.png'), dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()
