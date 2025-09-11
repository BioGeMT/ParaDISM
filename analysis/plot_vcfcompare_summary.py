#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 12
})

COLORS = {
    'mapper_snv': '#FF8C00',
    'bowtie_snv': '#008080',
    'mapper_indel': '#FFA64D',
    'bowtie_indel': '#33A6A6'
}

def read_metrics(csv_path):
    df = pd.read_csv(csv_path)
    df = df.set_index('Type')
    out = {}
    for t in ['SNV', 'INDEL']:
        if t in df.index:
            out[t] = {
                'precision': float(df.loc[t, 'Precision']),
                'recall': float(df.loc[t, 'Recall'])
            }
    return out

def main():
    ap = argparse.ArgumentParser(description='VCFCompare summary across iterations')
    ap.add_argument('--method', required=True, choices=['freebayes', 'bcftools'])
    ap.add_argument('--iterations', type=int, required=True)
    ap.add_argument('--root-dir', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    iters = []
    mapper_snv_p, mapper_snv_r = [], []
    bowtie_snv_p, bowtie_snv_r = [], []

    for i in range(1, args.iterations + 1):
        iter_dir = os.path.join(args.root_dir, f'iter_{i}')
        mapper_csv = os.path.join(iter_dir, 'vcfcompare_mapper.csv')
        bowtie_csv = os.path.join(iter_dir, 'vcfcompare_bowtie2.csv')
        if not (os.path.exists(mapper_csv) and os.path.exists(bowtie_csv)):
            continue
        m = read_metrics(mapper_csv)
        b = read_metrics(bowtie_csv)
        iters.append(i)
        mapper_snv_p.append(m.get('SNV', {}).get('precision', 0))
        mapper_snv_r.append(m.get('SNV', {}).get('recall', 0))
        bowtie_snv_p.append(b.get('SNV', {}).get('precision', 0))
        bowtie_snv_r.append(b.get('SNV', {}).get('recall', 0))

    if not iters:
        print('No iterations with VCFCompare CSVs found.')
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Precision
    ax1.plot(iters, mapper_snv_p, 'o-', color=COLORS['mapper_snv'], label='Mapper SNP', linewidth=3)
    ax1.plot(iters, bowtie_snv_p, 's-', color=COLORS['bowtie_snv'], label='Bowtie2 SNP', linewidth=3)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Precision')
    ax1.set_title('Precision across iterations')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.2)
    ax1.set_xticks(iters)
    ax1.set_yticks(np.linspace(0, 1, 6))

    # Recall
    ax2.plot(iters, mapper_snv_r, 'o-', color=COLORS['mapper_snv'], label='Mapper SNP', linewidth=3)
    ax2.plot(iters, bowtie_snv_r, 's-', color=COLORS['bowtie_snv'], label='Bowtie2 SNP', linewidth=3)
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('Recall')
    ax2.set_title('Recall across iterations')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1.2)
    ax2.set_xticks(iters)
    ax2.set_yticks(np.linspace(0, 1, 6))

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    fig.suptitle(f'{args.method.capitalize()} - VCFCompare Summary', fontsize=22)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()
