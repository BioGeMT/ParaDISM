#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

COLORS = ['#FF8C00', '#008080']  

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 24,
    'axes.titlesize': 28,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 14
})

def load_variant_metrics(error_rate):
    rate_str = f"{error_rate:.3f}".replace('.', '_')
    file_path = f"error_rate_results/err_{rate_str}/analysis/variant_calling_analysis/comprehensive_variant_metrics.tsv"
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    return pd.read_csv(file_path, sep='\t')

def calculate_variant_metrics(df):
    results = []
    
    gene_rows = df[df['Gene'] != 'AGGREGATED'].copy()
    
    for _, row in gene_rows.iterrows():
        # for variant calling, specificity = 1 - (FP / total_predicted) = precision
        tp1 = row['Method1_TP']
        fp1 = row['Method1_FP']
        total_predicted1 = tp1 + fp1
        specificity1 = 1.0 - (fp1 / total_predicted1) if total_predicted1 > 0 else 1.0
        
        tp2 = row['Method2_TP']
        fp2 = row['Method2_FP']
        total_predicted2 = tp2 + fp2
        specificity2 = 1.0 - (fp2 / total_predicted2) if total_predicted2 > 0 else 1.0
        
        results.append({
            'Gene': row['Gene'],
            'Precision': row['Method1_Precision'],
            'Recall': row['Method1_Recall'],
            'Specificity': specificity1,
            'Precision_Bowtie2': row['Method2_Precision'],
            'Recall_Bowtie2': row['Method2_Recall'],
            'Specificity_Bowtie2': specificity2
        })
    
    agg_row = df[df['Gene'] == 'AGGREGATED'].iloc[0]
    
    tp1_overall = agg_row['Method1_TP']
    fp1_overall = agg_row['Method1_FP']
    total_predicted1_overall = tp1_overall + fp1_overall
    specificity1_overall = 1.0 - (fp1_overall / total_predicted1_overall) if total_predicted1_overall > 0 else 1.0
    
    tp2_overall = agg_row['Method2_TP']
    fp2_overall = agg_row['Method2_FP']
    total_predicted2_overall = tp2_overall + fp2_overall
    specificity2_overall = 1.0 - (fp2_overall / total_predicted2_overall) if total_predicted2_overall > 0 else 1.0
    
    results.append({
        'Gene': 'Overall',
        'Precision': agg_row['Method1_Precision'],
        'Recall': agg_row['Method1_Recall'],
        'Specificity': specificity1_overall,
        'Precision_Bowtie2': agg_row['Method2_Precision'],
        'Recall_Bowtie2': agg_row['Method2_Recall'],
        'Specificity_Bowtie2': specificity2_overall
    })
    
    return results

def create_variant_plot(error_rate):
    df = load_variant_metrics(error_rate)
    
    if df is None:
        print(f"Could not load variant metrics for error rate {error_rate}")
        return
    
    print(f"\n=== Processing Variant Calling Error Rate {error_rate} ===")
    
    results = calculate_variant_metrics(df)
    
    combined_results = []
    for result in results:
        combined_results.append({
            'Gene': result['Gene'],
            'Precision_Mapper': result['Precision'],
            'Recall_Mapper': result['Recall'],
            'Specificity_Mapper': result['Specificity'],
            'Precision_Bowtie2': result['Precision_Bowtie2'],
            'Recall_Bowtie2': result['Recall_Bowtie2'],
            'Specificity_Bowtie2': result['Specificity_Bowtie2']
        })
    
    create_horizontal_plot(combined_results, f"Variant Calling Performance - Error Rate {error_rate*100:.1f}%", 
                          f"variant_calling_{str(error_rate).replace('.', '_')}.png")
    save_tsv_output(combined_results, error_rate)
    
    print(f"\n=== VARIANT CALLING RESULTS FOR ERROR RATE {error_rate} ===")
    for result in combined_results:
        print(f"{result['Gene']:<12} - Mapper: P={result['Precision_Mapper']:.3f}, R={result['Recall_Mapper']:.3f}, S={result['Specificity_Mapper']:.3f}")
        print(f"{'':>12}   Bowtie2: P={result['Precision_Bowtie2']:.3f}, R={result['Recall_Bowtie2']:.3f}, S={result['Specificity_Bowtie2']:.3f}")

def create_horizontal_plot(results, title, filename):
    metrics_df = pd.DataFrame(results)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=300)
    
    genes = metrics_df['Gene'].tolist()
    y_pos = np.arange(len(genes))
    bar_height = 0.35
    
    for i, metric in enumerate(['Precision', 'Recall', 'Specificity']):
        axes[i].barh(y_pos - bar_height/2, metrics_df[f'{metric}_Mapper'], bar_height, 
                    label='Mapper', color=COLORS[0])
        axes[i].barh(y_pos + bar_height/2, metrics_df[f'{metric}_Bowtie2'], bar_height, 
                    label='Bowtie2', color=COLORS[1])
        axes[i].set_title(metric, fontweight='bold', fontsize=20)
        axes[i].set_yticks(y_pos)
        axes[i].set_yticklabels(genes, fontweight='bold')
        axes[i].set_xlim(0, 1.1)
        axes[i].grid(axis='x', alpha=0.3)
        axes[i].invert_yaxis()
    
    fig.suptitle(title, fontsize=24, fontweight='bold')
    
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    
    plt.tight_layout()
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"Created plot: {filename}")

def save_tsv_output(results, error_rate):
    df = pd.DataFrame(results)
    filename = f"variant_calling_{str(error_rate).replace('.', '_')}.tsv"
    df.to_csv(filename, sep='\t', index=False)
    print(f"Created TSV: {filename}")

def main():
    error_rates = [0.001, 0.005, 0.010]
    for rate in error_rates:
        create_variant_plot(rate)

if __name__ == "__main__":
    main()