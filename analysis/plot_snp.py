#!/usr/bin/env python3
"""
Simple SNP calling performance visualization script (plot_snp.py)

Usage:
    python plot_snp.py <output_dir> <mapper_file.tsv> <bowtie2_file.tsv>
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Brighter blue and darker magenta
COLORS = ['#C71585', '#1E90FF']  # Medium Violet Red (darker magenta), Bright Blue (Dodger Blue)

def parse_tsv(file_path):
    """Parse TSV file into DataFrame and calculate accuracy if needed."""
    df = pd.read_csv(file_path, sep='\t')
    
    # Standardize Gene/Chromosome column name
    for col in df.columns:
        if col.lower() in ['chromosome', 'gene']:
            if col != 'Gene':
                df = df.rename(columns={col: 'Gene'})
            break
    
    # Calculate Accuracy if not present
    if 'Accuracy' not in df.columns:
        df['Accuracy'] = df['True_Positives'] / (df['True_Positives'] + df['False_Positives'] + df['False_Negatives'])
    
    return df

def generate_plot(datasets, names, output_path):
    """Generate a simple bar chart of performance metrics."""
    # Extract total metrics from each dataset
    total_metrics = []
    for i, df in enumerate(datasets):
        total_row = df[df['Gene'] == 'Total'].iloc[0]
        total_metrics.append({
            'Dataset': names[i],
            'Precision': total_row['Precision'],
            'Recall': total_row['Recall'],
            'F1_Score': total_row['F1_Score'],
            'Accuracy': total_row['Accuracy']
        })
    
    metrics_df = pd.DataFrame(total_metrics)
    
    # Set up the plot with much larger, bold fonts
    plt.rcParams.update({
        'font.size': 18,           # Base font size
        'font.weight': 'bold',     # Bold font for all text
        'axes.titleweight': 'bold',
        'axes.labelweight': 'bold',
        'axes.titlesize': 22
    })
    
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    
    # Define bar positions - final adjustment for 2 datasets
    bar_width = 0.12  # Maintain the current width
    bar_spacing = 0.14  # Space between bars within a group
    x = np.array([0, 0.4, 0.8, 1.2])  # Less distance between x-axis ticks
    
    # Plot bars for each dataset
    for i, dataset in enumerate(names):
        ax.bar(
            x + i * bar_spacing,  # Spacing between bars in the same group
            [metrics_df.loc[i, 'Precision'], metrics_df.loc[i, 'Recall'], 
             metrics_df.loc[i, 'F1_Score'], metrics_df.loc[i, 'Accuracy']], 
            width=bar_width, 
            label=dataset,
            color=COLORS[i]
        )
    
    # Configure plot
    ax.set_title('SNP Calling Performance', fontsize=22, fontweight='bold')
    ax.set_xticks(x + (bar_spacing/2))  # Center ticks between the two datasets
    ax.set_xticklabels(['Precision', 'Recall', 'F1 Score', 'Accuracy'], fontsize=20, fontweight='bold')
    ax.tick_params(axis='y', labelsize=20)
    ax.set_ylim(0, 0.5)
    
    # Position the legend with larger, bold font
    ax.legend(fontsize=18, prop={'weight': 'bold'})
    
    # Adjust margins
    ax.margins(y=0.1, x=0.01)
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)

def generate_metrics_tsv(datasets, names, output_path):
    """Generate consolidated metrics TSV."""
    dfs_with_dataset = []
    for df, name in zip(datasets, names):
        df_copy = df.copy()
        df_copy['Dataset'] = name
        dfs_with_dataset.append(df_copy)
    
    combined_df = pd.concat(dfs_with_dataset)
    cols = ['Gene', 'Dataset'] + [col for col in combined_df.columns if col not in ['Gene', 'Dataset']]
    combined_df = combined_df[cols]
    combined_df.to_csv(output_path, sep='\t', index=False, float_format='%.3f')

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <output_dir> <mapper_file.tsv> <bowtie2_file.tsv>")
        sys.exit(1)
    
    output_dir = Path(sys.argv[1])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Fixed dataset names with Mapper as the first dataset
    dataset_names = ["Mapper", "Bowtie2"]
    
    # Parse input files
    try:
        datasets = [
            parse_tsv(sys.argv[2]),  # Mapper (first input)
            parse_tsv(sys.argv[3])   # Bowtie2 (second input)
        ]
    except Exception as e:
        print(f"Error processing files: {e}")
        sys.exit(1)
    
    # Generate outputs
    generate_plot(datasets, dataset_names, output_dir / 'snp_calling_performance.png')
    generate_metrics_tsv(datasets, dataset_names, output_dir / 'all_metrics.tsv')
    print(f"Files saved to {output_dir}")

if __name__ == '__main__':
    main()