#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.metrics import precision_score, recall_score

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
    file_path = f"error_rate_results/err_{rate_str}/analysis/variant_calling_analysis/variant_calling_confusion_matrix.tsv"
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    return pd.read_csv(file_path, sep='\t')

def extract_overall_metrics(error_rates):
    results = []
    
    for rate in error_rates:
        df = load_variant_metrics(rate)
        if df is None:
            continue
            
        # Create y_true and y_pred arrays for sklearn
        y_true_method1 = []
        y_pred_method1 = []
        y_true_method2 = []
        y_pred_method2 = []
        
        for _, row in df.iterrows():
            gene = row['Gene']
            
            # Method 1 - Add TP as correct predictions
            y_true_method1.extend([gene] * int(row['Method1_TP']))
            y_pred_method1.extend([gene] * int(row['Method1_TP']))
            
            # Method 1 - Add FP as incorrect predictions
            y_true_method1.extend(['OTHER'] * int(row['Method1_FP']))
            y_pred_method1.extend([gene] * int(row['Method1_FP']))
            
            # Method 1 - Add FN as missed predictions
            y_true_method1.extend([gene] * int(row['Method1_FN']))
            y_pred_method1.extend(['OTHER'] * int(row['Method1_FN']))
            
            # Method 2 - Same logic
            y_true_method2.extend([gene] * int(row['Method2_TP']))
            y_pred_method2.extend([gene] * int(row['Method2_TP']))
            
            y_true_method2.extend(['OTHER'] * int(row['Method2_FP']))
            y_pred_method2.extend([gene] * int(row['Method2_FP']))
            
            y_true_method2.extend([gene] * int(row['Method2_FN']))
            y_pred_method2.extend(['OTHER'] * int(row['Method2_FN']))
        
        genes = df['Gene'].unique().tolist()
        
        # Calculate overall metrics using weighted average
        precision1 = precision_score(y_true_method1, y_pred_method1, labels=genes, average='weighted', zero_division=0)
        recall1 = recall_score(y_true_method1, y_pred_method1, labels=genes, average='weighted', zero_division=0)
        specificity1 = precision1  # For variant calling, specificity ≈ precision
        
        precision2 = precision_score(y_true_method2, y_pred_method2, labels=genes, average='weighted', zero_division=0)
        recall2 = recall_score(y_true_method2, y_pred_method2, labels=genes, average='weighted', zero_division=0)
        specificity2 = precision2  # For variant calling, specificity ≈ precision
        
        results.append({
            'Error_Rate': rate,
            'Precision_Mapper': precision1,
            'Recall_Mapper': recall1,
            'Specificity_Mapper': specificity1,
            'Precision_Bowtie2': precision2,
            'Recall_Bowtie2': recall2,
            'Specificity_Bowtie2': specificity2
        })
    
    return results

def create_error_rate_comparison_plot():
    error_rates = [0.001, 0.005, 0.010]
    results = extract_overall_metrics(error_rates)
    
    if not results:
        print("No results found for any error rate")
        return
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), dpi=300)
    
    # Error rate labels for y-axis
    error_labels = ['0.1%', '0.5%', '1.0%']
    y_pos = np.arange(len(error_labels))
    bar_height = 0.35
    
    # Plot each metric
    metrics = ['Precision', 'Recall', 'Specificity']
    
    for i, metric in enumerate(metrics):
        mapper_values = df[f'{metric}_Mapper'].values
        bowtie2_values = df[f'{metric}_Bowtie2'].values
        
        axes[i].barh(y_pos - bar_height/2, mapper_values, bar_height, 
                    label='Mapper', color=COLORS[0])
        axes[i].barh(y_pos + bar_height/2, bowtie2_values, bar_height, 
                    label='Bowtie2', color=COLORS[1])
        
        axes[i].set_title(metric, fontweight='bold', fontsize=20)
        axes[i].set_yticks(y_pos)
        axes[i].set_yticklabels(error_labels, fontweight='bold')
        axes[i].set_xlim(0, 1.1)
        axes[i].grid(axis='x', alpha=0.3)
        axes[i].invert_yaxis()
    
    # Add overall title
    fig.suptitle('Variant Calling Performance by Error Rate', fontsize=24, fontweight='bold')
    
    # Add legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    
    plt.tight_layout()
    
    # Save plot
    filename = "plots/variant_calling_error_rate_comparison.png"
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"Created plot: {filename}")
    
    # Print results
    print("\n=== VARIANT CALLING PERFORMANCE BY ERROR RATE ===")
    print("Error Rate | Metric      | Mapper  | Bowtie2")
    print("-" * 45)
    for result in results:
        rate_str = f"{result['Error_Rate']*100:.1f}%"
        print(f"{rate_str:<10} | Precision   | {result['Precision_Mapper']:.3f}   | {result['Precision_Bowtie2']:.3f}")
        print(f"{'':>10} | Recall      | {result['Recall_Mapper']:.3f}   | {result['Recall_Bowtie2']:.3f}")
        print(f"{'':>10} | Specificity | {result['Specificity_Mapper']:.3f}   | {result['Specificity_Bowtie2']:.3f}")
        print("-" * 45)

def main():
    create_error_rate_comparison_plot()

if __name__ == "__main__":
    main()