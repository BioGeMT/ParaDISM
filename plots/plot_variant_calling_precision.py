#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

def extract_variant_calling_precision():
    error_rates = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010]
    
    error_rate_list = []
    method1_precision = []
    method2_precision = []
    
    for rate in error_rates:
        rate_str = f"{rate:.3f}".replace('.', '_')
        file_path = f"error_rate_results/err_{rate_str}/analysis/variant_calling_analysis/comprehensive_variant_metrics.tsv"
        
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path, sep='\t')
                aggregated_row = df[df['Gene'] == 'AGGREGATED']
                
                if not aggregated_row.empty:
                    method1_prec = aggregated_row['Method1_Precision'].iloc[0]
                    method2_prec = aggregated_row['Method2_Precision'].iloc[0]
                    
                    error_rate_list.append(rate)
                    method1_precision.append(method1_prec)
                    method2_precision.append(method2_prec)
                    
                    print(f"Error rate {rate}: Method1 Precision = {method1_prec:.4f}, Method2 Precision = {method2_prec:.4f}")
                else:
                    print(f"Warning: No AGGREGATED row found in {file_path}")
                    
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
        else:
            print(f"File not found: {file_path}")
    
    return error_rate_list, method1_precision, method2_precision

def create_variant_calling_precision_plot():
    error_rates, method1_prec, method2_prec = extract_variant_calling_precision()
    
    if not error_rates:
        print("No data found to plot")
        return
    
    plt.rcParams.update({
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 18,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 14
    })
    
    plt.figure(figsize=(12, 8))
    
    plt.plot(error_rates, method1_prec, 'o-', color='#FF8C00', linewidth=2.5, 
             markersize=8, label='Mapper', markerfacecolor='#FF8C00')
    
    plt.plot(error_rates, method2_prec, 's-', color='#008080', linewidth=2.5, 
             markersize=8, label='Bowtie2', markerfacecolor='#008080')
    
    plt.xlabel('Error Rate', fontsize=16, fontweight='bold')
    plt.title('Variant Calling Precision vs Error Rate', fontsize=18, fontweight='bold')
    
    plt.xlim(0, 0.011)
    plt.ylim(0, 1.05)
    plt.grid(True, alpha=0.3, linestyle='--')
    
    plt.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    
    plt.gca().set_xticks(error_rates)
    plt.gca().set_xticklabels([f'{rate*100:.1f}%' for rate in error_rates])
    
    plt.gca().set_yticks(np.arange(0, 1.1, 0.1))
    plt.gca().set_yticklabels([f'{val*100:.0f}%' for val in np.arange(0, 1.1, 0.1)])
    
    for i, (rate, prec1, prec2) in enumerate(zip(error_rates, method1_prec, method2_prec)):
        plt.annotate(f'{prec1:.3f}', (rate, prec1), textcoords="offset points", 
                    xytext=(0,10), ha='center', fontsize=10, fontweight='bold',
                    color='#FF8C00')
        plt.annotate(f'{prec2:.3f}', (rate, prec2), textcoords="offset points", 
                    xytext=(0,-15), ha='center', fontsize=10, fontweight='bold',
                    color='#008080')
    
    plt.tight_layout()
    
    output_file = 'variant_calling_precision_vs_error_rate.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Variant calling precision plot saved: {output_file}")
    
    print("\n" + "="*60)
    print("VARIANT CALLING PRECISION SUMMARY")
    print("="*60)
    print(f"Error Rate Range: {min(error_rates):.1%} - {max(error_rates):.1%}")
    print(f"Custom Mapper Precision Range: {min(method1_prec):.3f} - {max(method1_prec):.3f}")
    print(f"Bowtie2 Precision Range: {min(method2_prec):.3f} - {max(method2_prec):.3f}")
    print(f"Average Custom Mapper Precision: {np.mean(method1_prec):.3f}")
    print(f"Average Bowtie2 Precision: {np.mean(method2_prec):.3f}")
    
    custom_mapper_drop = max(method1_prec) - min(method1_prec)
    bowtie2_drop = max(method2_prec) - min(method2_prec)
    print(f"Custom Mapper Precision Drop: {custom_mapper_drop:.3f}")
    print(f"Bowtie2 Precision Drop: {bowtie2_drop:.3f}")
    
    plt.show()
    
    return error_rates, method1_prec, method2_prec

if __name__ == "__main__":
    create_variant_calling_precision_plot()