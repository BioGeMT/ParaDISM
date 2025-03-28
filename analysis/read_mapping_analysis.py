import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
import numpy as np
from matplotlib.font_manager import FontProperties
import os
import argparse


COLORS = ['#C71585', '#1E90FF']  


plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 24,
    'axes.titlesize': 28,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 20
})

def parse_fastq(fastq_file):
    read_to_origin = {}
    with open(fastq_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                parts = line.strip().split()
                read_name = parts[0][1:]
                gene_info = parts[1].split(';')[0]
                read_to_origin[read_name] = gene_info
    return read_to_origin

def parse_sam(sam_file):
    read_to_mapped = {}
    with open(sam_file, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                parts = line.strip().split('\t')
                read_name = parts[0]
                mapped_gene = parts[2]
                if mapped_gene != '*':
                    read_to_mapped[read_name] = mapped_gene
    return read_to_mapped

def process_files(tsv_file, fastq_file, sam_file, output_folder, output_prefix):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Create full path for output files
    output_path = os.path.join(output_folder, output_prefix)
    
    # Read TSV without filtering NONE
    df_tsv = pd.read_csv(tsv_file, sep='\t')
    
    # Process FASTQ and SAM
    origin_dict = parse_fastq(fastq_file)
    sam_dict = parse_sam(sam_file)
    
    # Get the set of all read names from each file
    tsv_read_names = set(df_tsv['Read_Name'])
    fastq_read_names = set(origin_dict.keys())
    sam_read_names = set(sam_dict.keys())
    
    # Print counts to verify
    print(f"Reads in TSV: {len(tsv_read_names)}")
    print(f"Reads in FASTQ: {len(fastq_read_names)}")
    print(f"Reads in SAM: {len(sam_read_names)}")
    
    # Find reads that exist in all three files
    common_reads = tsv_read_names.intersection(fastq_read_names).intersection(sam_read_names)
    print(f"Common reads in all files: {len(common_reads)}")
    
    # Find all reads in both TSV and FASTQ (needed for analysis)
    valid_reads = tsv_read_names.intersection(fastq_read_names)
    print(f"Reads in both TSV and FASTQ: {len(valid_reads)}")
    
    # Process TSV data
    tsv_results = []
    for _, row in df_tsv.iterrows():
        read_name = row['Read_Name']
        if read_name in valid_reads:
            mapped_gene = row['Uniquely_Mapped']
            origin_gene = origin_dict[read_name]
            tsv_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 
                              'Mapped_Gene': mapped_gene})
    tsv_df = pd.DataFrame(tsv_results)
    
    # Process SAM data
    sam_results = []
    for read_name in valid_reads:
        origin_gene = origin_dict[read_name]
        # Use 'NONE' if the read is not in the SAM dictionary
        mapped_gene = sam_dict.get(read_name, 'NONE')
        sam_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 
                          'Mapped_Gene': mapped_gene})
    sam_df = pd.DataFrame(sam_results)
    
    # Save mapping TSVs
    tsv_df.to_csv(f'{output_path}_tsv_mapping.tsv', sep='\t', index=False)
    sam_df.to_csv(f'{output_path}_sam_mapping.tsv', sep='\t', index=False)
    
    # Create confusion matrices
    bold_props = FontProperties(weight='bold', size=16)
    for df, suffix in [(tsv_df, 'tsv'), (sam_df, 'sam')]:
        cm = pd.crosstab(df['Origin_Gene'], df['Mapped_Gene'])
        cm.to_csv(f'{output_path}_{suffix}_confusion_matrix.tsv', sep='\t')
        
        plt.figure(figsize=(12, 10))
        ax = sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
        
        # Y-axis
        plt.ylabel('Origin Genes', fontweight='bold')
        ax.set_yticklabels(ax.get_yticklabels(), fontproperties=bold_props)
        
    
        ax.xaxis.tick_top()
        if suffix == 'tsv':
            mapper_props = FontProperties(weight='bold', size=14)  
            ax.set_xticklabels(ax.get_xticklabels(), fontproperties=mapper_props, rotation=0)
        else:
            ax.set_xticklabels(ax.get_xticklabels(), fontproperties=bold_props, rotation=0)
        xlabel = 'Mapper Genes' if suffix == 'tsv' else 'Bowtie2 Genes'
        ax.set_xlabel(xlabel, fontweight='bold', labelpad=15)
        ax.xaxis.set_label_position('top')
        
        plt.savefig(f'{output_path}_{suffix}_confusion_matrix.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Calculate accuracies - treating NONE as any other category
    # Filter out 'NONE' from the list of genes for the plot
    all_genes = sorted(set(tsv_df['Origin_Gene'].unique()) | set(tsv_df['Mapped_Gene'].unique()) |
                      set(sam_df['Mapped_Gene'].unique()))
    
    # Keep 'NONE' in calculations but remove it from display
    display_genes = [gene for gene in all_genes if gene != 'NONE']
    
    tsv_accuracies = []
    sam_accuracies = []
    # Calculate accuracies for all genes (including NONE) for calculations
    for gene in all_genes:
        tsv_true = (tsv_df['Origin_Gene'] == gene).astype(int)
        tsv_pred = (tsv_df['Mapped_Gene'] == gene).astype(int)
        tsv_acc = accuracy_score(tsv_true, tsv_pred) if len(tsv_true) > 0 else 0
        
        sam_true = (sam_df['Origin_Gene'] == gene).astype(int)
        sam_pred = (sam_df['Mapped_Gene'] == gene).astype(int)
        sam_acc = accuracy_score(sam_true, sam_pred) if len(sam_true) > 0 else 0
        
        # Only add to displayed accuracies if not 'NONE'
        if gene != 'NONE':
            tsv_accuracies.append(tsv_acc)
            sam_accuracies.append(sam_acc)
    
    tsv_total_acc = accuracy_score(tsv_df['Origin_Gene'], tsv_df['Mapped_Gene'])
    sam_total_acc = accuracy_score(sam_df['Origin_Gene'], sam_df['Mapped_Gene'])
    display_genes.append('Total')
    tsv_accuracies.append(tsv_total_acc)
    sam_accuracies.append(sam_total_acc)
    
    # Plot accuracies
    plt.rcParams.update({'font.size': 18, 'font.weight': 'bold', 'axes.titlesize': 22})
    plt.figure(figsize=(10, 6), dpi=300)
    
    bar_width = 0.1
    bar_spacing = 0.12
    group_spacing = 0.15
    x_adjusted = np.linspace(0, len(display_genes) * (bar_width + group_spacing), len(display_genes))
    
    plt.bar(x_adjusted - bar_spacing/2, tsv_accuracies, width=bar_width, label='Mapper', color=COLORS[0])
    plt.bar(x_adjusted + bar_spacing/2, sam_accuracies, width=bar_width, label='Bowtie2', color=COLORS[1])
    
    plt.xticks(x_adjusted, display_genes, fontsize=13)
    plt.yticks(fontsize=20)
    plt.ylabel('Accuracy', fontsize=18, fontweight='bold')
    plt.ylim(min(min(tsv_accuracies), min(sam_accuracies)) * 0.95, 1.0)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.legend(fontsize=18, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
    plt.margins(y=0.1, x=0.01)
    plt.ylim(0.1, 1.0)
    plt.title('Accuracy for All Reads', fontsize=22, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_path}_accuracy_comparison.png', dpi=300)
    plt.close()
    
    # Additional plot: accuracy comparison for only reads that Mapper mapped (not NONE)
    mapped_reads_tsv = tsv_df[tsv_df['Mapped_Gene'] != 'NONE']
    
    # Get the corresponding reads from SAM data
    mapped_reads_sam = sam_df[sam_df['Read_Name'].isin(mapped_reads_tsv['Read_Name'])]
    
    # Get the set of genes present in the filtered datasets
    mapped_genes = sorted(set(mapped_reads_tsv['Origin_Gene'].unique()) | 
                         set(mapped_reads_tsv['Mapped_Gene'].unique()) |
                         set(mapped_reads_sam['Mapped_Gene'].unique()))
    
    # Calculate accuracies
    tsv_mapped_accuracies = []
    sam_mapped_accuracies = []
    for gene in mapped_genes:
        tsv_true = (mapped_reads_tsv['Origin_Gene'] == gene).astype(int)
        tsv_pred = (mapped_reads_tsv['Mapped_Gene'] == gene).astype(int)
        tsv_mapped_accuracies.append(accuracy_score(tsv_true, tsv_pred) if len(tsv_true) > 0 else 0)
        
        sam_true = (mapped_reads_sam['Origin_Gene'] == gene).astype(int)
        sam_pred = (mapped_reads_sam['Mapped_Gene'] == gene).astype(int)
        sam_mapped_accuracies.append(accuracy_score(sam_true, sam_pred) if len(sam_true) > 0 else 0)
    
    # Total accuracy
    tsv_mapped_total_acc = accuracy_score(mapped_reads_tsv['Origin_Gene'], mapped_reads_tsv['Mapped_Gene'])
    sam_mapped_total_acc = accuracy_score(mapped_reads_sam['Origin_Gene'], mapped_reads_sam['Mapped_Gene'])
    mapped_genes.append('Total')
    tsv_mapped_accuracies.append(tsv_mapped_total_acc)
    sam_mapped_accuracies.append(sam_mapped_total_acc)
    
    # Create plot
    plt.figure(figsize=(10, 6), dpi=300)
    x_mapped = np.linspace(0, len(mapped_genes) * (bar_width + group_spacing), len(mapped_genes))
    
    plt.bar(x_mapped - bar_spacing/2, tsv_mapped_accuracies, width=bar_width, label='Mapper', color=COLORS[0])
    plt.bar(x_mapped + bar_spacing/2, sam_mapped_accuracies, width=bar_width, label='Bowtie2', color=COLORS[1])
    
    plt.xticks(x_mapped, mapped_genes, fontsize=14)
    plt.yticks(fontsize=20)
    plt.ylabel('Accuracy', fontsize=18, fontweight='bold')
    plt.ylim(min(min(tsv_mapped_accuracies), min(sam_mapped_accuracies)) * 0.95, 1.0)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.legend(fontsize=18, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
    plt.title('Accuracy for Mapper-Mapped Reads Only', fontsize=22, fontweight='bold')
    plt.margins(y=0.1, x=0.01)
    plt.tight_layout()
    plt.savefig(f'{output_path}_mapped_only_accuracy_comparison.png', dpi=300)
    plt.close()
    
    # Return list of generated files
    generated_files = [
        f"{output_path}_tsv_mapping.tsv",
        f"{output_path}_sam_mapping.tsv",
        f"{output_path}_tsv_confusion_matrix.tsv",
        f"{output_path}_sam_confusion_matrix.tsv",
        f"{output_path}_tsv_confusion_matrix.png",
        f"{output_path}_sam_confusion_matrix.png",
        f"{output_path}_accuracy_comparison.png",
        f"{output_path}_mapped_only_accuracy_comparison.png"
    ]
    
    return generated_files

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Analyze gene mapping performance comparing Mapper and Bowtie2.')
    
    # Required arguments
    parser.add_argument('-t', '--tsv', required=True, help='Path to the TSV file with uniquely mapped reads')
    parser.add_argument('-f', '--fastq', required=True, help='Path to the FASTQ file with read origins')
    parser.add_argument('-s', '--sam', required=True, help='Path to the SAM file with mapped reads')
    
    # Optional arguments
    parser.add_argument('-o', '--output-dir', default='read_mapping_analysis', 
                        help='Output directory for analysis results (default: read_mapping_analysis)')
    parser.add_argument('-p', '--prefix', default='read_mapping_results',
                        help='Prefix for output files (default: read_mapping_results)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process files
    generated_files = process_files(
        args.tsv,
        args.fastq,
        args.sam,
        args.output_dir,
        args.prefix
    )
    
    print(f"\nFiles generated in folder '{args.output_dir}':")
    for file_path in generated_files:
        print(f"- {os.path.basename(file_path)}")

if __name__ == "__main__":
    main()