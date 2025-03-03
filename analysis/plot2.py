import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
import numpy as np
from matplotlib.font_manager import FontProperties

# Define color palette
COLORS = ['#C71585', '#1E90FF']  # Medium Violet Red for Mapper, Bright Blue for Bowtie2

# Set default font sizes
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

def process_files(tsv_file, fastq_file, sam_file, output_prefix):
    # Read and filter TSV
    df_tsv = pd.read_csv(tsv_file, sep='\t')
    mapped_df_tsv = df_tsv[df_tsv['Uniquely_Mapped'] != 'NONE']
    
    # Process FASTQ and SAM
    origin_dict = parse_fastq(fastq_file)
    sam_dict = parse_sam(sam_file)
    
    # Process TSV data
    tsv_results = []
    valid_reads = set()
    for _, row in mapped_df_tsv.iterrows():
        read_name = row['Read_Name']
        mapped_gene = row['Uniquely_Mapped']
        origin_gene = origin_dict.get(read_name, 'Unknown')
        if origin_gene != 'Unknown':
            tsv_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 
                              'Mapped_Gene': mapped_gene})
            valid_reads.add(read_name)
    tsv_df = pd.DataFrame(tsv_results)
    
    # Process SAM data
    sam_results = []
    for read_name in valid_reads:
        mapped_gene = sam_dict.get(read_name)
        if mapped_gene:
            origin_gene = origin_dict[read_name]
            sam_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 
                              'Mapped_Gene': mapped_gene})
    sam_df = pd.DataFrame(sam_results)
    
    # Save mapping TSVs
    tsv_df.to_csv(f'{output_prefix}_tsv_mapping.tsv', sep='\t', index=False)
    sam_df.to_csv(f'{output_prefix}_sam_mapping.tsv', sep='\t', index=False)
    
    # Create confusion matrices
    bold_props = FontProperties(weight='bold', size=16)
    for df, suffix in [(tsv_df, 'tsv'), (sam_df, 'sam')]:
        cm = pd.crosstab(df['Origin_Gene'], df['Mapped_Gene'])
        cm.to_csv(f'{output_prefix}_{suffix}_confusion_matrix.tsv', sep='\t')
        
        plt.figure(figsize=(12, 10))
        ax = sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
        
        # Y-axis
        plt.ylabel('Origin Genes', fontweight='bold')
        ax.set_yticklabels(ax.get_yticklabels(), fontproperties=bold_props)
        
        # X-axis (top, no tilt)
        ax.xaxis.tick_top()
        ax.set_xticklabels(ax.get_xticklabels(), fontproperties=bold_props, rotation=0)
        xlabel = 'Mapper Genes' if suffix == 'tsv' else 'Bowtie2 Genes'
        ax.set_xlabel(xlabel, fontweight='bold', labelpad=15)
        ax.xaxis.set_label_position('top')
        
        plt.savefig(f'{output_prefix}_{suffix}_confusion_matrix.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # Calculate accuracies
    all_genes = sorted(set(tsv_df['Origin_Gene'].unique()) | set(tsv_df['Mapped_Gene'].unique()) |
                      set(sam_df['Mapped_Gene'].unique()))
    
    tsv_accuracies = []
    sam_accuracies = []
    for gene in all_genes:
        tsv_true = (tsv_df['Origin_Gene'] == gene).astype(int)
        tsv_pred = (tsv_df['Mapped_Gene'] == gene).astype(int)
        tsv_accuracies.append(accuracy_score(tsv_true, tsv_pred) if len(tsv_true) > 0 else 0)
        
        sam_true = (sam_df['Origin_Gene'] == gene).astype(int)
        sam_pred = (sam_df['Mapped_Gene'] == gene).astype(int)
        sam_accuracies.append(accuracy_score(sam_true, sam_pred) if len(sam_true) > 0 else 0)
    
    tsv_total_acc = accuracy_score(tsv_df['Origin_Gene'], tsv_df['Mapped_Gene'])
    sam_total_acc = accuracy_score(sam_df['Origin_Gene'], sam_df['Mapped_Gene'])
    all_genes.append('Total')
    tsv_accuracies.append(tsv_total_acc)
    sam_accuracies.append(sam_total_acc)
    
    # Plot accuracies
    plt.rcParams.update({'font.size': 18, 'font.weight': 'bold', 'axes.titlesize': 22})
    plt.figure(figsize=(10, 6), dpi=300)
    
    bar_width = 0.1
    bar_spacing = 0.12
    group_spacing = 0.15
    x_adjusted = np.linspace(0, len(all_genes) * (bar_width + group_spacing), len(all_genes))
    
    plt.bar(x_adjusted - bar_spacing/2, tsv_accuracies, width=bar_width, label='Mapper', color=COLORS[0])
    plt.bar(x_adjusted + bar_spacing/2, sam_accuracies, width=bar_width, label='Bowtie2', color=COLORS[1])
    
    plt.xticks(x_adjusted, all_genes, fontsize=15)
    plt.yticks(fontsize=20)
    plt.ylabel('Accuracy', fontsize=18, fontweight='bold')
    plt.ylim(0.95, 1.0)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    plt.legend(fontsize=18, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
    plt.margins(y=0.1, x=0.01)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_accuracy_comparison.png', dpi=300)
    plt.close()

def main():
    tsv_file = '/Users/nucleotaid/dev/mapper/output/unique_mappings.tsv'
    fastq_file = '/Users/nucleotaid/dev/mapper/simulated_r1.fq'
    sam_file = '/Users/nucleotaid/dev/mapper/filtered.sam'
    output_prefix = 'gene_mapping_results'
    
    process_files(tsv_file, fastq_file, sam_file, output_prefix)
    print(f"Files generated:\n"
          f"- {output_prefix}_tsv_mapping.tsv\n"
          f"- {output_prefix}_sam_mapping.tsv\n"
          f"- {output_prefix}_tsv_confusion_matrix.tsv\n"
          f"- {output_prefix}_sam_confusion_matrix.tsv\n"
          f"- {output_prefix}_tsv_confusion_matrix.png\n"
          f"- {output_prefix}_sam_confusion_matrix.png\n"
          f"- {output_prefix}_accuracy_comparison.png")

if __name__ == "__main__":
    main()