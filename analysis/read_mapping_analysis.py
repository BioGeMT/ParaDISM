import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix
import os
import argparse

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

def calculate_metrics(df):
    all_genes = sorted([g for g in set(df['Origin_Gene'].unique()) | set(df['Mapped_Gene'].unique()) if g != 'NONE'])
    results = []
    
    for gene in all_genes:
        y_true = (df['Origin_Gene'] == gene).astype(int)
        y_pred = (df['Mapped_Gene'] == gene).astype(int)
        
        precision, recall, f1, _ = precision_recall_fscore_support(y_true, y_pred, average='binary', zero_division=0)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
        
        results.append({
            'Gene': gene,
            'TP': tp,
            'FP': fp,
            'FN': fn,
            'Precision': precision,
            'Recall': recall,
            'F1_Score': f1
        })
    
    y_true_multiclass = df['Origin_Gene'].values
    y_pred_multiclass = df['Mapped_Gene'].values
    
    agg_precision, agg_recall, agg_f1, _ = precision_recall_fscore_support(
        y_true_multiclass, y_pred_multiclass, average='weighted', zero_division=0)
    
    overall_accuracy = accuracy_score(y_true_multiclass, y_pred_multiclass)
    
    mapped_df = df[df['Mapped_Gene'] != 'NONE']
    mapped_accuracy = accuracy_score(mapped_df['Origin_Gene'], mapped_df['Mapped_Gene']) if len(mapped_df) > 0 else 0
    
    correct_predictions = (y_true_multiclass == y_pred_multiclass).sum()
    total_mapped = (y_pred_multiclass != 'NONE').sum()
    total_reads = len(y_true_multiclass)
    
    additional_metrics = {
        'correct_predictions': correct_predictions,
        'total_mapped': total_mapped,
        'total_reads': total_reads
    }
    
    return pd.DataFrame(results), agg_precision, agg_recall, agg_f1, overall_accuracy, mapped_accuracy, additional_metrics

COLORS = ['#FF8C00', '#008080']  

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 24,
    'axes.titlesize': 28,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 20
})

def create_aggregated_bar_plot(method1_metrics, method2_metrics, output_path):
    plt.rcParams.update({'font.size': 18, 'font.weight': 'bold', 'axes.titlesize': 22})
    plt.figure(figsize=(10, 6), dpi=300)
    
    metrics = ['Precision', 'Recall', 'Accuracy']
    mapper_values = [method1_metrics['agg_precision'], method1_metrics['agg_recall'], method1_metrics['mapped_accuracy']]
    bowtie2_values = [method2_metrics['agg_precision'], method2_metrics['agg_recall'], method2_metrics['mapped_accuracy']]
    
    bar_width = 0.1
    bar_spacing = 0.12
    group_spacing = 0.15
    x_adjusted = np.linspace(0, len(metrics) * (bar_width + group_spacing), len(metrics))
    
    plt.bar(x_adjusted - bar_spacing/2, mapper_values, width=bar_width,
            label='Mapper', color=COLORS[0])
    plt.bar(x_adjusted + bar_spacing/2, bowtie2_values, width=bar_width,
            label='Bowtie2', color=COLORS[1])
    
    plt.xticks(x_adjusted, metrics, fontsize=15)
    plt.yticks(fontsize=20)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    for i, (v1, v2) in enumerate(zip(mapper_values, bowtie2_values)):
        plt.text(x_adjusted[i] - bar_spacing/2, v1 + 0.03, f"{v1:.4f}", ha="center", fontsize=14)
        plt.text(x_adjusted[i] + bar_spacing/2, v2 + 0.03, f"{v2:.4f}", ha="center", fontsize=14)
    
    plt.legend(fontsize=18, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
    plt.ylim(0, 1.05)
    plt.margins(y=0.05, x=0.01)
    plt.tight_layout()
    
    plt.savefig(f'{output_path}_aggregated_metrics.png', dpi=300)
    plt.close()
    print(f"Bar plot saved: {output_path}_aggregated_metrics.png")

def save_comprehensive_tsv(mapper_per_gene, bowtie2_per_gene, mapper_metrics, bowtie2_metrics, output_path):
    merged = pd.merge(mapper_per_gene, bowtie2_per_gene, on='Gene', suffixes=('_Mapper', '_Bowtie2'))
    
    mapper_correct = mapper_metrics.get('correct_predictions', 0)
    mapper_total_mapped = mapper_metrics.get('total_mapped', 0)
    mapper_total_reads = mapper_metrics.get('total_reads', 100000)
    
    mapper_tp_agg = mapper_correct
    mapper_fp_agg = mapper_total_mapped - mapper_correct
    mapper_fn_agg = mapper_total_reads - mapper_total_mapped
    
    bowtie2_correct = bowtie2_metrics.get('correct_predictions', 0)
    bowtie2_total_mapped = bowtie2_metrics.get('total_mapped', 0)
    bowtie2_total_reads = bowtie2_metrics.get('total_reads', 100000)
    
    bowtie2_tp_agg = bowtie2_correct
    bowtie2_fp_agg = bowtie2_total_mapped - bowtie2_correct
    bowtie2_fn_agg = bowtie2_total_reads - bowtie2_total_mapped
    
    agg_row = pd.DataFrame({
        'Gene': ['AGGREGATED'],
        'TP_Mapper': [mapper_tp_agg],
        'FP_Mapper': [mapper_fp_agg],
        'FN_Mapper': [mapper_fn_agg],
        'Precision_Mapper': [mapper_metrics['agg_precision']],
        'Recall_Mapper': [mapper_metrics['agg_recall']],
        'F1_Score_Mapper': [mapper_metrics['agg_f1']],
        'TP_Bowtie2': [bowtie2_tp_agg],
        'FP_Bowtie2': [bowtie2_fp_agg],
        'FN_Bowtie2': [bowtie2_fn_agg],
        'Precision_Bowtie2': [bowtie2_metrics['agg_precision']],
        'Recall_Bowtie2': [bowtie2_metrics['agg_recall']],
        'F1_Score_Bowtie2': [bowtie2_metrics['agg_f1']]
    })
    
    acc_row = pd.DataFrame({
        'Gene': ['ACCURACY_SUMMARY'],
        'TP_Mapper': [f"Overall: {mapper_metrics['overall_accuracy']:.6f}"],
        'FP_Mapper': [f"Mapped: {mapper_metrics['mapped_accuracy']:.6f}"],
        'FN_Mapper': [''],
        'Precision_Mapper': [''],
        'Recall_Mapper': [''],
        'F1_Score_Mapper': [''],
        'TP_Bowtie2': [f"Overall: {bowtie2_metrics['overall_accuracy']:.6f}"],
        'FP_Bowtie2': [f"Mapped: {bowtie2_metrics['mapped_accuracy']:.6f}"],
        'FN_Bowtie2': [''],
        'Precision_Bowtie2': [''],
        'Recall_Bowtie2': [''],
        'F1_Score_Bowtie2': ['']
    })
    
    complete_data = pd.concat([merged, agg_row, acc_row], ignore_index=True)
    complete_data.to_csv(f'{output_path}_comprehensive_metrics.tsv', sep='\t', index=False)
    print(f"Comprehensive metrics saved: {output_path}_comprehensive_metrics.tsv")

def process_files(tsv_file, fastq_file, sam_file, output_folder, output_prefix):
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_prefix)
    
    df_tsv = pd.read_csv(tsv_file, sep='\t')
    origin_dict = parse_fastq(fastq_file)
    sam_dict = parse_sam(sam_file)
    
    tsv_read_names = set(df_tsv['Read_Name'])
    fastq_read_names = set(origin_dict.keys())
    valid_reads = tsv_read_names.intersection(fastq_read_names)
    
    print(f"Reads in TSV: {len(tsv_read_names)}")
    print(f"Reads in FASTQ: {len(fastq_read_names)}")
    print(f"Valid reads for analysis: {len(valid_reads)}")
    
    mapper_results = []
    for _, row in df_tsv.iterrows():
        read_name = row['Read_Name']
        if read_name in valid_reads:
            mapped_gene = row['Uniquely_Mapped']
            origin_gene = origin_dict[read_name]
            mapper_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 'Mapped_Gene': mapped_gene})
    mapper_df = pd.DataFrame(mapper_results)
    
    bowtie2_results = []
    for read_name in valid_reads:
        origin_gene = origin_dict[read_name]
        mapped_gene = sam_dict.get(read_name, 'NONE')
        bowtie2_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 'Mapped_Gene': mapped_gene})
    bowtie2_df = pd.DataFrame(bowtie2_results)
    
    mapper_cm = pd.crosstab(mapper_df['Origin_Gene'], mapper_df['Mapped_Gene'])
    bowtie2_cm = pd.crosstab(bowtie2_df['Origin_Gene'], bowtie2_df['Mapped_Gene'])
    
    mapper_cm.to_csv(f'{output_path}_mapper_confusion_matrix.tsv', sep='\t')
    bowtie2_cm.to_csv(f'{output_path}_bowtie2_confusion_matrix.tsv', sep='\t')
    
    mapper_per_gene, mapper_agg_prec, mapper_agg_rec, mapper_agg_f1, mapper_overall_acc, mapper_mapped_acc, mapper_additional = calculate_metrics(mapper_df)
    bowtie2_per_gene, bowtie2_agg_prec, bowtie2_agg_rec, bowtie2_agg_f1, bowtie2_overall_acc, bowtie2_mapped_acc, bowtie2_additional = calculate_metrics(bowtie2_df)
    
    mapper_metrics = {
        'agg_precision': mapper_agg_prec,
        'agg_recall': mapper_agg_rec,
        'agg_f1': mapper_agg_f1,
        'overall_accuracy': mapper_overall_acc,
        'mapped_accuracy': mapper_mapped_acc,
        **mapper_additional
    }
    
    bowtie2_metrics = {
        'agg_precision': bowtie2_agg_prec,
        'agg_recall': bowtie2_agg_rec,
        'agg_f1': bowtie2_agg_f1,
        'overall_accuracy': bowtie2_overall_acc,
        'mapped_accuracy': bowtie2_mapped_acc,
        **bowtie2_additional
    }
    
    create_aggregated_bar_plot(mapper_metrics, bowtie2_metrics, output_path)
    save_comprehensive_tsv(mapper_per_gene, bowtie2_per_gene, mapper_metrics, bowtie2_metrics, output_path)
    
    print(f"\n=== AGGREGATED METRICS SUMMARY ===")
    print(f"Mapper - Precision: {mapper_agg_prec:.6f}, Recall: {mapper_agg_rec:.6f}, F1: {mapper_agg_f1:.6f}")
    print(f"Bowtie2 - Precision: {bowtie2_agg_prec:.6f}, Recall: {bowtie2_agg_rec:.6f}, F1: {bowtie2_agg_f1:.6f}")

def main():
    parser = argparse.ArgumentParser(description='Read mapping analysis with aggregated metrics')
    parser.add_argument('-t', '--tsv', required=True, help='TSV file with Mapper results')
    parser.add_argument('-f', '--fastq', required=True, help='FASTQ file with read origins')
    parser.add_argument('-s', '--sam', required=True, help='SAM file with Bowtie2 results')
    parser.add_argument('-o', '--output-dir', default='mapping_results', help='Output directory')
    parser.add_argument('-p', '--prefix', default='read_mapping', help='Output file prefix')
    
    args = parser.parse_args()
    
    process_files(args.tsv, args.fastq, args.sam, args.output_dir, args.prefix)

if __name__ == "__main__":
    main()