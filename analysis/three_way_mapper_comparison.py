import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix
import os
import argparse
import csv

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
    """Calculate precision, recall, and accuracy for each gene using sklearn"""
    all_genes = sorted([g for g in set(df['Origin_Gene'].unique()) | set(df['Mapped_Gene'].unique()) if g != 'NONE'])
    
    results = []
    
    # Per-gene metrics
    for gene in all_genes:
        # Create binary labels for this gene (one-vs-rest)
        y_true = (df['Origin_Gene'] == gene).astype(int)
        y_pred = (df['Mapped_Gene'] == gene).astype(int)
        
        # Calculate metrics using sklearn
        precision, recall, f1, _ = precision_recall_fscore_support(y_true, y_pred, average='binary', zero_division=0)
        
        # Get TP, FP, FN from confusion matrix
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
    
    # Use multiclass approach for aggregated metrics
    y_true_multiclass = df['Origin_Gene'].values
    y_pred_multiclass = df['Mapped_Gene'].values
    
    # Calculate aggregated metrics using weighted average
    agg_precision, agg_recall, agg_f1, _ = precision_recall_fscore_support(
        y_true_multiclass, y_pred_multiclass, average='weighted', zero_division=0)
    
    # Overall accuracy: all correct predictions / total reads
    overall_accuracy = accuracy_score(y_true_multiclass, y_pred_multiclass)
    
    # Mapped reads accuracy: correctly mapped reads / total mapped reads (excluding NONE)
    mapped_df = df[df['Mapped_Gene'] != 'NONE']
    mapped_accuracy = accuracy_score(mapped_df['Origin_Gene'], mapped_df['Mapped_Gene']) if len(mapped_df) > 0 else 0
    
    # Calculate correct predictions and total mapped for aggregated TP/FP/FN
    correct_predictions = (y_true_multiclass == y_pred_multiclass).sum()
    total_mapped = (y_pred_multiclass != 'NONE').sum()
    total_reads = len(y_true_multiclass)
    
    # Return additional metrics
    additional_metrics = {
        'correct_predictions': correct_predictions,
        'total_mapped': total_mapped,
        'total_reads': total_reads
    }
    
    return pd.DataFrame(results), agg_precision, agg_recall, agg_f1, overall_accuracy, mapped_accuracy, additional_metrics

COLORS = ['#FF8C00', '#008080', '#9932CC']  # Orange, Teal, Purple

plt.rcParams.update({
    'font.size': 18,
    'axes.labelsize': 20,
    'axes.titlesize': 24,
    'xtick.labelsize': 14,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

def create_three_way_bar_plot(mapper_metrics, mapper_snp_metrics, bowtie2_metrics, output_path):
    """Create three-way comparison bar plot"""
    
    plt.figure(figsize=(12, 7), dpi=300)
    
    # Data for the 3 groups
    metrics = ['Precision', 'Recall', 'Accuracy']
    mapper_values = [mapper_metrics['agg_precision'], mapper_metrics['agg_recall'], mapper_metrics['mapped_accuracy']]
    mapper_snp_values = [mapper_snp_metrics['agg_precision'], mapper_snp_metrics['agg_recall'], mapper_snp_metrics['mapped_accuracy']]
    bowtie2_values = [bowtie2_metrics['agg_precision'], bowtie2_metrics['agg_recall'], bowtie2_metrics['mapped_accuracy']]
    
    # Create bar positions
    bar_width = 0.08
    bar_spacing = 0.1
    group_spacing = 0.2
    x_adjusted = np.linspace(0, len(metrics) * (bar_width + group_spacing), len(metrics))
    
    # Create bars
    plt.bar(x_adjusted - bar_spacing, mapper_values, width=bar_width,
            label='Mapper (Full)', color=COLORS[0])
    plt.bar(x_adjusted, mapper_snp_values, width=bar_width,
            label='Mapper (SNP-only)', color=COLORS[2])
    plt.bar(x_adjusted + bar_spacing, bowtie2_values, width=bar_width,
            label='Bowtie2', color=COLORS[1])
    
    # Formatting
    plt.xticks(x_adjusted, metrics, fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on bars
    for i, (v1, v2, v3) in enumerate(zip(mapper_values, mapper_snp_values, bowtie2_values)):
        plt.text(x_adjusted[i] - bar_spacing, v1 + 0.02, f"{v1:.4f}", ha="center", fontsize=12)
        plt.text(x_adjusted[i], v2 + 0.02, f"{v2:.4f}", ha="center", fontsize=12)
        plt.text(x_adjusted[i] + bar_spacing, v3 + 0.02, f"{v3:.4f}", ha="center", fontsize=12)
    
    plt.legend(fontsize=14, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)
    plt.ylim(0, 1.1)
    plt.margins(y=0.1, x=0.02)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'{output_path}_three_way_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Three-way comparison plot saved: {output_path}_three_way_comparison.png")

def save_three_way_comprehensive_tsv(mapper_per_gene, mapper_snp_per_gene, bowtie2_per_gene, 
                                    mapper_metrics, mapper_snp_metrics, bowtie2_metrics, output_path):
    """Save comprehensive TSV with all three methods"""
    
    # Merge per-gene data
    merged = pd.merge(mapper_per_gene, mapper_snp_per_gene, on='Gene', suffixes=('_Mapper_Full', '_Mapper_SNP'))
    merged = pd.merge(merged, bowtie2_per_gene, on='Gene')
    merged = merged.rename(columns={col: col + '_Bowtie2' if not col.endswith(('_Mapper_Full', '_Mapper_SNP')) and col != 'Gene' else col for col in merged.columns})
    
    # Calculate aggregated metrics for each method
    def get_aggregated_metrics(metrics_dict):
        correct = metrics_dict.get('correct_predictions', 0)
        total_mapped = metrics_dict.get('total_mapped', 0)
        total_reads = metrics_dict.get('total_reads', 100000)
        
        tp_agg = correct
        fp_agg = total_mapped - correct
        fn_agg = total_reads - total_mapped
        
        return tp_agg, fp_agg, fn_agg
    
    mapper_tp, mapper_fp, mapper_fn = get_aggregated_metrics(mapper_metrics)
    mapper_snp_tp, mapper_snp_fp, mapper_snp_fn = get_aggregated_metrics(mapper_snp_metrics)
    bowtie2_tp, bowtie2_fp, bowtie2_fn = get_aggregated_metrics(bowtie2_metrics)
    
    # Add aggregated summary rows
    agg_row = pd.DataFrame({
        'Gene': ['AGGREGATED'],
        'TP_Mapper_Full': [mapper_tp],
        'FP_Mapper_Full': [mapper_fp],
        'FN_Mapper_Full': [mapper_fn],
        'Precision_Mapper_Full': [mapper_metrics['agg_precision']],
        'Recall_Mapper_Full': [mapper_metrics['agg_recall']],
        'F1_Score_Mapper_Full': [mapper_metrics['agg_f1']],
        'TP_Mapper_SNP': [mapper_snp_tp],
        'FP_Mapper_SNP': [mapper_snp_fp],
        'FN_Mapper_SNP': [mapper_snp_fn],
        'Precision_Mapper_SNP': [mapper_snp_metrics['agg_precision']],
        'Recall_Mapper_SNP': [mapper_snp_metrics['agg_recall']],
        'F1_Score_Mapper_SNP': [mapper_snp_metrics['agg_f1']],
        'TP_Bowtie2': [bowtie2_tp],
        'FP_Bowtie2': [bowtie2_fp],
        'FN_Bowtie2': [bowtie2_fn],
        'Precision_Bowtie2': [bowtie2_metrics['agg_precision']],
        'Recall_Bowtie2': [bowtie2_metrics['agg_recall']],
        'F1_Score_Bowtie2': [bowtie2_metrics['agg_f1']]
    })
    
    # Accuracy summary
    acc_row = pd.DataFrame({
        'Gene': ['ACCURACY_SUMMARY'],
        'TP_Mapper_Full': [f"Overall: {mapper_metrics['overall_accuracy']:.6f}"],
        'FP_Mapper_Full': [f"Mapped: {mapper_metrics['mapped_accuracy']:.6f}"],
        'FN_Mapper_Full': [''],
        'Precision_Mapper_Full': [''],
        'Recall_Mapper_Full': [''],
        'F1_Score_Mapper_Full': [''],
        'TP_Mapper_SNP': [f"Overall: {mapper_snp_metrics['overall_accuracy']:.6f}"],
        'FP_Mapper_SNP': [f"Mapped: {mapper_snp_metrics['mapped_accuracy']:.6f}"],
        'FN_Mapper_SNP': [''],
        'Precision_Mapper_SNP': [''],
        'Recall_Mapper_SNP': [''],
        'F1_Score_Mapper_SNP': [''],
        'TP_Bowtie2': [f"Overall: {bowtie2_metrics['overall_accuracy']:.6f}"],
        'FP_Bowtie2': [f"Mapped: {bowtie2_metrics['mapped_accuracy']:.6f}"],
        'FN_Bowtie2': [''],
        'Precision_Bowtie2': [''],
        'Recall_Bowtie2': [''],
        'F1_Score_Bowtie2': ['']
    })
    
    # Combine all data
    complete_data = pd.concat([merged, agg_row, acc_row], ignore_index=True)
    
    # Save to TSV
    complete_data.to_csv(f'{output_path}_three_way_comprehensive_metrics.tsv', sep='\t', index=False)
    
    print(f"Three-way comprehensive metrics saved: {output_path}_three_way_comprehensive_metrics.tsv")

def process_three_way_comparison(mapper_tsv, mapper_snp_tsv, fastq_file, sam_file, output_folder, output_prefix):
    """Main processing function for three-way comparison"""
    
    # Create output folder
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_prefix)
    
    # Parse input files
    df_mapper = pd.read_csv(mapper_tsv, sep='\t')
    df_mapper_snp = pd.read_csv(mapper_snp_tsv, sep='\t')
    origin_dict = parse_fastq(fastq_file)
    sam_dict = parse_sam(sam_file)
    
    # Find valid reads (common to all methods)
    mapper_read_names = set(df_mapper['Read_Name'])
    mapper_snp_read_names = set(df_mapper_snp['Read_Name'])
    fastq_read_names = set(origin_dict.keys())
    valid_reads = mapper_read_names.intersection(mapper_snp_read_names).intersection(fastq_read_names)
    
    print(f"Reads in Full Mapper: {len(mapper_read_names)}")
    print(f"Reads in SNP-only Mapper: {len(mapper_snp_read_names)}")
    print(f"Reads in FASTQ: {len(fastq_read_names)}")
    print(f"Valid reads for analysis: {len(valid_reads)}")
    
    # Process Full Mapper data
    mapper_results = []
    for _, row in df_mapper.iterrows():
        read_name = row['Read_Name']
        if read_name in valid_reads:
            mapped_gene = row['Uniquely_Mapped']
            origin_gene = origin_dict[read_name]
            mapper_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 'Mapped_Gene': mapped_gene})
    mapper_df = pd.DataFrame(mapper_results)
    
    # Process SNP-only Mapper data
    mapper_snp_results = []
    for _, row in df_mapper_snp.iterrows():
        read_name = row['Read_Name']
        if read_name in valid_reads:
            mapped_gene = row['Uniquely_Mapped']
            origin_gene = origin_dict[read_name]
            mapper_snp_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 'Mapped_Gene': mapped_gene})
    mapper_snp_df = pd.DataFrame(mapper_snp_results)
    
    # Process Bowtie2 data
    bowtie2_results = []
    for read_name in valid_reads:
        origin_gene = origin_dict[read_name]
        mapped_gene = sam_dict.get(read_name, 'NONE')
        bowtie2_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 'Mapped_Gene': mapped_gene})
    bowtie2_df = pd.DataFrame(bowtie2_results)
    
    # Create confusion matrices
    mapper_cm = pd.crosstab(mapper_df['Origin_Gene'], mapper_df['Mapped_Gene'])
    mapper_snp_cm = pd.crosstab(mapper_snp_df['Origin_Gene'], mapper_snp_df['Mapped_Gene'])
    bowtie2_cm = pd.crosstab(bowtie2_df['Origin_Gene'], bowtie2_df['Mapped_Gene'])
    
    mapper_cm.to_csv(f'{output_path}_mapper_full_confusion_matrix.tsv', sep='\t')
    mapper_snp_cm.to_csv(f'{output_path}_mapper_snp_confusion_matrix.tsv', sep='\t')
    bowtie2_cm.to_csv(f'{output_path}_bowtie2_confusion_matrix.tsv', sep='\t')
    
    # Calculate metrics
    mapper_per_gene, mapper_agg_prec, mapper_agg_rec, mapper_agg_f1, mapper_overall_acc, mapper_mapped_acc, mapper_additional = calculate_metrics(mapper_df)
    mapper_snp_per_gene, mapper_snp_agg_prec, mapper_snp_agg_rec, mapper_snp_agg_f1, mapper_snp_overall_acc, mapper_snp_mapped_acc, mapper_snp_additional = calculate_metrics(mapper_snp_df)
    bowtie2_per_gene, bowtie2_agg_prec, bowtie2_agg_rec, bowtie2_agg_f1, bowtie2_overall_acc, bowtie2_mapped_acc, bowtie2_additional = calculate_metrics(bowtie2_df)
    
    # Package metrics for plotting
    mapper_metrics = {
        'agg_precision': mapper_agg_prec,
        'agg_recall': mapper_agg_rec,
        'agg_f1': mapper_agg_f1,
        'overall_accuracy': mapper_overall_acc,
        'mapped_accuracy': mapper_mapped_acc,
        **mapper_additional
    }
    
    mapper_snp_metrics = {
        'agg_precision': mapper_snp_agg_prec,
        'agg_recall': mapper_snp_agg_rec,
        'agg_f1': mapper_snp_agg_f1,
        'overall_accuracy': mapper_snp_overall_acc,
        'mapped_accuracy': mapper_snp_mapped_acc,
        **mapper_snp_additional
    }
    
    bowtie2_metrics = {
        'agg_precision': bowtie2_agg_prec,
        'agg_recall': bowtie2_agg_rec,
        'agg_f1': bowtie2_agg_f1,
        'overall_accuracy': bowtie2_overall_acc,
        'mapped_accuracy': bowtie2_mapped_acc,
        **bowtie2_additional
    }
    
    # Create three-way comparison plot
    create_three_way_bar_plot(mapper_metrics, mapper_snp_metrics, bowtie2_metrics, output_path)
    
    # Save comprehensive TSV
    save_three_way_comprehensive_tsv(mapper_per_gene, mapper_snp_per_gene, bowtie2_per_gene, 
                                    mapper_metrics, mapper_snp_metrics, bowtie2_metrics, output_path)
    
    # Print summary
    print(f"\n=== THREE-WAY COMPARISON SUMMARY ===")
    print(f"Mapper (Full)    - Precision: {mapper_agg_prec:.6f}, Recall: {mapper_agg_rec:.6f}, F1: {mapper_agg_f1:.6f}, Accuracy: {mapper_mapped_acc:.6f}")
    print(f"Mapper (SNP-only) - Precision: {mapper_snp_agg_prec:.6f}, Recall: {mapper_snp_agg_rec:.6f}, F1: {mapper_snp_agg_f1:.6f}, Accuracy: {mapper_snp_mapped_acc:.6f}")
    print(f"Bowtie2          - Precision: {bowtie2_agg_prec:.6f}, Recall: {bowtie2_agg_rec:.6f}, F1: {bowtie2_agg_f1:.6f}, Accuracy: {bowtie2_mapped_acc:.6f}")

# Variant calling comparison functions
def soft_normalize_variant(pos, ref, alt):
    """Normalize variant representation for comparison"""
    pos = int(pos)
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt

    prefix_len = 0
    max_prefix = min(len(ref), len(alt))
    while prefix_len < max_prefix and ref[prefix_len] == alt[prefix_len]:
        prefix_len += 1

    if prefix_len > 0:
        pos = pos + prefix_len - 1
        ref = ref[prefix_len-1:]
        alt = alt[prefix_len-1:]
    
    suffix_len = 0
    max_suffix = min(len(ref), len(alt)) - 1
    while suffix_len < max_suffix and ref[-(suffix_len+1)] == alt[-(suffix_len+1)]:
        suffix_len += 1
    if suffix_len > 0:
        ref = ref[:-suffix_len] if suffix_len < len(ref) else ref
        alt = alt[:-suffix_len] if suffix_len < len(alt) else alt

    return pos, ref, alt

def parse_vcf_for_variants(filename):
    """Parse VCF file and return variant dictionary"""
    records_by_gene = {}
    total_variants_by_gene = {}
    headers = None
    total_lines = 0
    total_variants = 0
    
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                headers = line.split("\t")
                chrom_index = headers.index("#CHROM")
                pos_index = headers.index("POS")
                ref_index = headers.index("REF")
                alt_index = headers.index("ALT")
                continue
            if headers is None:
                continue

            total_lines += 1
            parts = line.split("\t")
            gene = parts[chrom_index]
            pos = parts[pos_index]
            ref = parts[ref_index]
            alt = parts[alt_index]
            
            if gene not in total_variants_by_gene:
                total_variants_by_gene[gene] = 0
            total_variants_by_gene[gene] += 1
            total_variants += 1
            
            norm_pos, norm_ref, norm_alt = soft_normalize_variant(pos, ref, alt)
            key = (int(norm_pos), norm_ref, norm_alt)
            
            if gene not in records_by_gene:
                records_by_gene[gene] = {}
            
            records_by_gene[gene][key] = True
    
    print(f"File {filename}: {total_lines} total lines, {total_variants} variants processed")
    for gene, count in total_variants_by_gene.items():
        print(f"  Gene {gene}: {count} variants")
        
    return records_by_gene, total_variants_by_gene

def merge_method_vcfs(file_list):
    """Merge multiple VCF files into single dictionary"""
    merged = {}
    total_by_gene = {}
    
    for filename in file_list:
        gene_dict, gene_counts = parse_vcf_for_variants(filename)
        for gene, variants in gene_dict.items():
            if gene not in merged:
                merged[gene] = {}
            merged[gene].update(variants)
            if gene not in total_by_gene:
                total_by_gene[gene] = 0
            total_by_gene[gene] += gene_counts.get(gene, 0)
            
    return merged, total_by_gene

def compute_variant_confusion_matrix(gt_variants, method_variants):
    """Compute confusion matrix for variant calling"""
    all_variants = set(gt_variants.keys()) | set(method_variants.keys())
    
    y_true = [1 if variant in gt_variants else 0 for variant in all_variants]
    y_pred = [1 if variant in method_variants else 0 for variant in all_variants]
    
    precision, recall, f1, _ = precision_recall_fscore_support(y_true, y_pred, average='binary', zero_division=0)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    
    return tp, fp, fn, precision, recall, f1

def aggregate_variant_metrics(gt_dict, method_dict):
    """Aggregate variant metrics across all genes"""
    all_true_labels = []
    all_pred_labels = []
    
    for gene, gt_variants in gt_dict.items():
        method_variants = method_dict.get(gene, {})
        all_variants = set(gt_variants.keys()) | set(method_variants.keys())
        
        y_true = [1 if variant in gt_variants else 0 for variant in all_variants]
        y_pred = [1 if variant in method_variants else 0 for variant in all_variants]
        
        all_true_labels.extend(y_true)
        all_pred_labels.extend(y_pred)
    
    precision, recall, f1, _ = precision_recall_fscore_support(
        all_true_labels, all_pred_labels, average='binary', zero_division=0)
    
    tn, fp, fn, tp = confusion_matrix(all_true_labels, all_pred_labels, labels=[0, 1]).ravel()
    
    return tp, fp, fn, precision, recall, f1

def create_three_way_variant_bar_plot(method1_metrics, method2_metrics, method3_metrics, output_path):
    """Create three-way variant calling comparison bar plot"""
    
    plt.figure(figsize=(12, 7), dpi=300)
    
    # Data for the 3 groups
    metrics = ['Precision', 'Recall', 'F1-Score']
    method1_values = [method1_metrics['precision'], method1_metrics['recall'], method1_metrics['f1']]
    method2_values = [method2_metrics['precision'], method2_metrics['recall'], method2_metrics['f1']]
    method3_values = [method3_metrics['precision'], method3_metrics['recall'], method3_metrics['f1']]
    
    # Create bar positions
    bar_width = 0.08
    bar_spacing = 0.1
    group_spacing = 0.2
    x_adjusted = np.linspace(0, len(metrics) * (bar_width + group_spacing), len(metrics))
    
    # Create bars
    plt.bar(x_adjusted - bar_spacing, method1_values, width=bar_width,
            label='Mapper (Full)', color=COLORS[0])
    plt.bar(x_adjusted, method2_values, width=bar_width,
            label='Mapper (SNP-only)', color=COLORS[2])
    plt.bar(x_adjusted + bar_spacing, method3_values, width=bar_width,
            label='Bowtie2', color=COLORS[1])
    
    # Formatting
    plt.xticks(x_adjusted, metrics, fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value labels on bars
    for i, (v1, v2, v3) in enumerate(zip(method1_values, method2_values, method3_values)):
        plt.text(x_adjusted[i] - bar_spacing, v1 + 0.02, f"{v1:.4f}", ha="center", fontsize=12)
        plt.text(x_adjusted[i], v2 + 0.02, f"{v2:.4f}", ha="center", fontsize=12)
        plt.text(x_adjusted[i] + bar_spacing, v3 + 0.02, f"{v3:.4f}", ha="center", fontsize=12)
    
    plt.legend(fontsize=14, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=3)
    plt.ylim(0, 1.1)
    plt.margins(y=0.1, x=0.02)
    plt.tight_layout()
    
    # Save plot
    plt.savefig(f'{output_path}_three_way_variant_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Three-way variant comparison plot saved: {output_path}_three_way_variant_comparison.png")

def process_three_way_variant_comparison(gt_vcf, method1_vcfs, method2_vcfs, method3_vcf, output_folder, output_prefix):
    """Process three-way variant calling comparison"""
    
    # Create output folder
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_prefix)
    
    # Parse ground truth
    print("Parsing ground truth VCF...")
    gt_dict, gt_total_by_gene = parse_vcf_for_variants(gt_vcf)
    
    # Parse method files
    print("Parsing Method 1 VCFs (Mapper Full)...")
    method1_dict, m1_total_by_gene = merge_method_vcfs(method1_vcfs)
    
    print("Parsing Method 2 VCFs (Mapper SNP-only)...")
    method2_dict, m2_total_by_gene = merge_method_vcfs(method2_vcfs)
    
    print("Parsing Method 3 VCF (Bowtie2)...")
    method3_dict, m3_total_by_gene = parse_vcf_for_variants(method3_vcf)
    
    # Calculate aggregated metrics
    tp1_agg, fp1_agg, fn1_agg, prec1_agg, rec1_agg, f1_1_agg = aggregate_variant_metrics(gt_dict, method1_dict)
    tp2_agg, fp2_agg, fn2_agg, prec2_agg, rec2_agg, f1_2_agg = aggregate_variant_metrics(gt_dict, method2_dict)
    tp3_agg, fp3_agg, fn3_agg, prec3_agg, rec3_agg, f1_3_agg = aggregate_variant_metrics(gt_dict, method3_dict)
    
    method1_agg = {'tp': tp1_agg, 'fp': fp1_agg, 'fn': fn1_agg, 'precision': prec1_agg, 'recall': rec1_agg, 'f1': f1_1_agg}
    method2_agg = {'tp': tp2_agg, 'fp': fp2_agg, 'fn': fn2_agg, 'precision': prec2_agg, 'recall': rec2_agg, 'f1': f1_2_agg}
    method3_agg = {'tp': tp3_agg, 'fp': fp3_agg, 'fn': fn3_agg, 'precision': prec3_agg, 'recall': rec3_agg, 'f1': f1_3_agg}
    
    # Create plots
    create_three_way_variant_bar_plot(method1_agg, method2_agg, method3_agg, output_path)
    
    # Print summary
    print(f"\n=== THREE-WAY VARIANT CALLING COMPARISON SUMMARY ===")
    print(f"Mapper (Full)    - Precision: {prec1_agg:.6f}, Recall: {rec1_agg:.6f}, F1: {f1_1_agg:.6f}")
    print(f"Mapper (SNP-only) - Precision: {prec2_agg:.6f}, Recall: {rec2_agg:.6f}, F1: {f1_2_agg:.6f}")
    print(f"Bowtie2          - Precision: {prec3_agg:.6f}, Recall: {rec3_agg:.6f}, F1: {f1_3_agg:.6f}")

def main():
    parser = argparse.ArgumentParser(description='Three-way mapper comparison: Full Mapper vs SNP-only Mapper vs Bowtie2')
    parser.add_argument('-t1', '--mapper-tsv', required=True, help='TSV file with Full Mapper results')
    parser.add_argument('-t2', '--mapper-snp-tsv', required=True, help='TSV file with SNP-only Mapper results')
    parser.add_argument('-f', '--fastq', required=True, help='FASTQ file with read origins')
    parser.add_argument('-s', '--sam', required=True, help='SAM file with Bowtie2 results')
    parser.add_argument('-o', '--output-dir', default='three_way_comparison', help='Output directory')
    parser.add_argument('-p', '--prefix', default='three_way_mapper_comparison', help='Output file prefix')
    
    # Variant calling comparison arguments
    parser.add_argument('--variant-comparison', action='store_true', help='Enable variant calling comparison')
    parser.add_argument('--gt-vcf', help='Ground truth VCF file')
    parser.add_argument('--method1-vcfs', nargs='+', help='Method 1 VCF files (Mapper Full)')
    parser.add_argument('--method2-vcfs', nargs='+', help='Method 2 VCF files (Mapper SNP-only)')
    parser.add_argument('--method3-vcf', help='Method 3 VCF file (Bowtie2)')
    
    args = parser.parse_args()
    
    # Run read mapping comparison
    process_three_way_comparison(args.mapper_tsv, args.mapper_snp_tsv, args.fastq, args.sam, args.output_dir, args.prefix)
    
    # Run variant calling comparison if requested
    if args.variant_comparison:
        if not all([args.gt_vcf, args.method1_vcfs, args.method2_vcfs, args.method3_vcf]):
            print("ERROR: For variant comparison, you must provide --gt-vcf, --method1-vcfs, --method2-vcfs, and --method3-vcf")
            return
        
        print("\n" + "="*50)
        print("STARTING VARIANT CALLING COMPARISON")
        print("="*50)
        
        process_three_way_variant_comparison(args.gt_vcf, args.method1_vcfs, args.method2_vcfs, 
                                           args.method3_vcf, args.output_dir, args.prefix)

if __name__ == "__main__":
    main()