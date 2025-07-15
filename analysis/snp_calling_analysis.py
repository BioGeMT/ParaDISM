#!/usr/bin/env python3
import argparse
import csv
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix

COLORS = ['#FF8C00', '#008080']  

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 24,
    'axes.titlesize': 28,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 20
})

def soft_normalize_variant(pos, ref, alt):
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

def count_mapped_reads_by_gene(tsv_file, sam_file, fastq_file):
    read_to_origin = {}
    with open(fastq_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                parts = line.strip().split()
                read_name = parts[0][1:]
                gene_info = parts[1].split(';')[0]
                read_to_origin[read_name] = gene_info
    
    df_tsv = pd.read_csv(tsv_file, sep='\t')
    mapper_counts = {}
    
    for _, row in df_tsv.iterrows():
        read_name = row['Read_Name']
        mapped_gene = row['Uniquely_Mapped']
        
        if read_name in read_to_origin:
            origin_gene = read_to_origin[read_name]
            
            if origin_gene not in mapper_counts:
                mapper_counts[origin_gene] = 0
            
            if mapped_gene != 'NONE' and mapped_gene == origin_gene:
                mapper_counts[origin_gene] += 1
    
    bowtie2_counts = {}
    with open(sam_file, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                parts = line.strip().split('\t')
                read_name = parts[0]
                mapped_gene = parts[2]
                
                if read_name in read_to_origin:
                    origin_gene = read_to_origin[read_name]
                    
                    if origin_gene not in bowtie2_counts:
                        bowtie2_counts[origin_gene] = 0
                    
                    if mapped_gene != '*' and mapped_gene == origin_gene:
                        bowtie2_counts[origin_gene] += 1
    
    return mapper_counts, bowtie2_counts

def parse_vcf(filename):
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

def merge_method1_vcfs(file_list):
    merged = {}
    total_by_gene = {}
    
    for filename in file_list:
        gene_dict, gene_counts = parse_vcf(filename)
        for gene, variants in gene_dict.items():
            if gene not in merged:
                merged[gene] = {}
            merged[gene].update(variants)
            if gene not in total_by_gene:
                total_by_gene[gene] = 0
            total_by_gene[gene] += gene_counts.get(gene, 0)
            
    return merged, total_by_gene

def compute_confusion_matrix_sklearn(gt_variants, method_variants):
    all_variants = set(gt_variants.keys()) | set(method_variants.keys())
    
    y_true = [1 if variant in gt_variants else 0 for variant in all_variants]
    y_pred = [1 if variant in method_variants else 0 for variant in all_variants]
    
    precision, recall, f1, _ = precision_recall_fscore_support(y_true, y_pred, average='binary', zero_division=0)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    
    return tp, fp, fn, precision, recall, f1

def aggregate_metrics_sklearn(gt_dict, method_dict):
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

def write_tsv(filename, header, rows):
    with open(filename, "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

def create_bar_plot(metrics, method1_vals, method2_vals, title, output_file):
    plt.rcParams.update({'font.size': 18, 'font.weight': 'bold', 'axes.titlesize': 22})
    plt.figure(figsize=(10, 6), dpi=300)
    
    bar_width = 0.1
    bar_spacing = 0.12
    group_spacing = 0.15
    x_adjusted = np.linspace(0, len(metrics) * (bar_width + group_spacing), len(metrics))
    
    plt.bar(x_adjusted - bar_spacing/2, method1_vals, width=bar_width, label='Mapper', color=COLORS[0])
    plt.bar(x_adjusted + bar_spacing/2, method2_vals, width=bar_width, label='Bowtie2', color=COLORS[1])
    
    plt.xticks(x_adjusted, metrics, fontsize=15)
    plt.yticks(fontsize=20)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    for i, (v1, v2) in enumerate(zip(method1_vals, method2_vals)):
        plt.text(x_adjusted[i] - bar_spacing/2, v1 + 0.03, f"{v1:.4f}", ha="center", fontsize=14)
        plt.text(x_adjusted[i] + bar_spacing/2, v2 + 0.03, f"{v2:.4f}", ha="center", fontsize=14)
    
    plt.legend(fontsize=18, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
    plt.ylim(0, 1.1)
    plt.margins(y=0.1, x=0.01)
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Plot saved to: {output_file}")

def create_precision_vs_mapped_scatter(per_gene_metrics, read_counts_mapper, read_counts_bowtie2, output_file):
    plt.figure(figsize=(10, 8))
    
    mapper_color = COLORS[0]
    bowtie2_color = COLORS[1]
    
    genes = []
    mapper_precision = []
    mapper_reads_mapped = []
    bowtie2_precision = []
    bowtie2_reads_mapped = []
    
    for gene_data in per_gene_metrics:
        gene = gene_data['gene']
        genes.append(gene)
        mapper_precision.append(gene_data['mapper_precision'])
        bowtie2_precision.append(gene_data['bowtie2_precision'])
        
        mapper_reads_mapped.append(read_counts_mapper.get(gene, 0))
        bowtie2_reads_mapped.append(read_counts_bowtie2.get(gene, 0) // 2)
    
    plt.scatter(mapper_reads_mapped, mapper_precision, 
               c=mapper_color, s=100, alpha=0.8, label='Mapper', marker='o')
    plt.scatter(bowtie2_reads_mapped, bowtie2_precision, 
               c=bowtie2_color, s=100, alpha=0.8, label='Bowtie2', marker='s')
    
    for i, gene in enumerate(genes):
        if mapper_reads_mapped[i] > 0:
            plt.annotate(gene, (mapper_reads_mapped[i], mapper_precision[i]), 
                        xytext=(5, 5), textcoords='offset points', fontsize=12, 
                        color=mapper_color, fontweight='bold')
        if bowtie2_reads_mapped[i] > 0:
            plt.annotate(gene, (bowtie2_reads_mapped[i], bowtie2_precision[i]), 
                        xytext=(5, -15), textcoords='offset points', fontsize=12, 
                        color=bowtie2_color, fontweight='bold')
    
    plt.ylim(-0.05, 1.05)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=16, loc='lower right')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Precision vs mapped reads scatter plot saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Compare two variant-calling methods (Method 1: per gene VCFs; Method 2: single VCF) against a common ground truth VCF. Evaluates all variant detection performance."
    )
    parser.add_argument("ground_truth", help="Path to the ground truth VCF file")
    parser.add_argument("method1_files", nargs="+", help="List of VCF files for Method 1 (one per gene)")
    parser.add_argument("--method2", required=True, help="Path to the Method 2 VCF file")
    parser.add_argument("--tsv", required=True, help="TSV file with read mapping results (Mapper)")
    parser.add_argument("--sam", required=True, help="SAM file with read mapping results (Bowtie2)")
    parser.add_argument("--fastq", required=True, help="FASTQ file with read origins")
    parser.add_argument("--output-dir", default="output", help="Output directory for TSVs and plots")
    args = parser.parse_args()

    # Create the output directory if it doesn't exist.
    os.makedirs(args.output_dir, exist_ok=True)

    # Count mapped reads by gene for scatter plot
    print("Counting mapped reads by gene...")
    read_counts_mapper, read_counts_bowtie2 = count_mapped_reads_by_gene(args.tsv, args.sam, args.fastq)

    # Parse VCF files.
    print("Parsing ground truth file...")
    gt_dict, gt_total_by_gene = parse_vcf(args.ground_truth)
    
    print("Parsing method1 files...")
    method1_dict, m1_total_by_gene = merge_method1_vcfs(args.method1_files)
    
    print("Parsing method2 file...")
    method2_dict, m2_total_by_gene = parse_vcf(args.method2)

    # Calculate total variants
    gt_total_all = sum(gt_total_by_gene.values())
    m1_total_all = sum(m1_total_by_gene.values())
    m2_total_all = sum(m2_total_by_gene.values())
    
    # Process each gene and collect data for comprehensive TSV and scatter plot
    all_genes = sorted(set(list(gt_dict.keys()) + list(method1_dict.keys()) + list(method2_dict.keys())))
    per_gene_metrics = []
    comprehensive_rows = []
    
    for gene in all_genes:
        gt_variants = gt_dict.get(gene, {})
        m1_variants = method1_dict.get(gene, {})
        m2_variants = method2_dict.get(gene, {})

        # Get variant counts for this gene
        gt_total = gt_total_by_gene.get(gene, 0)
        m1_total = m1_total_by_gene.get(gene, 0)
        m2_total = m2_total_by_gene.get(gene, 0)

        # Compute confusion matrices using sklearn
        TP1, FP1, FN1, precision1, recall1, f1_1 = compute_confusion_matrix_sklearn(gt_variants, m1_variants)
        TP2, FP2, FN2, precision2, recall2, f1_2 = compute_confusion_matrix_sklearn(gt_variants, m2_variants)
        
        # Calculate accuracy (TP / total variants for this gene)
        total1 = len(set(gt_variants.keys()) | set(m1_variants.keys()))
        total2 = len(set(gt_variants.keys()) | set(m2_variants.keys()))
        accuracy1 = TP1 / total1 if total1 > 0 else 0
        accuracy2 = TP2 / total2 if total2 > 0 else 0

        # Store data for scatter plot (keep variant precision, but use read counts)
        per_gene_metrics.append({
            'gene': gene,
            'mapper_precision': precision1,
            'bowtie2_precision': precision2
        })
        
        # Store data for comprehensive TSV
        comprehensive_rows.append([
            gene, gt_total, m1_total, m2_total,
            TP1, FP1, FN1, f"{precision1:.6f}", f"{recall1:.6f}", f"{f1_1:.6f}", f"{accuracy1:.6f}",
            TP2, FP2, FN2, f"{precision2:.6f}", f"{recall2:.6f}", f"{f1_2:.6f}", f"{accuracy2:.6f}"
        ])

    # Compute aggregated metrics using sklearn
    m1_TP, m1_FP, m1_FN, m1_precision, m1_recall, m1_f1 = aggregate_metrics_sklearn(gt_dict, method1_dict)
    m2_TP, m2_FP, m2_FN, m2_precision, m2_recall, m2_f1 = aggregate_metrics_sklearn(gt_dict, method2_dict)
    
    # Calculate accuracy for aggregated metrics
    m1_total = m1_TP + m1_FP + m1_FN + (len(set().union(*[set(variants.keys()) for variants in gt_dict.values()])) - m1_TP - m1_FP - m1_FN)
    m2_total = m2_TP + m2_FP + m2_FN + (len(set().union(*[set(variants.keys()) for variants in gt_dict.values()])) - m2_TP - m2_FP - m2_FN)
    m1_accuracy = m1_TP / m1_total if m1_total > 0 else 0
    m2_accuracy = m2_TP / m2_total if m2_total > 0 else 0

    # Write comprehensive TSV with per-gene and aggregated data
    comprehensive_tsv = os.path.join(args.output_dir, "comprehensive_variant_metrics.tsv")
    header = ["Gene", "GT_Variants", "Method1_Variants", "Method2_Variants",
              "Method1_TP", "Method1_FP", "Method1_FN", "Method1_Precision", "Method1_Recall", "Method1_F1", "Method1_Accuracy",
              "Method2_TP", "Method2_FP", "Method2_FN", "Method2_Precision", "Method2_Recall", "Method2_F1", "Method2_Accuracy"]
    
    # Add aggregated row
    comprehensive_rows.append([
        "AGGREGATED", gt_total_all, m1_total_all, m2_total_all,
        m1_TP, m1_FP, m1_FN, f"{m1_precision:.6f}", f"{m1_recall:.6f}", f"{m1_f1:.6f}", f"{m1_accuracy:.6f}",
        m2_TP, m2_FP, m2_FN, f"{m2_precision:.6f}", f"{m2_recall:.6f}", f"{m2_f1:.6f}", f"{m2_accuracy:.6f}"
    ])
    
    write_tsv(comprehensive_tsv, header, comprehensive_rows)
    print(f"Comprehensive variant metrics written to: {comprehensive_tsv}")

    # Create aggregated bar plot
    agg_plot = os.path.join(args.output_dir, "aggregated_variant_metrics.png")
    metrics = ["Precision", "Recall", "F1-Score"]
    m1_vals = [m1_precision, m1_recall, m1_f1]
    m2_vals = [m2_precision, m2_recall, m2_f1]
    create_bar_plot(metrics, m1_vals, m2_vals,
                    title="",
                    output_file=agg_plot)

    # Create precision vs mapped reads scatter plot
    scatter_plot = os.path.join(args.output_dir, "precision_vs_mapped_reads_scatter.png")
    create_precision_vs_mapped_scatter(per_gene_metrics, read_counts_mapper, read_counts_bowtie2, scatter_plot)
    
    # Print summary
    print(f"\n=== AGGREGATED METRICS SUMMARY ===")
    print(f"Mapper - Precision: {m1_precision:.6f}, Recall: {m1_recall:.6f}, F1: {m1_f1:.6f}")
    print(f"Bowtie2 - Precision: {m2_precision:.6f}, Recall: {m2_recall:.6f}, F1: {m2_f1:.6f}")

if __name__ == "__main__":
    main()