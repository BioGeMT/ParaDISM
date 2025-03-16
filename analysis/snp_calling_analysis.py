#!/usr/bin/env python3
import argparse
import csv
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.font_manager import FontProperties

# Define color palette
COLORS = ['#FF8C00', '#008080']  # Dark Orange for Method1/Mapper, Teal for Method2/Bowtie2

# Set default font sizes
plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 24,
    'axes.titlesize': 28,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 20
})

def soft_normalize_variant(pos, ref, alt):
    """
    Soft normalize a variant by trimming common prefixes and suffixes (keeping one base)
    without iterative left alignment.
    """
    pos = int(pos)
    # For SNPs, no normalization is needed.
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt

    # Compute common prefix length.
    prefix_len = 0
    max_prefix = min(len(ref), len(alt))
    while prefix_len < max_prefix and ref[prefix_len] == alt[prefix_len]:
        prefix_len += 1

    # Retain one base from the common prefix.
    if prefix_len > 0:
        pos = pos + prefix_len - 1
        ref = ref[prefix_len-1:]
        alt = alt[prefix_len-1:]
    
    # Compute common suffix length (leave at least one base).
    suffix_len = 0
    max_suffix = min(len(ref), len(alt)) - 1
    while suffix_len < max_suffix and ref[-(suffix_len+1)] == alt[-(suffix_len+1)]:
        suffix_len += 1
    if suffix_len > 0:
        ref = ref[:-suffix_len] if suffix_len < len(ref) else ref
        alt = alt[:-suffix_len] if suffix_len < len(alt) else alt

    return pos, ref, alt

def parse_vcf(filename):
    """
    Parse a VCF file and return a dictionary grouping variants by gene.
    (Assumes the gene is given in the CHROM column.)
    """
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
                format_index = headers.index("FORMAT") if "FORMAT" in headers else -1
                continue
            if headers is None:
                continue

            total_lines += 1
            parts = line.split("\t")
            gene = parts[chrom_index]
            pos = parts[pos_index]
            ref = parts[ref_index]
            alt = parts[alt_index]
            
            # Check if we need to consider genotype
            is_called_variant = True
            if format_index != -1 and len(parts) > format_index + 1:
                format_fields = parts[format_index].split(":")
                sample_fields = parts[format_index + 1].split(":")
                if "GT" in format_fields:
                    gt_index = format_fields.index("GT")
                    genotype = sample_fields[gt_index]
                    is_called_variant = genotype != "0/0"
            
            # Only count actual variants (not 0/0) in the total
            if is_called_variant:
                if gene not in total_variants_by_gene:
                    total_variants_by_gene[gene] = 0
                total_variants_by_gene[gene] += 1
                total_variants += 1
            
            # Normalize the variant
            norm_pos, norm_ref, norm_alt = soft_normalize_variant(pos, ref, alt)
            key = (int(norm_pos), norm_ref, norm_alt)
            
            if gene not in records_by_gene:
                records_by_gene[gene] = {}
            
            # Only store actual variants (not 0/0)
            if is_called_variant:
                records_by_gene[gene][key] = True
    
    print(f"File {filename}: {total_lines} total lines, {total_variants} actual variants processed")
    for gene, count in total_variants_by_gene.items():
        print(f"  Gene {gene}: {count} variants")
        
    return records_by_gene, total_variants_by_gene

def merge_method1_vcfs(file_list):
    """
    Given a list of VCF files (each for one gene) for Method 1,
    merge their parsed dictionaries into a single dictionary grouped by gene.
    """
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

def compute_confusion_matrix(gt_variants, method_variants):
    """
    Compute confusion matrix for all variant types.
    
    Returns: TP, FP, FN, total_variants.
    """
    # Get sets of all variant positions
    gt_keys = set(gt_variants.keys())
    method_keys = set(method_variants.keys())
    
    # True Positives: variants in both sets
    TP = len(gt_keys & method_keys)
    
    # False Negatives: variants in ground truth but not in method
    FN = len(gt_keys - method_keys)
    
    # False Positives: variants in method but not in ground truth
    FP = len(method_keys - gt_keys)
    
    # Total evaluated variants (for accuracy calculation)
    total_variants = len(gt_keys | method_keys)
    
    return TP, FP, FN, total_variants

def aggregate_metrics(gt_dict, method_dict):
    """
    Aggregate confusion matrix counts over all genes.
    """
    agg_TP = agg_FP = agg_FN = total_variants = 0
    for gene, gt_variants in gt_dict.items():
        method_variants = method_dict.get(gene, {})
        TP, FP, FN, total = compute_confusion_matrix(gt_variants, method_variants)
        agg_TP += TP
        agg_FP += FP
        agg_FN += FN
        total_variants += total
    return agg_TP, agg_FP, agg_FN, total_variants

def write_tsv(filename, header, rows):
    """
    Write a TSV file with the given header and rows.
    """
    with open(filename, "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

def create_bar_plot(metrics, method1_vals, method2_vals, title, output_file):
    """
    Create a grouped bar plot comparing method1 and method2 for the given metrics.
    """
    plt.rcParams.update({'font.size': 18, 'font.weight': 'bold', 'axes.titlesize': 22})
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Create bar positions
    bar_width = 0.1
    bar_spacing = 0.12
    group_spacing = 0.15
    x_adjusted = np.linspace(0, len(metrics) * (bar_width + group_spacing), len(metrics))
    
    # Create bars
    plt.bar(x_adjusted - bar_spacing/2, method1_vals, width=bar_width, label='Mapper', color=COLORS[0])
    plt.bar(x_adjusted + bar_spacing/2, method2_vals, width=bar_width, label='Bowtie2', color=COLORS[1])
    
    # Set up axes and labels
    plt.xticks(x_adjusted, metrics, fontsize=15)
    plt.yticks(fontsize=20)
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)
    
    # Add value annotations on top of bars
    for i, (v1, v2) in enumerate(zip(method1_vals, method2_vals)):
        plt.text(x_adjusted[i] - bar_spacing/2, v1 + 0.03, f"{v1:.2f}", ha="center", fontsize=14)
        plt.text(x_adjusted[i] + bar_spacing/2, v2 + 0.03, f"{v2:.2f}", ha="center", fontsize=14)
    
    # Add legend and title with consistent styling
    plt.legend(fontsize=18, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
    plt.title(title, fontsize=22, fontweight='bold')
    
    # Set y-axis range from 0 to 1.1 to give more space at the top
    plt.ylim(0, 1.1)
    
    # Adjust margins and layout
    plt.margins(y=0.1, x=0.01)
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Plot saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Compare two variant-calling methods (Method 1: per gene VCFs; Method 2: single VCF) against a common ground truth VCF. Evaluates all variant detection performance."
    )
    parser.add_argument("ground_truth", help="Path to the ground truth VCF file")
    parser.add_argument("method1_files", nargs="+", help="List of VCF files for Method 1 (one per gene)")
    parser.add_argument("--method2", required=True, help="Path to the Method 2 VCF file")
    parser.add_argument("--output-dir", default="output", help="Output directory for TSVs and plots")
    args = parser.parse_args()

    # Create the output directory if it doesn't exist.
    os.makedirs(args.output_dir, exist_ok=True)

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
    
    # Write total variant counts TSV
    variant_counts_tsv = os.path.join(args.output_dir, "total_variant_counts.tsv")
    variant_header = ["Gene", "Ground Truth Variants", "Mapper Variants", "Bowtie2 Variants"]
    variant_rows = []
    
    for gene in sorted(set(list(gt_total_by_gene.keys()) + list(m1_total_by_gene.keys()) + list(m2_total_by_gene.keys()))):
        gt_count = gt_total_by_gene.get(gene, 0)
        m1_count = m1_total_by_gene.get(gene, 0)
        m2_count = m2_total_by_gene.get(gene, 0)
        variant_rows.append([gene, gt_count, m1_count, m2_count])
    
    write_tsv(variant_counts_tsv, variant_header, variant_rows)
    print(f"Total variant counts written to: {variant_counts_tsv}")

    # Process each gene
    all_genes = sorted(set(list(gt_dict.keys()) + list(method1_dict.keys()) + list(method2_dict.keys())))
    
    # Process metrics for each gene
    for gene in all_genes:
        gt_variants = gt_dict.get(gene, {})
        m1_variants = method1_dict.get(gene, {})
        m2_variants = method2_dict.get(gene, {})

        # Get variant counts for this gene
        gt_total = gt_total_by_gene.get(gene, 0)
        m1_total = m1_total_by_gene.get(gene, 0)
        m2_total = m2_total_by_gene.get(gene, 0)

        # Compute confusion matrices for ALL variants.
        TP1, FP1, FN1, total1 = compute_confusion_matrix(gt_variants, m1_variants)
        TP2, FP2, FN2, total2 = compute_confusion_matrix(gt_variants, m2_variants)

        # Calculate ALL variants metrics.
        precision1 = TP1 / (TP1 + FP1) if (TP1 + FP1) > 0 else 0
        recall1 = TP1 / (TP1 + FN1) if (TP1 + FN1) > 0 else 0
        accuracy1 = TP1 / total1 if total1 > 0 else 0

        precision2 = TP2 / (TP2 + FP2) if (TP2 + FP2) > 0 else 0
        recall2 = TP2 / (TP2 + FN2) if (TP2 + FN2) > 0 else 0
        accuracy2 = TP2 / total2 if total2 > 0 else 0

        # Save per-gene results for ALL variant metrics.
        gene_tsv = os.path.join(args.output_dir, f"gene_{gene}_metrics.tsv")
        header = ["Method", "Total Variants", "GT Variants", "TP", "FP", "FN", "Precision", "Recall", "Accuracy"]
        rows = [
            ["Mapper", m1_total, gt_total, TP1, FP1, FN1, f"{precision1:.3f}", f"{recall1:.3f}", f"{accuracy1:.3f}"],
            ["Bowtie2", m2_total, gt_total, TP2, FP2, FN2, f"{precision2:.3f}", f"{recall2:.3f}", f"{accuracy2:.3f}"]
        ]
        write_tsv(gene_tsv, header, rows)
        print(f"Gene {gene}: variant metrics written to: {gene_tsv}")
        
        # Create per-gene ALL variants metrics plot.
        gene_plot = os.path.join(args.output_dir, f"gene_{gene}_metrics.png")
        metrics = ["Precision", "Recall", "Accuracy"]
        m1_vals = [precision1, recall1, accuracy1]
        m2_vals = [precision2, recall2, accuracy2]
        create_bar_plot(metrics, m1_vals, m2_vals,
                        title=f"SNP Calling Performance for {gene}",
                        output_file=gene_plot)

    # Compute aggregated metrics for ALL variants.
    m1_TP, m1_FP, m1_FN, m1_total = aggregate_metrics(gt_dict, method1_dict)
    m2_TP, m2_FP, m2_FN, m2_total = aggregate_metrics(gt_dict, method2_dict)

    m1_precision = m1_TP / (m1_TP + m1_FP) if (m1_TP + m1_FP) > 0 else 0
    m1_recall = m1_TP / (m1_TP + m1_FN) if (m1_TP + m1_FN) > 0 else 0
    m1_accuracy = m1_TP / m1_total if m1_total > 0 else 0

    m2_precision = m2_TP / (m2_TP + m2_FP) if (m2_TP + m2_FP) > 0 else 0
    m2_recall = m2_TP / (m2_TP + m2_FN) if (m2_TP + m2_FN) > 0 else 0
    m2_accuracy = m2_TP / m2_total if m2_total > 0 else 0

    # Write aggregated ALL variants metrics TSV.
    agg_tsv = os.path.join(args.output_dir, "aggregated_metrics.tsv")
    header = ["Method", "Total Variants", "GT Variants", "TP", "FP", "FN", "Precision", "Recall", "Accuracy"]
    rows = [
        ["Mapper", m1_total_all, gt_total_all, m1_TP, m1_FP, m1_FN, f"{m1_precision:.3f}", f"{m1_recall:.3f}", f"{m1_accuracy:.3f}"],
        ["Bowtie2", m2_total_all, gt_total_all, m2_TP, m2_FP, m2_FN, f"{m2_precision:.3f}", f"{m2_recall:.3f}", f"{m2_accuracy:.3f}"]
    ]
    write_tsv(agg_tsv, header, rows)
    print(f"Aggregated variant metrics written to: {agg_tsv}")

    # Create aggregated ALL variants metrics plot.
    agg_plot = os.path.join(args.output_dir, "aggregated_metrics.png")
    metrics = ["Precision", "Recall", "Accuracy"]
    m1_vals = [m1_precision, m1_recall, m1_accuracy]
    m2_vals = [m2_precision, m2_recall, m2_accuracy]
    create_bar_plot(metrics, m1_vals, m2_vals,
                    title="Aggregated SNP Calling Performance",
                    output_file=agg_plot)

if __name__ == "__main__":
    main()