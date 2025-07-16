#!/usr/bin/env python3
import argparse
import csv
import os
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix

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

def write_tsv(filename, header, rows):
    with open(filename, "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(
        description="Compare two variant-calling methods against ground truth VCF - confusion matrices only"
    )
    parser.add_argument("ground_truth", help="Path to the ground truth VCF file")
    parser.add_argument("method1_files", nargs="+", help="List of VCF files for Method 1 (one per gene)")
    parser.add_argument("--method2", required=True, help="Path to the Method 2 VCF file")
    parser.add_argument("--tsv", required=True, help="TSV file with read mapping results (Mapper)")
    parser.add_argument("--sam", required=True, help="SAM file with read mapping results (Bowtie2)")
    parser.add_argument("--fastq", required=True, help="FASTQ file with read origins")
    parser.add_argument("--output-dir", default="output", help="Output directory for confusion matrices")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Parsing ground truth file...")
    gt_dict, gt_total_by_gene = parse_vcf(args.ground_truth)
    
    print("Parsing method1 files...")
    method1_dict, m1_total_by_gene = merge_method1_vcfs(args.method1_files)
    
    print("Parsing method2 file...")
    method2_dict, m2_total_by_gene = parse_vcf(args.method2)

    all_genes = sorted(set(list(gt_dict.keys()) + list(method1_dict.keys()) + list(method2_dict.keys())))
    confusion_matrix_rows = []
    
    for gene in all_genes:
        gt_variants = gt_dict.get(gene, {})
        m1_variants = method1_dict.get(gene, {})
        m2_variants = method2_dict.get(gene, {})

        TP1, FP1, FN1, precision1, recall1, f1_1 = compute_confusion_matrix_sklearn(gt_variants, m1_variants)
        TP2, FP2, FN2, precision2, recall2, f1_2 = compute_confusion_matrix_sklearn(gt_variants, m2_variants)
        
        confusion_matrix_rows.append([
            gene,
            TP1, FP1, FN1,
            TP2, FP2, FN2
        ])

    confusion_matrix_tsv = os.path.join(args.output_dir, "variant_calling_confusion_matrix.tsv")
    header = ["Gene", "Method1_TP", "Method1_FP", "Method1_FN", "Method2_TP", "Method2_FP", "Method2_FN"]
    
    write_tsv(confusion_matrix_tsv, header, confusion_matrix_rows)
    print(f"Variant calling confusion matrix saved: {confusion_matrix_tsv}")

if __name__ == "__main__":
    main()