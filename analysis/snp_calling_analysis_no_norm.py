#!/usr/bin/env python3
import argparse
import csv
import os
import pandas as pd
from sklearn.metrics import confusion_matrix


def parse_vcf(filename):
    """Parse VCF without any normalization or ALT splitting.
    Keys are taken exactly as (#CHROM, POS, REF, ALT) as strings (POS cast to int).
    Returns: (records_by_gene: dict[gene]->set[(pos, ref, alt)], counts_by_gene: dict[gene]->int)
    """
    records_by_gene = {}
    total_variants_by_gene = {}
    headers = None
    total_lines = 0
    total_variants = 0

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("##"):
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
            pos = int(parts[pos_index])
            ref = parts[ref_index]
            alt = parts[alt_index]  # keep ALT exactly as-is (including multi-allelic string if present)

            if gene not in total_variants_by_gene:
                total_variants_by_gene[gene] = 0
            total_variants_by_gene[gene] += 1
            total_variants += 1

            key = (pos, ref, alt)
            if gene not in records_by_gene:
                records_by_gene[gene] = set()
            records_by_gene[gene].add(key)

    print(f"File {filename}: {total_lines} total lines, {total_variants} variants processed (no normalization)")
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
                merged[gene] = set()
            merged[gene].update(variants)
            if gene not in total_by_gene:
                total_by_gene[gene] = 0
            total_by_gene[gene] += gene_counts.get(gene, 0)

    return merged, total_by_gene


def counts_from_sets(gt_set, method_set):
    """Given ground truth and method variant sets, compute TP/FP/FN/TN via confusion matrix.
    Universe is union of variants. True variant if in gt_set; predicted if in method_set.
    """
    all_variants = gt_set | method_set
    if not all_variants:
        return 0, 0, 0, 0

    y_true = [1 if v in gt_set else 0 for v in all_variants]
    y_pred = [1 if v in method_set else 0 for v in all_variants]

    # confusion_matrix with labels [0,1] returns [[TN, FP], [FN, TP]]
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return tp, fp, fn, tn


def write_tsv(filename, header, rows):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Compare two variant-calling methods against ground truth VCF without normalization (exact tuple match)."
    )
    parser.add_argument("ground_truth", help="Path to the ground truth VCF file")
    parser.add_argument("method1_files", nargs="+", help="List of VCF files for Method 1 (one per gene)")
    parser.add_argument("--method2", required=True, help="Path to the Method 2 VCF file")
    # Unused but accepted for CLI compatibility
    parser.add_argument("--tsv", required=False, help="TSV file (ignored)")
    parser.add_argument("--sam", required=False, help="SAM file (ignored)")
    parser.add_argument("--fastq", required=False, help="FASTQ file (ignored)")
    parser.add_argument("--output-dir", default="output", help="Output directory for confusion matrices")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Parsing ground truth file (no normalization)...")
    gt_dict, _ = parse_vcf(args.ground_truth)

    print("Parsing method1 files (no normalization)...")
    method1_dict, _ = merge_method1_vcfs(args.method1_files)

    print("Parsing method2 file (no normalization)...")
    method2_dict, _ = parse_vcf(args.method2)

    all_genes = sorted(set(list(gt_dict.keys()) + list(method1_dict.keys()) + list(method2_dict.keys())))
    confusion_matrix_rows = []

    for gene in all_genes:
        gt_variants = gt_dict.get(gene, set())
        m1_variants = method1_dict.get(gene, set())
        m2_variants = method2_dict.get(gene, set())

        TP1, FP1, FN1, TN1 = counts_from_sets(gt_variants, m1_variants)
        TP2, FP2, FN2, TN2 = counts_from_sets(gt_variants, m2_variants)

        confusion_matrix_rows.append([
            gene,
            TP1, FP1, FN1, TN1,
            TP2, FP2, FN2, TN2
        ])

    confusion_matrix_tsv = os.path.join(args.output_dir, "variant_calling_confusion_matrix.tsv")
    header = [
        "Gene",
        "Method1_TP", "Method1_FP", "Method1_FN", "Method1_TN",
        "Method2_TP", "Method2_FP", "Method2_FN", "Method2_TN"
    ]

    write_tsv(confusion_matrix_tsv, header, confusion_matrix_rows)
    print(f"Variant calling confusion matrix saved (no normalization): {confusion_matrix_tsv}")


if __name__ == "__main__":
    main()

