#!/usr/bin/env python3

import sys
import csv
import os
from collections import defaultdict

def parse_vcf(vcf_file):
    """Parse VCF file and extract CHROM, POS, REF, ALT into a dictionary by CHROM."""
    variants_by_chrom = defaultdict(set)
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt = fields[:5]
            variant = (int(pos), ref, alt)  # Store as (position, ref, alt)
            variants_by_chrom[chrom].add(variant)
    return variants_by_chrom

def compare_vcf(ground_truth_variants, called_variants):
    """Compare ground truth variants with called variants for each chromosome."""
    results = {}

    all_chroms = set(ground_truth_variants.keys()) | set(called_variants.keys())

    for chrom in all_chroms:
        gt_variants = ground_truth_variants.get(chrom, set())
        call_variants = called_variants.get(chrom, set())

        true_positives = gt_variants & call_variants
        false_negatives = gt_variants - call_variants
        false_positives = call_variants - gt_variants

        results[chrom] = (true_positives, false_negatives, false_positives)
    
    return results

def calculate_metrics(tp, fn, fp):
    """Calculate precision, recall, and F1-score."""
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1

def write_tsv(output_file, all_results):
    """Write aggregated results with precision, recall, and F1-score to a TSV file."""
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['Chromosome', 'True_Positives', 'False_Negatives', 'False_Positives', 'Precision', 'Recall', 'F1_Score'])

        aggregated_results = defaultdict(lambda: {'tp': set(), 'fn': set(), 'fp': set()})

        # Aggregate results across all VCFs
        for results in all_results.values():
            for chrom, (tp, fn, fp) in results.items():
                aggregated_results[chrom]['tp'].update(tp)
                aggregated_results[chrom]['fn'].update(fn)
                aggregated_results[chrom]['fp'].update(fp)

        total_tp, total_fn, total_fp = 0, 0, 0

        for chrom in sorted(aggregated_results.keys()):
            tp = len(aggregated_results[chrom]['tp'])
            fn = len(aggregated_results[chrom]['fn'])
            fp = len(aggregated_results[chrom]['fp'])

            precision, recall, f1 = calculate_metrics(tp, fn, fp)

            writer.writerow([chrom, tp, fn, fp, f"{precision:.3f}", f"{recall:.3f}", f"{f1:.3f}"])

            total_tp += tp
            total_fn += fn
            total_fp += fp

        # Total metrics
        total_precision, total_recall, total_f1 = calculate_metrics(total_tp, total_fn, total_fp)

        writer.writerow(['Total', total_tp, total_fn, total_fp, f"{total_precision:.3f}", f"{total_recall:.3f}", f"{total_f1:.3f}"])

def main():
    if len(sys.argv) != 4:
        print("Usage: python compare_vcf_folder.py <ground_truth.vcf> <vcf_folder> <output.tsv>")
        sys.exit(1)

    ground_truth_file = sys.argv[1]
    vcf_folder = sys.argv[2]
    output_file = sys.argv[3]

    if not os.path.isfile(ground_truth_file):
        print(f"Error: Ground truth file '{ground_truth_file}' not found.")
        sys.exit(1)

    if not os.path.isdir(vcf_folder):
        print(f"Error: Folder '{vcf_folder}' not found.")
        sys.exit(1)

    ground_truth_variants = parse_vcf(ground_truth_file)
    print(f"Loaded ground truth variants: {sum(len(v) for v in ground_truth_variants.values())}")

    vcf_files = [f for f in os.listdir(vcf_folder) if f.endswith('.vcf')]
    if len(vcf_files) == 0:
        print(f"No VCF files found in folder '{vcf_folder}'.")
        sys.exit(1)

    all_results = {}

    for vcf_file in vcf_files:
        file_path = os.path.join(vcf_folder, vcf_file)
        called_variants = parse_vcf(file_path)
        print(f"Processing: {vcf_file} ({sum(len(v) for v in called_variants.values())} variants)")

        results = compare_vcf(ground_truth_variants, called_variants)
        all_results[vcf_file] = results

    write_tsv(output_file, all_results)
    print(f"Results written to '{output_file}'")

if __name__ == "__main__":
    main()
