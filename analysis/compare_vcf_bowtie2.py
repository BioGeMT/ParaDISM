#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict
import os

def parse_vcf(vcf_file):
    """Parse VCF file and extract CHROM, POS, REF, ALT into a dictionary by CHROM."""
    variants_by_chrom = defaultdict(set)
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos, _, ref, alt = fields[:5]
            variant = (int(pos), ref, alt)  # Store as tuple (pos, ref, alt)
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

def write_tsv(output_file, results):
    """Write results with precision, recall, and F1-score to a TSV file."""
    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['Chromosome', 'True_Positives', 'False_Negatives', 'False_Positives', 'Precision', 'Recall', 'F1_Score'])

        total_tp, total_fn, total_fp = 0, 0, 0

        for chrom, (tp, fn, fp) in sorted(results.items()):
            tp_count, fn_count, fp_count = len(tp), len(fn), len(fp)
            precision, recall, f1 = calculate_metrics(tp_count, fn_count, fp_count)

            writer.writerow([chrom, tp_count, fn_count, fp_count, f"{precision:.3f}", f"{recall:.3f}", f"{f1:.3f}"])

            total_tp += tp_count
            total_fn += fn_count
            total_fp += fp_count

        # Overall totals and metrics
        total_precision, total_recall, total_f1 = calculate_metrics(total_tp, total_fn, total_fp)
        writer.writerow(['Total', total_tp, total_fn, total_fp, f"{total_precision:.3f}", f"{total_recall:.3f}", f"{total_f1:.3f}"])

def main():
    if len(sys.argv) != 4:
        print("Usage: python script_compare.py <ground_truth.vcf> <called_vcf> <output.tsv>")
        sys.exit(1)

    ground_truth_file = sys.argv[1]
    called_vcf_file = sys.argv[2]
    output_file = sys.argv[3]

    if not os.path.isfile(ground_truth_file):
        print(f"Error: Ground truth file '{ground_truth_file}' not found.")
        sys.exit(1)
    if not os.path.isfile(called_vcf_file):
        print(f"Error: Called VCF file '{called_vcf_file}' not found.")
        sys.exit(1)

    # Parse VCFs
    ground_truth_variants = parse_vcf(ground_truth_file)
    called_variants = parse_vcf(called_vcf_file)

    print(f"Ground truth variants loaded: {sum(len(v) for v in ground_truth_variants.values())}")
    print("Ground truth chromosomes:", list(ground_truth_variants.keys()))
    print(f"Called variants loaded: {sum(len(v) for v in called_variants.values())}")
    print("Called chromosomes:", list(called_variants.keys()))

    # Compare
    results = compare_vcf(ground_truth_variants, called_variants)

    # Write to TSV
    write_tsv(output_file, results)
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()
