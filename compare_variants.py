#!/usr/bin/env python3
"""
Compare ParaDISM variants to GIAB ground truth and calculate Precision/Recall/Specificity/F1.

Calculates TP, FP, FN, Precision, Recall, Specificity, F1 for:
- Baseline (iteration_1) for all 3 aligners
- Final outputs (full ParaDISM) for all 3 aligners
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple

try:
    from Bio.Seq import Seq
except ImportError:
    print("Error: BioPython not found. Install with: conda install biopython", file=sys.stderr)
    sys.exit(1)


def load_coordinate_mapping(mapping_file: Path) -> Dict[str, Dict]:
    """Load coordinate mapping from TSV file."""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            gene = parts[0]
            chr16_start = int(parts[1])
            chr16_end = int(parts[2])
            strand = parts[4]
            mapping[gene] = {
                'chr16_start': chr16_start,
                'chr16_end': chr16_end,
                'strand': strand,
                'length': chr16_end - chr16_start + 1
            }
    return mapping


def parse_vcf_line(line: str) -> Dict:
    """Parse a VCF line into a dictionary."""
    parts = line.strip().split('\t')
    if len(parts) < 5:
        return None
    
    return {
        'chrom': parts[0],
        'pos': int(parts[1]),
        'id': parts[2],
        'ref': parts[3],
        'alt': parts[4],
        'qual': parts[5],
        'filter': parts[6],
        'info': parts[7] if len(parts) > 7 else ''
    }


def chr16_to_gene_coord(chr16_pos: int, gene: str, mapping: Dict[str, Dict]) -> int:
    """Convert chr16 coordinate to gene coordinate."""
    if gene not in mapping:
        return None
    
    gene_start = mapping[gene]['chr16_start']
    gene_end = mapping[gene]['chr16_end']
    
    if gene_start <= chr16_pos <= gene_end:
        return chr16_pos - gene_start + 1
    
    return None


def determine_gene_from_chr16(chr16_pos: int, mapping: Dict[str, Dict]) -> str:
    """Determine which gene a chr16 position belongs to."""
    for gene, info in mapping.items():
        if info['chr16_start'] <= chr16_pos <= info['chr16_end']:
            return gene
    return None


def load_giab_variants(giab_vcf: Path, mapping: Dict[str, Dict]) -> Set[Tuple]:
    """Load GIAB variants as (gene, gene_pos, ref, alt) tuples."""
    variants = set()
    
    open_func = gzip.open if giab_vcf.suffix == '.gz' else open
    mode = 'rt' if giab_vcf.suffix == '.gz' else 'r'
    
    with open_func(giab_vcf, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            variant = parse_vcf_line(line)
            if not variant:
                continue
            
            chr16_pos = variant['pos']
            gene = determine_gene_from_chr16(chr16_pos, mapping)
            
            if not gene:
                continue
            
            gene_pos = chr16_to_gene_coord(chr16_pos, gene, mapping)
            if gene_pos is None:
                continue
            
            ref = variant['ref'].upper()
            alt = variant['alt'].upper()
            
            # Note: ref.fa sequences are stored in plus-strand orientation (not reverse-complemented)
            # So we don't need to reverse-complement GIAB variants
            # Both GIAB VCF and ref.fa use plus-strand coordinates
            
            variants.add((gene, gene_pos, ref, alt))
    
    return variants


def load_paradism_variants(paradism_vcf: Path) -> Set[Tuple]:
    """Load ParaDISM variants as (gene, gene_pos, ref, alt) tuples."""
    variants = set()
    
    open_func = gzip.open if paradism_vcf.suffix == '.gz' else open
    mode = 'rt' if paradism_vcf.suffix == '.gz' else 'r'
    
    with open_func(paradism_vcf, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            variant = parse_vcf_line(line)
            if not variant:
                continue
            
            gene = variant['chrom']  # ParaDISM uses gene name as chrom
            gene_pos = variant['pos']
            ref = variant['ref'].upper()
            alt = variant['alt'].upper()
            
            # Don't normalize ref/alt - keep original order for matching
            variants.add((gene, gene_pos, ref, alt))
    
    return variants


def compare_variants(truth: Set[Tuple], called: Set[Tuple], total_positions: int = None, debug: bool = False) -> Dict[str, float]:
    """Compare called variants to truth set.
    
    Args:
        truth: Set of truth variants (gene, pos, ref, alt)
        called: Set of called variants (gene, pos, ref, alt)
        total_positions: Total number of positions that could have variants (for specificity)
        debug: Print debug information
    """
    tp_set = truth & called  # True Positives: in both sets
    fp_set = called - truth  # False Positives: in called but not truth
    fn_set = truth - called  # False Negatives: in truth but not called
    
    tp = len(tp_set)
    fp = len(fp_set)
    fn = len(fn_set)
    
    if debug:
        print(f"  TP variants: {sorted(tp_set)[:5]}...", file=sys.stderr)
        print(f"  FP variants: {sorted(fp_set)[:5]}...", file=sys.stderr)
        print(f"  FN variants: {sorted(fn_set)[:5]}...", file=sys.stderr)
    
    # Standard metrics
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    
    # Specificity: TN / (TN + FP)
    # TN = True Negatives = positions that are NOT variants in truth AND NOT called
    # For variant calling, we approximate this as:
    # - If we know total positions: TN = total_positions - TP - FP - FN
    # - Otherwise, we use: Specificity = 1 - (FP / (FP + estimated_TN))
    #   where estimated_TN = number of positions in called genes that don't have variants
    
    # Calculate specificity
    # Approach: TN = positions in gene regions without variants
    # We estimate this by considering all called positions
    if total_positions:
        tn = total_positions - tp - fp - fn
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    else:
        # Estimate: assume most positions are negatives
        # Use a conservative estimate based on gene lengths
        # For now, calculate without TN (set to None)
        specificity = None
    
    # Alternative: Calculate specificity using only called positions
    # Specificity = (positions called correctly as non-variant) / (all non-variant positions)
    # Since we don't have explicit non-variant positions, we use:
    # Specificity ≈ 1 - (FP rate) = 1 - (FP / (FP + estimated_non_variant_positions))
    # For practical purposes, we'll calculate it as: TN / (TN + FP) where TN is estimated
    
    # Estimate TN: positions in genes that could have variants but don't
    # We'll use the number of unique positions in called variants as a proxy
    called_positions = {(g, p) for g, p, _, _ in called}
    truth_positions = {(g, p) for g, p, _, _ in truth}
    
    # TN = positions that are NOT in truth AND NOT in called (but could be)
    # This is hard to calculate without knowing all possible positions
    # For now, we'll use a simplified metric:
    # Specificity = 1 - (FP / (FP + TP + FN))  # FP rate among all variants
    # Or better: Specificity = TN / (TN + FP) where we estimate TN
    
    # Simplified specificity calculation:
    # Assume that positions not in truth are negatives
    # TN ≈ positions in called genes that are not variants
    # We'll calculate it as: Specificity = 1 - (FP / total_called_positions)
    # But this isn't quite right either
    
    # Most practical approach for variant calling:
    # Calculate specificity based on the assumption that most positions are non-variant
    # Specificity = (correctly identified non-variants) / (all non-variants)
    # Since we only have variant positions, we'll use:
    # Specificity = 1 - (FP / (FP + estimated_TN))
    
    # For now, let's use a simpler metric:
    # Negative Predictive Value (NPV) = TN / (TN + FN)
    # But we don't have TN either
    
    # Best approach: Calculate specificity using gene lengths
    # TN = total positions in genes - TP - FP - FN
    # But we need total positions per gene
    
    return {
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'specificity': specificity,
        'total_truth': len(truth),
        'total_called': len(called),
        'tp_set': tp_set,
        'fp_set': fp_set,
        'fn_set': fn_set
    }


def calculate_specificity_with_gene_lengths(truth: Set[Tuple], called: Set[Tuple], mapping: Dict[str, Dict]) -> float:
    """Calculate specificity using gene lengths to estimate true negatives."""
    # Get all positions that could have variants (all positions in genes)
    total_positions = sum(info['length'] for info in mapping.values())
    
    # Get variant positions
    truth_positions = {(g, p) for g, p, _, _ in truth}
    called_positions = {(g, p) for g, p, _, _ in called}
    
    # TN = positions that are NOT variants in truth AND NOT called
    # This is an approximation since we don't know all possible variant positions
    # But we can estimate: TN ≈ total_positions - unique_variant_positions
    
    # Count unique positions with variants
    all_variant_positions = truth_positions | called_positions
    
    # Estimate TN: positions without variants
    # This is approximate since not all positions can have variants
    # But for large genes, most positions are non-variant
    tn_estimate = total_positions - len(all_variant_positions)
    
    # FP positions
    fp_positions = called_positions - truth_positions
    fp_count = len(fp_positions)
    
    # Specificity = TN / (TN + FP)
    if (tn_estimate + fp_count) > 0:
        specificity = tn_estimate / (tn_estimate + fp_count)
    else:
        specificity = 0.0
    
    return specificity


def main():
    parser = argparse.ArgumentParser(
        description="Compare ParaDISM variants to GIAB ground truth"
    )
    parser.add_argument('--giab-vcf', required=True, type=Path,
                        help='GIAB truth VCF file (chr16 coordinates)')
    parser.add_argument('--mapping', required=True, type=Path,
                        help='Coordinate mapping TSV file')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='ParaDISM output directory (giab_hg002_output)')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output TSV file with comparison results')
    
    args = parser.parse_args()
    
    # Load coordinate mapping
    mapping = load_coordinate_mapping(args.mapping)
    print(f"Loaded mapping for {len(mapping)} genes", file=sys.stderr)
    
    # Load GIAB truth variants
    print(f"Loading GIAB truth variants from {args.giab_vcf}...", file=sys.stderr)
    truth_variants = load_giab_variants(args.giab_vcf, mapping)
    print(f"Loaded {len(truth_variants)} GIAB truth variants", file=sys.stderr)
    
    # Compare for each aligner and method
    results = []
    aligners = ['minimap2', 'bwa-mem2', 'bowtie2']
    methods = [
        ('baseline', 'baseline/{aligner}/variant_calling/variants.vcf.gz'),
        ('final', '{aligner}/final_outputs/variant_calling/variants.vcf.gz')
    ]
    
    for aligner in aligners:
        for method_name, vcf_pattern in methods:
            vcf_path = args.output_dir / vcf_pattern.format(aligner=aligner)
            
            if not vcf_path.exists():
                print(f"Warning: {vcf_path} not found, skipping", file=sys.stderr)
                continue
            
            print(f"Loading {method_name} variants for {aligner}...", file=sys.stderr)
            called_variants = load_paradism_variants(vcf_path)
            print(f"  Loaded {len(called_variants)} variants", file=sys.stderr)
            
            # Debug: show gene distribution
            truth_genes = defaultdict(int)
            for g, _, _, _ in truth_variants:
                truth_genes[g] += 1
            called_genes = defaultdict(int)
            for g, _, _, _ in called_variants:
                called_genes[g] += 1
            print(f"  Truth genes: {dict(truth_genes)}", file=sys.stderr)
            print(f"  Called genes: {dict(called_genes)}", file=sys.stderr)
            
            comparison = compare_variants(truth_variants, called_variants, debug=(aligner=='minimap2' and method_name=='final'))
            
            # Calculate specificity using gene lengths
            specificity = calculate_specificity_with_gene_lengths(truth_variants, called_variants, mapping)
            comparison['specificity'] = specificity
            
            comparison['aligner'] = aligner
            comparison['method'] = method_name
            results.append(comparison)
    
    # Write results
    with open(args.output, 'w') as f:
        f.write("aligner\tmethod\ttp\tfp\tfn\ttotal_truth\ttotal_called\tprecision\trecall\tspecificity\tf1\n")
        for r in results:
            f.write(f"{r['aligner']}\t{r['method']}\t{r['tp']}\t{r['fp']}\t{r['fn']}\t"
                   f"{r['total_truth']}\t{r['total_called']}\t"
                   f"{r['precision']:.4f}\t{r['recall']:.4f}\t{r['specificity']:.4f}\t{r['f1']:.4f}\n")
    
    # Print summary
    print("\n=== Comparison Results ===", file=sys.stderr)
    print("Align        Method     TP     FP     FN     Precision  Recall     Specificity  F1", file=sys.stderr)
    print("-" * 85, file=sys.stderr)
    for r in results:
        print(f"{r['aligner']:12} {r['method']:8} {r['tp']:5}  {r['fp']:5}  {r['fn']:5}  "
              f"{r['precision']:8.4f}  {r['recall']:8.4f}  {r['specificity']:11.4f}  {r['f1']:6.4f}", file=sys.stderr)
    
    print(f"\nResults written to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
