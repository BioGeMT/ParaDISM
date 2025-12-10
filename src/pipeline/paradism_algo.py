import argparse
import sys
from typing import Dict, List
from collections import defaultdict
from Bio import AlignIO
from Bio.Align.sam import AlignmentIterator


def load_msa(msa_fasta_path: str):
    """
    Returns:
        msa: MultipleSeqAlignment object (for base lookup)
        ref_to_msa: {gene: {ref_pos (1-based): msa_column}}
        gene_names: list of gene names
    """
    msa = AlignIO.read(msa_fasta_path, 'fasta')
    gene_names = [record.id for record in msa]
    msa_length = msa.get_alignment_length()

    ref_to_msa = {gene: {} for gene in gene_names}
    ref_positions = {gene: 0 for gene in gene_names}

    for msa_col in range(msa_length):
        for i, gene in enumerate(gene_names):
            base = msa[i, msa_col]
            if base != '-':
                ref_positions[gene] += 1
                ref_to_msa[gene][ref_positions[gene]] = msa_col

    return msa, ref_to_msa, gene_names


def process_read():
    pass


def process_sam(sam_path):
    alignments = AlignmentIterator(sam_path)


def main():
    parser = argparse.ArgumentParser(
        description='ParaDISM: Paralog-aware read assignment using SNPs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sam', required=True,
                        help='Input SAM/BAM file')
    parser.add_argument('--msa', required=True,
                        help='Input MSA FASTA file')
    parser.add_argument('--output', required=True,
                        help='Output file for unique mapping results')
    parser.add_argument('--single-end', action='store_true',
                        help='Treat input as single-end reads')

    args = parser.parse_args()

    msa, ref_to_msa, gene_names = load_msa(args.msa)


if __name__ == '__main__':
    import time
    start_time = time.perf_counter()

    main()

    end_time = time.perf_counter()
    print(f"ParaDISM execution time: {end_time - start_time:.6f} seconds", file=sys.stderr)
