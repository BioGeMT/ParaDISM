import argparse
import sys
from typing import Dict, List, Tuple
from collections import defaultdict
from Bio.Align.sam import AlignmentIterator
from Bio import AlignIO


def load_msa_mapping(msa_filepath: str) -> Tuple[Dict, Dict, List]:
    sequence_positions = {}
    sequence_bases = {}
    gene_names = []

    with open(msa_filepath, 'r') as f:
        header = next(f).strip().split('\t')

        for col in header:
            if col.endswith('_Position') and col != 'MSA_Position':
                gene_name = col.replace('_Position', '')
                gene_names.append(gene_name)
                sequence_positions[gene_name] = {}

        for line in f:
            fields = line.strip().split('\t')
            msa_pos = int(fields[0])
            sequence_bases[msa_pos] = {}

            for gene in gene_names:
                pos_idx = header.index(f'{gene}_Position')
                base_idx = header.index(f'{gene}_Base')

                pos = fields[pos_idx]
                base = fields[base_idx]
                sequence_bases[msa_pos][gene] = base

                if pos != '-':
                    sequence_positions[gene][int(pos)] = msa_pos

    return sequence_positions, sequence_bases, gene_names


def process_read():
    pass

def process_sam(sam_path, msa, msa_maps):
    alignments = AlignmentIterator(sam_path)
    qname_to_c1 = defaultdict(False)  # does the read pair pass c1?
    qname_to_c2 = defaultdict(True) # does the read pair pass c2?  
    for read in alignments:
        qname = read.sequences[0].id 
        c1, c2 = process_read(read, msa, msa_maps)
        if c1:
            qname_to_c1[qname] = c1
        if not c2:
            qname_to_c2[qname] = c2
            


def main():
    parser = argparse.ArgumentParser(
        description='ParaDISM: Paralog-aware read assignment using SNPs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sam', required=True,
                        help='Input SAM/BAM file')
    parser.add_argument('--msa', required=True,
                        help='Input MSA mapping TSV file')
    parser.add_argument('--output', required=True,
                        help='Output file for unique mapping results')
    parser.add_argument('--single-end', action='store_true',
                        help='Treat input as single-end reads')

    args = parser.parse_args()

    msa = AlignIO(args.msa, 'fasta')
    
    sequence_positions, sequence_bases, gene_names = load_msa_mapping(args.msa)


if __name__ == '__main__':
    import time
    start_time = time.perf_counter()

    main()

    end_time = time.perf_counter()
    print(f"ParaDISM execution time: {end_time - start_time:.6f} seconds", file=sys.stderr)
