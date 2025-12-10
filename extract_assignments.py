#!/usr/bin/env python3
"""
Extract assignments TSV from ParaDISM pipeline output.
This script reads the SAM and MSA files created during pipeline execution
and outputs assignments.tsv without re-running the full pipeline.
"""

import argparse
import sys
from pathlib import Path

SRC_DIR = Path(__file__).parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pipeline.paradism_algo import load_msa, process_sam_to_dict
from pipeline.paradism_algo_NEW import load_msa as load_msa_new, process_sam_to_dict_simple as process_sam_to_dict_new


def write_assignments_tsv(assignments: dict[str, str], output_path: Path):
    """Write assignments dictionary to TSV file."""
    with open(output_path, 'w') as f:
        f.write("read_name\tgene_assignment\n")
        for read_name, gene in sorted(assignments.items()):
            f.write(f"{read_name}\t{gene}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Extract assignments TSV from ParaDISM pipeline output',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sam', required=True,
                        help='Input SAM file (from pipeline output)')
    parser.add_argument('--msa', required=True,
                        help='Input MSA file (from pipeline output)')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path')
    parser.add_argument('--version', choices=['original', 'NEW'], default='original',
                        help='Which algorithm version to use')

    args = parser.parse_args()

    sam_path = Path(args.sam)
    msa_path = Path(args.msa)
    output_path = Path(args.output)

    if not sam_path.exists():
        print(f"Error: SAM file not found: {sam_path}", file=sys.stderr)
        sys.exit(1)
    
    if not msa_path.exists():
        print(f"Error: MSA file not found: {msa_path}", file=sys.stderr)
        sys.exit(1)

    # Load MSA and process SAM
    if args.version == 'NEW':
        msa_obj, seq_to_aln, gene_names = load_msa_new(str(msa_path))
        assignments = process_sam_to_dict_new(str(sam_path), msa_obj, seq_to_aln, gene_names)
    else:
        msa_obj, seq_to_aln, gene_names = load_msa(str(msa_path))
        assignments = process_sam_to_dict(str(sam_path), msa_obj, seq_to_aln, gene_names)

    # Write TSV
    write_assignments_tsv(assignments, output_path)
    print(f"Assignments written to: {output_path}", file=sys.stderr)
    print(f"Total reads: {len(assignments)}", file=sys.stderr)
    print(f"Assigned reads: {sum(1 for g in assignments.values() if g != 'NONE')}", file=sys.stderr)


if __name__ == '__main__':
    main()
