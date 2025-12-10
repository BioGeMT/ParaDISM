#!/usr/bin/env python3
"""
Extract assignments TSV from final FASTQ outputs.
This script reads the gene-specific FASTQ files created by the pipeline
and reconstructs the assignments dictionary.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO


def extract_assignments_from_fastq(fastq_dir: Path) -> dict[str, str]:
    """
    Extract assignments from gene-specific FASTQ files.
    
    Returns:
        dict: {read_name: gene_assignment} where gene_assignment is gene name or "NONE"
    """
    assignments = {}
    
    # Find all FASTQ files in the directory
    fastq_files = list(fastq_dir.glob("*.fq")) + list(fastq_dir.glob("*.fastq"))
    
    if not fastq_files:
        print(f"Warning: No FASTQ files found in {fastq_dir}", file=sys.stderr)
        return {}
    
    # Process each gene-specific FASTQ file
    for fastq_file in fastq_files:
        # Extract gene name from filename (e.g., "prefix_GENE.fq" or "GENE.fq")
        gene_name = fastq_file.stem
        if "_" in gene_name:
            # Remove prefix if present
            gene_name = gene_name.split("_", 1)[1] if gene_name.split("_", 1)[0] else gene_name
        
        # Read all reads from this FASTQ file
        for record in SeqIO.parse(fastq_file, "fastq"):
            # Normalize read ID (remove /1 or /2 suffix)
            read_id = record.id
            if read_id.endswith("/1") or read_id.endswith("/2"):
                base_read_id = read_id[:-2]
            else:
                base_read_id = read_id
            
            # Assign this read to the gene
            if base_read_id in assignments:
                # Already assigned - check for conflicts
                if assignments[base_read_id] != gene_name:
                    # Conflict: read appears in multiple gene files
                    assignments[base_read_id] = "NONE"
            else:
                assignments[base_read_id] = gene_name
    
    return assignments


def write_assignments_tsv(assignments: dict[str, str], output_path: Path):
    """Write assignments dictionary to TSV file."""
    with open(output_path, 'w') as f:
        f.write("read_name\tgene_assignment\n")
        for read_name, gene in sorted(assignments.items()):
            f.write(f"{read_name}\t{gene}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Extract assignments TSV from final FASTQ outputs',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--fastq-dir', required=True,
                        help='Directory containing gene-specific FASTQ files (e.g., final_outputs/{prefix}_fastq)')
    parser.add_argument('--output', required=True,
                        help='Output TSV file path')
    parser.add_argument('--original-fastq', required=False,
                        help='Original R1 FASTQ file (to include NONE reads)')
    parser.add_argument('--original-fastq-r2', required=False,
                        help='Original R2 FASTQ file (optional, for paired-end)')

    args = parser.parse_args()

    fastq_dir = Path(args.fastq_dir)
    output_path = Path(args.output)

    if not fastq_dir.exists():
        print(f"Error: FASTQ directory not found: {fastq_dir}", file=sys.stderr)
        sys.exit(1)

    # Extract assignments from FASTQ files
    assignments = extract_assignments_from_fastq(fastq_dir)
    
    # If original FASTQ files provided, add NONE reads
    if args.original_fastq:
        original_r1 = Path(args.original_fastq)
        if original_r1.exists():
            for record in SeqIO.parse(original_r1, "fastq"):
                read_id = record.id
                if read_id.endswith("/1") or read_id.endswith("/2"):
                    base_read_id = read_id[:-2]
                else:
                    base_read_id = read_id
                
                if base_read_id not in assignments:
                    assignments[base_read_id] = "NONE"
            
            # Also check R2 if provided
            if args.original_fastq_r2:
                original_r2 = Path(args.original_fastq_r2)
                if original_r2.exists():
                    for record in SeqIO.parse(original_r2, "fastq"):
                        read_id = record.id
                        if read_id.endswith("/1") or read_id.endswith("/2"):
                            base_read_id = read_id[:-2]
                        else:
                            base_read_id = read_id
                        
                        if base_read_id not in assignments:
                            assignments[base_read_id] = "NONE"

    # Write TSV
    write_assignments_tsv(assignments, output_path)
    print(f"Assignments written to: {output_path}", file=sys.stderr)
    print(f"Total reads: {len(assignments)}", file=sys.stderr)
    print(f"Assigned reads: {sum(1 for g in assignments.values() if g != 'NONE')}", file=sys.stderr)
    print(f"NONE reads: {sum(1 for g in assignments.values() if g == 'NONE')}", file=sys.stderr)


if __name__ == '__main__':
    main()
