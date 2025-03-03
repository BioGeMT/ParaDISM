#!/usr/bin/env python3
"""
FASTQ Filter Script

This script filters paired-end FASTQ files (R1 and R2) based on a TSV file that indicates 
whether each read is uniquely mapped or not. Only reads that don't have "NONE" in the 
Uniquely_Mapped column will be kept.

Usage:
    python filter_fastq.py --r1 input_r1.fastq --r2 input_r2.fastq --tsv mapping.tsv --out_prefix filtered

Arguments:
    --r1: Path to the R1 FASTQ file
    --r2: Path to the R2 FASTQ file
    --tsv: Path to the TSV file with mapping information
    --out_prefix: Prefix for output filtered FASTQ files (will create {prefix}_r1.fastq and {prefix}_r2.fastq)
"""

import argparse
import re
import sys

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Filter FASTQ files based on mapping status from TSV file.')
    parser.add_argument('--r1', required=True, help='Path to R1 FASTQ file')
    parser.add_argument('--r2', required=True, help='Path to R2 FASTQ file')
    parser.add_argument('--tsv', required=True, help='Path to TSV file with mapping information')
    parser.add_argument('--out_prefix', required=True, help='Prefix for output filtered FASTQ files')
    
    return parser.parse_args()

def read_mapping_status(tsv_file):
    """
    Read the TSV file and create a dictionary mapping read names to their mapping status.
    Returns a set of read names that should be kept (not mapped as "NONE").
    """
    reads_to_keep = set()
    
    with open(tsv_file, 'r') as f:
        # Skip header line
        next(f)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                read_name = parts[0]
                mapping_status = parts[1]
                
                # Keep reads that don't have "NONE" as mapping status
                if mapping_status != "NONE":
                    reads_to_keep.add(read_name)
    
    return reads_to_keep

def filter_fastq_file(input_fastq, output_fastq, reads_to_keep):
    """Filter a FASTQ file to keep only reads in the reads_to_keep set."""
    with open(input_fastq, 'r') as infile, open(output_fastq, 'w') as outfile:
        while True:
            # Read the four lines that make up a FASTQ record
            header = infile.readline().strip()
            if not header:
                break  # End of file
                
            sequence = infile.readline().strip()
            plus_line = infile.readline().strip()
            quality = infile.readline().strip()
            
            # Extract read name from the header
            # The format is @Read0 PKD1P1; 10514-10797, so we need to get "Read0"
            match = re.match(r'^@(\S+)', header)
            if match:
                read_name = match.group(1)
                
                # If this read should be kept, write it to the output file
                if read_name in reads_to_keep:
                    outfile.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")

def main():
    args = parse_arguments()
    
    # Get the set of reads to keep (those not mapped as "NONE")
    reads_to_keep = read_mapping_status(args.tsv)
    
    print(f"Found {len(reads_to_keep)} reads to keep out of those in the TSV file.")
    
    # Filter R1 FASTQ file
    r1_output = f"{args.out_prefix}_r1.fastq"
    filter_fastq_file(args.r1, r1_output, reads_to_keep)
    
    # Filter R2 FASTQ file
    r2_output = f"{args.out_prefix}_r2.fastq"
    filter_fastq_file(args.r2, r2_output, reads_to_keep)
    
    print(f"Filtered FASTQ files have been written to {r1_output} and {r2_output}")

if __name__ == "__main__":
    main()