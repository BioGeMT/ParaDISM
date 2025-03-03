#!/usr/bin/env python3
import sys
import argparse
import re

def clean_fasta(input_file, output_file):
    """
    Format FASTA file to have each sequence on a single line:
    1. Remove all whitespace from sequence lines
    2. Combine all sequence lines for each entry into a single line
    3. Preserve header lines intact
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_sequence = ""
        current_header = None
        
        for line in infile:
            line = line.strip()
            
            # Check if this is a header line
            if line.startswith('>'):
                # If we have a previous sequence to write out
                if current_header is not None:
                    outfile.write(current_header + '\n')
                    outfile.write(current_sequence + '\n')
                
                # Start new sequence
                current_header = line
                current_sequence = ""
            else:
                # Remove all whitespace and add to current sequence
                current_sequence += re.sub(r'\s', '', line)
        
        # Write the last sequence
        if current_header is not None:
            outfile.write(current_header + '\n')
            outfile.write(current_sequence + '\n')

def main():
    parser = argparse.ArgumentParser(description='Format FASTA files to have each sequence on a single line.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output formatted FASTA file')
    
    args = parser.parse_args()
    
    try:
        clean_fasta(args.input, args.output)
        print(f"Formatted FASTA file saved to {args.output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()