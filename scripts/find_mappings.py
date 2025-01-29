import pandas as pd
import argparse
import os
import glob
import re

def natural_sort_key(s):
    """Convert string with numbers into list of strings and integers for proper numerical sorting."""
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def find_unique_mapping_in_file(df, genes):
    """Check if read in this file uniquely maps to one gene."""
    # Split into + and - strands
    plus_rows = df[df['Read_Name'].str.contains(r'\+')]
    minus_rows = df[df['Read_Name'].str.contains(r'-')]
    
    # Check each gene
    mapped_genes = []
    for gene in genes:
        # Check both possible cases
        case1 = (
            # At least one + row with C1_Pass=1 and C2_Fail=0
            any((plus_rows[f'{gene}_C1_Pass'] == 1) & (plus_rows[f'{gene}_C2_Fail'] == 0)) and
            # At least one - row with C2_Fail=0
            any(minus_rows[f'{gene}_C2_Fail'] == 0)
        )
        
        case2 = (
            # At least one - row with C1_Pass=1 and C2_Fail=0
            any((minus_rows[f'{gene}_C1_Pass'] == 1) & (minus_rows[f'{gene}_C2_Fail'] == 0)) and
            # At least one + row with C2_Fail=0
            any(plus_rows[f'{gene}_C2_Fail'] == 0)
        )
        
        if case1 or case2:
            mapped_genes.append(gene)
    
    # If exactly one gene meets the criteria, it's unique
    return mapped_genes[0] if len(mapped_genes) == 1 else "NONE"

def find_unique_mappings(input_folder, output_file):
    """Process each results file in the input folder."""
    # Get all results TSV files
    tsv_files = glob.glob(os.path.join(input_folder, '*_results.tsv'))
    if not tsv_files:
        print(f"No *_results.tsv files found in {input_folder}")
        return

    with open(output_file, 'w') as f:
        # Write header line
        f.write('Read_Name\tUniquely_Mapped\n')

        for file in tsv_files:
            try:
                # Get read name from filename
                read_name = os.path.basename(file).replace('_results.tsv', '')

                # Read TSV file
                df = pd.read_csv(file, sep='\t')

                # Get list of genes (looking at C1_Pass columns)
                genes = [col.replace('_C1_Pass', '') for col in df.columns if col.endswith('_C1_Pass')]

                # Find unique mapping
                unique_gene = find_unique_mapping_in_file(df, genes)

                # Write result to file incrementally
                f.write(f'{read_name}\t{unique_gene}\n')
                print(f"Processed: {read_name}")

            except Exception as e:
                print(f"Error processing {file}: {str(e)}")
                continue

    # Sort the output file numerically by Read_Name
    df = pd.read_csv(output_file, sep='\t')
    df = df.sort_values(by='Read_Name', key=lambda x: x.map(natural_sort_key))
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Results written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Find uniquely mapped reads from folder of refinement results')
    parser.add_argument('input_folder', help='Input folder containing *_results.tsv files')
    parser.add_argument('output_file', help='Output TSV file for results')
    
    args = parser.parse_args()
    
    find_unique_mappings(args.input_folder, args.output_file)

if __name__ == "__main__":
    main()