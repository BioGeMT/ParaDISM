import pandas as pd
import argparse
import os
import glob
import re


def find_unique_mapping_in_file(df, genes):
    # Check mapped genes and determine uniqueness
    plus_rows = df[df['Read_Name'].str.contains(r'\+')]
    minus_rows = df[df['Read_Name'].str.contains(r'-')]
    
    mapped_genes = []
    for gene in genes:
        case1 = (
            any((plus_rows[f'{gene}_C1_Pass'] == 1) & (plus_rows[f'{gene}_C2_Fail'] == 0)) and
            any(minus_rows[f'{gene}_C2_Fail'] == 0)
        )
        
        case2 = (
            any((minus_rows[f'{gene}_C1_Pass'] == 1) & (minus_rows[f'{gene}_C2_Fail'] == 0)) and
            any(plus_rows[f'{gene}_C2_Fail'] == 0)
        )
        
        if case1 or case2:
            mapped_genes.append(gene)
    
    return mapped_genes[0] if len(mapped_genes) == 1 else "NONE"

def natural_sort_key(s):
    # Numerical natural sorting
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', s)]


def find_unique_mappings(input_folder, output_file):
    # Process each results file in the input folder
    tsv_files = glob.glob(os.path.join(input_folder, '*_results.tsv'))
    if not tsv_files:
        print(f"No *_results.tsv files found in {input_folder}")
        return

    total_files = len(tsv_files)
    counter = 0

    with open(output_file, 'w') as f:
        f.write('Read_Name\tUniquely_Mapped\n')

        for file in tsv_files:
            try:
                read_name = os.path.basename(file).replace('_results.tsv', '')
                df = pd.read_csv(file, sep='\t')
                genes = [col.replace('_C1_Pass', '') for col in df.columns if col.endswith('_C1_Pass')]
                unique_gene = find_unique_mapping_in_file(df, genes)
                f.write(f'{read_name}\t{unique_gene}\n')
                counter += 1
                print(f"Processed {counter}/{total_files} files", end='\r')

            except Exception as e:
                print(f"Error processing {file}: {str(e)}")
                continue

    # Sort the output file numerically
    df = pd.read_csv(output_file, sep='\t')
    df['Read_Name'] = df['Read_Name'].astype(str)
    df['Sort_Key'] = df['Read_Name'].apply(natural_sort_key)
    df.sort_values(by='Sort_Key', inplace=True)
    df.drop('Sort_Key', axis=1, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Find uniquely mapped reads from folder of refinement results')
    parser.add_argument('--input_folder', required=True, help='Input folder containing *_results.tsv files')
    parser.add_argument('--output_file', required=True, help='Output TSV file for unique mapping results')
    
    args = parser.parse_args()
    find_unique_mappings(args.input_folder, args.output_file)

if __name__ == "__main__":
    main()