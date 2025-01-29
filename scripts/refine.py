import os
import pandas as pd
import glob
import argparse

def find_unique_mapping(read, msa):
   genes = list(msa.keys())
   sequence_length = len(read)

   match_matrix = {}
   for gene in genes:
       match_matrix[gene] = [1 if msa[gene][i].upper() == read[i].upper() else 0 
                           for i in range(sequence_length)]

   condition1_results = {}
   condition2_fail_pos = {}
   
   for gene in genes:
       c1_passed = False
       c1_position = None
       
       for i in range(sequence_length):
           target_match = match_matrix[gene][i]
           other_matches = [match_matrix[other][i] for other in genes if other != gene]
           
           if target_match == 1 and sum(other_matches) == 0:
               c1_passed = True
               c1_position = i
               break
               
       condition1_results[gene] = (c1_passed, c1_position)
       
       c2_fail_pos = None
       for i in range(sequence_length):
           if match_matrix[gene][i] == 0:
               other_matches = [match_matrix[other][i] for other in genes if other != gene]
               if sum(other_matches) > 0:
                   c2_fail_pos = i
                   break
                   
       condition2_fail_pos[gene] = c2_fail_pos

   return condition1_results, condition2_fail_pos

def process_tsv_file(file_path):
    
   df = pd.read_csv(file_path, sep='\t', dtype={'MSA_Position': int})
   read_name = os.path.basename(file_path)
   gene_cols = [col for col in df.columns if col.endswith('_Base') and not col.startswith('Read')]
   msa = {}
   
   read_seq = ''.join(df['Read_Base'])
   
   for col in gene_cols:
       gene = col.replace('_Base', '')
       msa[gene] = ''.join(df[col])
   
   c1_results, c2_fail_pos = find_unique_mapping(read_seq, msa)
   results = {'Read_Name': read_name}
   pos_map = dict(zip(df.index, df['MSA_Position']))
   
   for gene in sorted(msa.keys()):
       passed, position = c1_results[gene]
       results[f'{gene}_C1_Pass'] = 1 if passed else 0
       results[f'{gene}_C1_Position'] = int(pos_map[position]) if passed else None
       
       fail_pos = c2_fail_pos[gene]
       results[f'{gene}_C2_Fail'] = 1 if fail_pos is not None else 0
       results[f'{gene}_C2_Fail_Position'] = int(pos_map[fail_pos]) if fail_pos is not None else None
   
   return results

def process_folder(folder_path, output_file):
   tsv_files = glob.glob(os.path.join(folder_path, '*.tsv'))
   
   if not tsv_files:
       print(f"No TSV files found in {folder_path}")
       return
   
   all_results = []
   for file_path in sorted(tsv_files):
       try:
           results = process_tsv_file(file_path)
           all_results.append(results)
       except Exception as e:
           print(f"Error processing {file_path}: {str(e)}")
   
   if not all_results:
       print(f"No results generated for {folder_path}")
       return
       
   results_df = pd.DataFrame(all_results)
   columns = ['Read_Name']
   genes = sorted(set(col.split('_')[0] for col in results_df.columns 
                 if col != 'Read_Name' and not col.startswith('Read_')))
   
   for gene in genes:
       columns.extend([f'{gene}_C1_Pass', f'{gene}_C2_Fail'])
   
   for gene in genes:
       columns.extend([f'{gene}_C1_Position', f'{gene}_C2_Fail_Position'])
   
   results_df = results_df[columns]
   os.makedirs(os.path.dirname(output_file), exist_ok=True)
   results_df.to_csv(output_file, sep='\t', index=False, na_rep='')

def process_multiple_folders(base_dir, output_dir):
   folders = [f.path for f in os.scandir(base_dir) if f.is_dir()]
   if not folders:
       folders = [base_dir]
   
   total_folders = len(folders)
   counter = 0
   
   for folder in sorted(folders):
       folder_name = os.path.basename(folder)
       output_file = os.path.join(output_dir, f"{folder_name}_results.tsv")
       process_folder(folder, output_file)
       counter += 1
       print(f"Folders processed: {counter}/{total_folders}", end='\r')
   
   print("\nComplete")

def main():
   parser = argparse.ArgumentParser(description='Process read mapping files from multiple folders')
   parser.add_argument('--reads_dir', required=True, help='Base directory containing folders with TSV files')
   parser.add_argument('--output_dir', required=True, help='Output directory for results files')
   
   args = parser.parse_args()
   os.makedirs(args.output_dir, exist_ok=True)
   process_multiple_folders(args.reads_dir, args.output_dir)

if __name__ == "__main__":
   main()