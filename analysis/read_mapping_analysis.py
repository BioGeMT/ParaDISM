import pandas as pd
import os
import argparse

def parse_fastq(fastq_file):
    read_to_origin = {}
    with open(fastq_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                parts = line.strip().split()
                read_name = parts[0][1:]
                gene_info = parts[1].split(';')[0]
                read_to_origin[read_name] = gene_info
    return read_to_origin

def parse_sam(sam_file):
    read_to_mapped = {}
    with open(sam_file, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                parts = line.strip().split('\t')
                read_name = parts[0]
                mapped_gene = parts[2]
                if mapped_gene != '*':
                    read_to_mapped[read_name] = mapped_gene
    return read_to_mapped

def process_files(tsv_file, fastq_file, sam_file, output_folder, output_prefix):
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_prefix)
    
    df_tsv = pd.read_csv(tsv_file, sep='\t')
    origin_dict = parse_fastq(fastq_file)
    sam_dict = parse_sam(sam_file)
    
    tsv_read_names = set(df_tsv['Read_Name'])
    fastq_read_names = set(origin_dict.keys())
    valid_reads = tsv_read_names.intersection(fastq_read_names)
    
    print(f"Reads in TSV: {len(tsv_read_names)}")
    print(f"Reads in FASTQ: {len(fastq_read_names)}")
    print(f"Valid reads for analysis: {len(valid_reads)}")
    
    mapper_results = []
    for _, row in df_tsv.iterrows():
        read_name = row['Read_Name']
        if read_name in valid_reads:
            mapped_gene = row['Uniquely_Mapped']
            origin_gene = origin_dict[read_name]
            mapper_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 'Mapped_Gene': mapped_gene})
    mapper_df = pd.DataFrame(mapper_results)
    
    bowtie2_results = []
    for read_name in valid_reads:
        origin_gene = origin_dict[read_name]
        mapped_gene = sam_dict.get(read_name, 'NONE')
        bowtie2_results.append({'Read_Name': read_name, 'Origin_Gene': origin_gene, 'Mapped_Gene': mapped_gene})
    bowtie2_df = pd.DataFrame(bowtie2_results)
    
    mapper_cm = pd.crosstab(mapper_df['Origin_Gene'], mapper_df['Mapped_Gene'])
    bowtie2_cm = pd.crosstab(bowtie2_df['Origin_Gene'], bowtie2_df['Mapped_Gene'])
    
    mapper_cm.to_csv(f'{output_path}_mapper_confusion_matrix.tsv', sep='\t')
    bowtie2_cm.to_csv(f'{output_path}_bowtie2_confusion_matrix.tsv', sep='\t')
    
    print(f"Confusion matrices saved:")
    print(f"  Mapper: {output_path}_mapper_confusion_matrix.tsv")
    print(f"  Bowtie2: {output_path}_bowtie2_confusion_matrix.tsv")

def main():
    parser = argparse.ArgumentParser(description='Read mapping analysis - confusion matrices only')
    parser.add_argument('-t', '--tsv', required=True, help='TSV file with Mapper results')
    parser.add_argument('-f', '--fastq', required=True, help='FASTQ file with read origins')
    parser.add_argument('-s', '--sam', required=True, help='SAM file with Bowtie2 results')
    parser.add_argument('-o', '--output-dir', default='mapping_results', help='Output directory')
    parser.add_argument('-p', '--prefix', default='read_mapping', help='Output file prefix')
    
    args = parser.parse_args()
    
    process_files(args.tsv, args.fastq, args.sam, args.output_dir, args.prefix)

if __name__ == "__main__":
    main()