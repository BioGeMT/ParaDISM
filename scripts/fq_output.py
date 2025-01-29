import os
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def process_files(tsv_path, r1_path, r2_path, output_dir):
    # Read TSV and preserve order of mapped reads
    ordered_reads = []  # Read names in TSV order (non-NONE entries)
    gene_mapping = {}   # {read_name: gene}

    with open(tsv_path) as tsv:
        header = tsv.readline()  # Skip header
        for line in tsv:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            read_name = parts[0]
            gene = parts[1] if len(parts) > 1 else "NONE"
            if gene != "NONE":
                ordered_reads.append(read_name)
                gene_mapping[read_name] = gene

    # Load FASTQ records into dictionaries for quick lookup
    r1_records = SeqIO.to_dict(SeqIO.parse(r1_path, "fastq"))
    r2_records = SeqIO.to_dict(SeqIO.parse(r2_path, "fastq"))

    # Organize records by gene with TSV order and +/- pairing
    gene_collection = defaultdict(list)
    
    for read_name in ordered_reads:
        gene = gene_mapping[read_name]
        
        # Add R1 (+) record if exists
        if read_name in r1_records:
            r1 = r1_records[read_name]
            gene_collection[gene].append(SeqRecord(
                r1.seq,
                id=f"{read_name}+",
                description="",
                letter_annotations=r1.letter_annotations
            ))
        
        # Add R2 (-) record if exists
        if read_name in r2_records:
            r2 = r2_records[read_name]
            gene_collection[gene].append(SeqRecord(
                r2.seq,
                id=f"{read_name}-",
                description="",
                letter_annotations=r2.letter_annotations
            ))

    # Write output files to the specified directory
    for gene, records in gene_collection.items():
        output_file = os.path.join(output_dir, f"{gene}.fq")
        with open(output_file, "w") as out_f:
            SeqIO.write(records, out_f, "fastq")
        print(f"Created {output_file} with {len(records)} reads")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Organize sequencing reads into gene-specific FASTQ files based on TSV mappings.")
    parser.add_argument("--tsv", required=True, help="Input TSV file containing read-to-gene mappings.")
    parser.add_argument("--r1", required=True, help="Input forward reads (R1) FASTQ file.")
    parser.add_argument("--r2", required=True, help="Input reverse reads (R2) FASTQ file.")
    parser.add_argument("-o", "--output-dir", default="mapped_fastq", help="Output directory for gene-specific FASTQ files. Defaults to current directory.")
    
    args = parser.parse_args()
    
    # Ensure the output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    process_files(args.tsv, args.r1, args.r2, args.output_dir)