import os
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess

def process_files(tsv_path, r1_path, r2_path, output_dir):
    """
    Organize sequencing reads into gene-specific FASTQ files based on TSV mappings.
    
    Returns a list of genes processed.
    """
    # Read TSV and preserve order of mapped reads
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
                gene_mapping[read_name] = gene

    # Load FASTQ records into dictionaries for quick lookup
    r1_records = SeqIO.to_dict(SeqIO.parse(r1_path, "fastq"))
    r2_records = SeqIO.to_dict(SeqIO.parse(r2_path, "fastq"))

    # Organize all records by gene in a single collection
    gene_collection = defaultdict(list)
    
    for read_name, gene in gene_mapping.items():
        # Add R1 record if exists
        if read_name in r1_records:
            r1 = r1_records[read_name]
            gene_collection[gene].append(SeqRecord(
                r1.seq,
                id=f"{read_name}/1",  # Add suffix to distinguish R1
                description="",
                letter_annotations=r1.letter_annotations
            ))
        
        # Add R2 record if exists
        if read_name in r2_records:
            r2 = r2_records[read_name]
            gene_collection[gene].append(SeqRecord(
                r2.seq,
                id=f"{read_name}/2",  # Add suffix to distinguish R2
                description="",
                letter_annotations=r2.letter_annotations
            ))

    processed_genes = []
    
    # Write output files to the specified directory, one file per gene
    for gene, records in gene_collection.items():
        if gene == "NONE":
            continue
            
        # Only process genes that have reads
        if not records:
            continue
            
        output_file = os.path.join(output_dir, f"{gene}.fq")
        
        with open(output_file, "w") as out_f:
            SeqIO.write(records, out_f, "fastq")
                
        print(f"Created {gene}.fq with {len(records)} reads")
        processed_genes.append(gene)
        
    return processed_genes

def create_bam_files(genes, ref_fasta, fastq_dir, output_dir):
    ref_db = {}
    for record in SeqIO.parse(ref_fasta, "fasta"):
        gene_name = record.id.split()[0]
        ref_db[gene_name] = str(record.seq)

    for gene in genes:
        if gene not in ref_db:
            print(f"Skipping {gene}: no reference sequence found")
            continue

        print(f"Processing {gene} for alignment...")
        fastq_file = os.path.join(fastq_dir, f"{gene}.fq")
        tmp_ref = os.path.join(output_dir, f"{gene}_ref.fa")
        bt2_index_base = os.path.join(output_dir, f"{gene}_index")
        sam_path = os.path.join(output_dir, f"{gene}.sam")
        bam_path = os.path.join(output_dir, f"{gene}.bam")
        sorted_bam = os.path.join(output_dir, f"{gene}.sorted.bam")

        # Verify files exist
        if not os.path.exists(fastq_file):
            print(f"Error: FASTQ file {fastq_file} does not exist")
            continue
        if not os.path.exists(tmp_ref):
            with open(tmp_ref, "w") as f:
                f.write(f">{gene}\n{ref_db[gene]}\n")
        else:
            print(f"Warning: {tmp_ref} already exists, overwriting")

        # Build index
        build_cmd = f"bowtie2-build {tmp_ref} {bt2_index_base}"
        print(f"Building bowtie2 index: {build_cmd}")
        try:
            result = subprocess.run(build_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error building index for {gene}: {e.stderr}")
            continue

        # Align with bowtie2
        bowtie_cmd = f"bowtie2 --very-sensitive --end-to-end --interleaved {fastq_file} -x {bt2_index_base} -S {sam_path}"
        print(f"Running alignment: {bowtie_cmd}")
        try:
            result = subprocess.run(bowtie_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            print(result.stdout)
            if result.stderr:
                print(f"bowtie2 stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"Error in bowtie2 alignment for {gene}: {e.stderr}")
            continue

        # Convert to BAM
        try:
            subprocess.run(f"samtools view -b {sam_path} > {bam_path}", shell=True, check=True)
            subprocess.run(f"samtools sort -o {sorted_bam} {bam_path}", shell=True, check=True)
            subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error in SAM/BAM conversion for {gene}: {e}")
            continue

        # Cleanup
        os.remove(tmp_ref)
        os.remove(sam_path)
        os.remove(bam_path)
        for ext in ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']:
            index_file = f"{bt2_index_base}{ext}"
            if os.path.exists(index_file):
                os.remove(index_file)

def main(tsv_file, r1_file, r2_file, ref_fasta, fastq_dir, bam_dir):
    # Ensure the output directories exist
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(bam_dir, exist_ok=True)
    
    # Step 1: Process files to create gene-specific FASTQs
    processed_genes = process_files(tsv_file, r1_file, r2_file, fastq_dir)
    
    # Step 2: Create BAM files from the gene-specific FASTQs
    create_bam_files(processed_genes, ref_fasta, fastq_dir, bam_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create gene-specific FASTQ and BAM files from read mappings.')
    parser.add_argument('--tsv', required=True, help='Input TSV file mapping read names to genes')
    parser.add_argument('--r1', required=True, help='Forward reads FASTQ file')
    parser.add_argument('--r2', required=True, help='Reverse reads FASTQ file')
    parser.add_argument('--ref', required=True, dest='ref_fasta', help='Reference sequences FASTA file')
    parser.add_argument('--fastq-dir', default='mapped_fastq', help='Output directory for gene-specific FASTQ files')
    parser.add_argument('--bam-dir', default='bam_files', help='Output directory for BAM files')
    
    args = parser.parse_args()
    
    main(
        tsv_file=args.tsv,
        r1_file=args.r1,
        r2_file=args.r2,
        ref_fasta=args.ref_fasta,
        fastq_dir=args.fastq_dir,
        bam_dir=args.bam_dir
    )