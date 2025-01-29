import os
import argparse
from Bio import SeqIO
from collections import defaultdict
import subprocess

def main(tsv_file, r1_file, r2_file, ref_fasta, output_dir):
    # Read TSV to map reads to genes (excluding NONE)
    read_to_gene = {}
    with open(tsv_file) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2 or parts[1] == "NONE":
                continue
            read_name, gene = parts[0], parts[1]
            read_to_gene[read_name] = gene

    # Read reference FASTA (gene names derived from first word of header)
    ref_db = {}
    for record in SeqIO.parse(ref_fasta, "fasta"):
        gene_name = record.id.split()[0]  # Use first part of header as gene name
        ref_db[gene_name] = str(record.seq)

    # Load all R1/R2 reads into memory for quick access
    r1_reads = SeqIO.to_dict(SeqIO.parse(r1_file, "fastq"))
    r2_reads = SeqIO.to_dict(SeqIO.parse(r2_file, "fastq"))

    # Group reads by gene
    gene_groups = defaultdict(list)
    for read_name, gene in read_to_gene.items():
        gene_groups[gene].append(read_name)

    # Prepare output directory
    os.makedirs(output_dir, exist_ok=True)

    # Process each gene
    for gene, read_names in gene_groups.items():
        if gene not in ref_db:
            print(f"Skipping {gene}: no reference sequence found")
            continue

        print(f"Processing {gene}...")

        # Create temporary reference FASTA for the gene
        tmp_ref = os.path.join(output_dir, f"{gene}_ref.fa")
        with open(tmp_ref, "w") as f:
            f.write(f">{gene}\n{ref_db[gene]}\n")

        # Extract R1/R2 reads for this gene
        tmp_r1 = os.path.join(output_dir, f"{gene}_R1.fq")
        tmp_r2 = os.path.join(output_dir, f"{gene}_R2.fq")
        
        with open(tmp_r1, "w") as f1, open(tmp_r2, "w") as f2:
            for name in read_names:
                if name in r1_reads:
                    SeqIO.write(r1_reads[name], f1, "fastq")
                if name in r2_reads:
                    SeqIO.write(r2_reads[name], f2, "fastq")

        # Skip if no reads found
        if os.path.getsize(tmp_r1) == 0 and os.path.getsize(tmp_r2) == 0:
            print(f"No reads for {gene}")
            os.remove(tmp_ref)
            continue

        # Align with minimap2 and create sorted BAM
        sam_path = os.path.join(output_dir, f"{gene}.sam")
        bam_path = os.path.join(output_dir, f"{gene}.bam")
        sorted_bam = os.path.join(output_dir, f"{gene}.sorted.bam")

        # Run minimap2 (paired-end alignment)
        minimap_cmd = f"minimap2 -ax sr {tmp_ref} {tmp_r1} {tmp_r2} > {sam_path}"
        subprocess.run(minimap_cmd, shell=True, check=True)

        # Convert SAM to sorted BAM
        subprocess.run(f"samtools view -b {sam_path} > {bam_path}", shell=True, check=True)
        subprocess.run(f"samtools sort -o {sorted_bam} {bam_path}", shell=True, check=True)
        subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)

        # Cleanup temporary files
        os.remove(tmp_ref)
        os.remove(tmp_r1)
        os.remove(tmp_r2)
        os.remove(sam_path)
        os.remove(bam_path)

    print("Done! Sorted BAM files are in:", output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create gene-specific BAM files from read mappings.')
    parser.add_argument('--tsv', required=True, help='Input TSV file mapping read names to genes')
    parser.add_argument('--r1', required=True, help='Forward reads FASTQ file')
    parser.add_argument('--r2', required=True, help='Reverse reads FASTQ file')
    parser.add_argument('--ref', required=True, dest='ref_fasta', help='Reference sequences FASTA file')
    parser.add_argument('--output', default='bam_files', help='Output directory for BAM files (default: gene_bams)')
    
    args = parser.parse_args()
    
    main(
        tsv_file=args.tsv,
        r1_file=args.r1,
        r2_file=args.r2,
        ref_fasta=args.ref_fasta,
        output_dir=args.output
    )