import os
import sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess

def process_files(tsv_path, r1_path, r2_path, output_dir):
    """
    Organize sequencing reads into gene-specific FASTQ files based on TSV mappings.
    Supports both paired-end (R1+R2) and single-end (R1 only) modes.

    Args:
        r2_path: Optional R2 file path (None for single-end)

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
                gene_mapping[read_name]  = gene

    # Load FASTQ records into dictionaries for quick lookup
    r1_records = SeqIO.to_dict(SeqIO.parse(r1_path, "fastq"))
    r2_records = SeqIO.to_dict(SeqIO.parse(r2_path, "fastq")) if r2_path else {}

    # Organize all records by gene in a single collection
    gene_collection = defaultdict(list)
    is_paired = r2_path is not None

    for read_name, gene in gene_mapping.items():
        # Add R1 record if exists
        if read_name in r1_records:
            r1 = r1_records[read_name]
            if is_paired:
                # Paired-end: add /1 suffix
                gene_collection[gene].append(SeqRecord(
                    r1.seq,
                    id=f"{read_name}/1",
                    description="",
                    letter_annotations=r1.letter_annotations
                ))
            else:
                # Single-end: no suffix
                gene_collection[gene].append(SeqRecord(
                    r1.seq,
                    id=read_name,
                    description="",
                    letter_annotations=r1.letter_annotations
                ))

        # Add R2 record if exists (only for paired-end)
        if is_paired and read_name in r2_records:
            r2 = r2_records[read_name]
            gene_collection[gene].append(SeqRecord(
                r2.seq,
                id=f"{read_name}/2",
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

        # Verbose output to stderr (goes to log)
        print(f"Created {gene}.fq with {len(records)} reads", file=sys.stderr)
        processed_genes.append(gene)

    return processed_genes

def create_bam_files(
    genes,
    ref_fasta,
    fastq_dir,
    output_dir,
    aligner='bowtie2',
    threads=4,
    minimap2_profile='short',
    is_paired=False,
):
    ref_db = {}
    for record in SeqIO.parse(ref_fasta, "fasta"):
        gene_name = record.id.split()[0]
        ref_db[gene_name] = str(record.seq)

    total_genes = len(genes)
    for idx, gene in enumerate(genes, 1):
        # Progress counter to stdout (visible under spinner)
        sys.stdout.write(f"\rProcessing genes: {idx}/{total_genes} ({gene})")
        sys.stdout.flush()

        if gene not in ref_db:
            print(f"Skipping {gene}: no reference sequence found", file=sys.stderr)
            continue

        # Verbose output to stderr (goes to log)
        print(f"Processing {gene} for alignment...", file=sys.stderr)
        fastq_file = os.path.join(fastq_dir, f"{gene}.fq")
        tmp_ref = os.path.join(output_dir, f"{gene}_ref.fa")
        index_base = os.path.join(output_dir, f"{gene}_index")
        sam_path = os.path.join(output_dir, f"{gene}.sam")
        bam_path = os.path.join(output_dir, f"{gene}.bam")
        sorted_bam = os.path.join(output_dir, f"{gene}.sorted.bam")

        # Verify files exist
        if not os.path.exists(fastq_file):
            print(f"Error: FASTQ file {fastq_file} does not exist", file=sys.stderr)
            continue
        if not os.path.exists(tmp_ref):
            with open(tmp_ref, "w") as f:
                f.write(f">{gene}\n{ref_db[gene]}\n")
        else:
            print(f"Warning: {tmp_ref} already exists, overwriting", file=sys.stderr)

        # Build index and align based on aligner choice
        if aligner == 'bowtie2':
            # Build bowtie2 index
            build_cmd = f"bowtie2-build {tmp_ref} {index_base}"
            print(f"Building bowtie2 index: {build_cmd}", file=sys.stderr)
            try:
                result = subprocess.run(build_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if result.stdout:
                    print(result.stdout, file=sys.stderr)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"Error building bowtie2 index for {gene}: {e.stderr}", file=sys.stderr)
                continue

            # Align with bowtie2
            if is_paired:
                align_cmd = f"bowtie2 --very-sensitive --interleaved {fastq_file} -x {index_base} -S {sam_path}"
            else:
                align_cmd = f"bowtie2 --very-sensitive -x {index_base} -U {fastq_file} -S {sam_path}"
            print(f"Running alignment: {align_cmd}", file=sys.stderr)
            try:
                result = subprocess.run(align_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if result.stdout:
                    print(result.stdout, file=sys.stderr)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"Error in bowtie2 alignment for {gene}: {e.stderr}", file=sys.stderr)
                continue

        elif aligner == 'bwa-mem2':
            # Build bwa-mem2 index
            build_cmd = f"bwa-mem2 index -p {index_base} {tmp_ref}"
            print(f"Building bwa-mem2 index: {build_cmd}", file=sys.stderr)
            try:
                result = subprocess.run(build_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if result.stdout:
                    print(result.stdout, file=sys.stderr)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"Error building bwa-mem2 index for {gene}: {e.stderr}", file=sys.stderr)
                continue

            # Align with bwa-mem2
            if is_paired:
                align_cmd = f"bwa-mem2 mem -t {threads} -p {index_base} {fastq_file} > {sam_path}"
            else:
                align_cmd = f"bwa-mem2 mem -t {threads} {index_base} {fastq_file} > {sam_path}"
            print(f"Running alignment: {align_cmd}", file=sys.stderr)
            try:
                result = subprocess.run(align_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if result.stdout:
                    print(result.stdout, file=sys.stderr)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"Error in bwa-mem2 alignment for {gene}: {e.stderr}", file=sys.stderr)
                continue

        elif aligner == 'minimap2':
            # Determine minimap2 preset
            preset_map = {
                'short': 'sr',
                'pacbio-hifi': 'map-hifi',
                'pacbio-clr': 'map-pb',
                'ont-q20': 'lr:hq',
                'ont-standard': 'map-ont',
            }
            preset = preset_map.get(minimap2_profile, 'sr')

            # Build minimap2 index
            index_mmi = f"{index_base}.mmi"
            build_cmd = f"minimap2 -x {preset} -d {index_mmi} {tmp_ref}"
            print(f"Building minimap2 ({minimap2_profile}) index: {build_cmd}", file=sys.stderr)
            try:
                result = subprocess.run(build_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if result.stdout:
                    print(result.stdout, file=sys.stderr)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"Error building minimap2 index for {gene}: {e.stderr}", file=sys.stderr)
                continue

            # Align with minimap2 (paired or single-end)
            align_cmd = f"minimap2 -ax {preset} --MD -t {threads} {index_mmi} {fastq_file} > {sam_path}"
            print(f"Running alignment: {align_cmd}", file=sys.stderr)
            try:
                result = subprocess.run(align_cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if result.stdout:
                    print(result.stdout, file=sys.stderr)
                if result.stderr:
                    print(result.stderr, file=sys.stderr)
            except subprocess.CalledProcessError as e:
                print(f"Error in minimap2 alignment for {gene}: {e.stderr}", file=sys.stderr)
                continue

        else:
            print(f"Error: Unsupported aligner '{aligner}'", file=sys.stderr)
            continue

        # Convert to BAM
        try:
            subprocess.run(f"samtools view -b {sam_path} > {bam_path}", shell=True, check=True)
            subprocess.run(f"samtools sort -o {sorted_bam} {bam_path}", shell=True, check=True)
            subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error in SAM/BAM conversion for {gene}: {e}", file=sys.stderr)
            continue

        # Cleanup index files based on aligner
        os.remove(tmp_ref)
        os.remove(sam_path)
        os.remove(bam_path)

        if aligner == 'bowtie2':
            for ext in ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']:
                index_file = f"{index_base}{ext}"
                if os.path.exists(index_file):
                    os.remove(index_file)
        elif aligner == 'bwa-mem2':
            for ext in ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']:
                index_file = f"{index_base}{ext}"
                if os.path.exists(index_file):
                    os.remove(index_file)
        elif aligner == 'minimap2':
            index_mmi = f"{index_base}.mmi"
            if os.path.exists(index_mmi):
                os.remove(index_mmi)

    # Print final newline to complete progress counter
    sys.stdout.write(f"\rProcessing genes: {total_genes}/{total_genes} (complete)\n")
    sys.stdout.flush()

def main(
    tsv_file,
    r1_file,
    r2_file,
    ref_fasta,
    fastq_dir,
    bam_dir,
    aligner='bowtie2',
    threads=4,
    minimap2_profile='short',
):
    # Ensure the output directories exist
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(bam_dir, exist_ok=True)

    # Step 1: Process files to create gene-specific FASTQs
    processed_genes = process_files(tsv_file, r1_file, r2_file, fastq_dir)

    # Step 2: Create BAM files from the gene-specific FASTQs
    create_bam_files(
        processed_genes,
        ref_fasta,
        fastq_dir,
        bam_dir,
        aligner,
        threads,
        minimap2_profile,
        is_paired=r2_file is not None,
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create gene-specific FASTQ and BAM files from read mappings (supports both paired-end and single-end).')
    parser.add_argument('--tsv', required=True, help='Input TSV file mapping read names to genes')
    parser.add_argument('--r1', required=True, help='Forward reads FASTQ file (or single-end reads)')
    parser.add_argument('--r2', required=False, default=None, help='Reverse reads FASTQ file (optional, for paired-end)')
    parser.add_argument('--ref', required=True, dest='ref_fasta', help='Reference sequences FASTA file')
    parser.add_argument('--fastq-dir', default='mapped_fastq', help='Output directory for gene-specific FASTQ files')
    parser.add_argument('--bam-dir', default='bam_files', help='Output directory for BAM files')
    parser.add_argument('--aligner', default='bowtie2', choices=['bowtie2', 'bwa-mem2', 'minimap2'], help='Aligner to use (default: bowtie2)')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads for alignment (default: 4)')
    parser.add_argument('--minimap2-profile', default='short',
                        choices=['short', 'pacbio-hifi', 'pacbio-clr', 'ont-q20', 'ont-standard'],
                        help='Minimap2 profile (default: short)')

    args = parser.parse_args()

    main(
        tsv_file=args.tsv,
        r1_file=args.r1,
        r2_file=args.r2,
        ref_fasta=args.ref_fasta,
        fastq_dir=args.fastq_dir,
        bam_dir=args.bam_dir,
        aligner=args.aligner,
        threads=args.threads,
        minimap2_profile=args.minimap2_profile
    )
