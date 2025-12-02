import os
import sys
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import subprocess

def _base_read_id(read_id: str) -> str:
    """Normalise read ID to match SAM/TSV names.

    Aligners like BWA-MEM2 often strip trailing '/1' or '/2' from FASTQ
    read IDs when writing SAM records. Our TSV uses the SAM-style IDs,
    so we normalise FASTQ IDs by removing a single trailing '/1' or '/2'.
    """
    if read_id.endswith("/1") or read_id.endswith("/2"):
        return read_id[:-2]
    return read_id


def process_files(tsv_path, r1_path, r2_path, output_dir, prefix=""):
    """Split reads into per-gene FASTQs (paired or single-end)."""
    # Read TSV and preserve order of mapped reads
    gene_mapping = {}  # {read_name: (gene, strand_info)}

    with open(tsv_path) as tsv:
        header = tsv.readline()  # Skip header
        for line in tsv:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            read_name = parts[0]
            gene = parts[1] if len(parts) > 1 else "NONE"

            if gene == "NONE":
                continue

            # Parse strand information from gene label
            # Format: "GENE" or "GENE (plus)" or "GENE (minus)"
            if " (plus)" in gene:
                base_gene = gene.replace(" (plus)", "")
                strand_info = "plus"
            elif " (minus)" in gene:
                base_gene = gene.replace(" (minus)", "")
                strand_info = "minus"
            else:
                base_gene = gene
                strand_info = "both"

            gene_mapping[read_name] = (base_gene, strand_info)

    # Load FASTQ records into dictionaries for quick lookup.
    # Keys are normalised to match the SAM/TSV read IDs (no /1 or /2).
    r1_records = {}
    for rec in SeqIO.parse(r1_path, "fastq"):
        r1_records[_base_read_id(rec.id)] = rec

    r2_records = {}
    if r2_path:
        for rec in SeqIO.parse(r2_path, "fastq"):
            r2_records[_base_read_id(rec.id)] = rec

    # Organize all records by gene in a single collection
    gene_collection = defaultdict(list)
    is_paired = r2_path is not None

    for read_name, (gene, strand_info) in gene_mapping.items():
        # Determine which mates to include based on strand_info
        # plus = R1 only, minus = R2 only, both = both mates
        include_r1 = strand_info in ("plus", "both")
        include_r2 = strand_info in ("minus", "both")

        # Add R1 record if exists and should be included
        if include_r1:
            r1 = r1_records.get(read_name)
            if r1 is not None:
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

        # Add R2 record if exists and should be included (only for paired-end)
        if is_paired and include_r2:
            r2 = r2_records.get(read_name)
            if r2 is not None:
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

        # Add prefix to gene filename
        filename = f"{prefix}_{gene}.fq" if prefix else f"{gene}.fq"
        output_file = os.path.join(output_dir, filename)

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
    aligner='bwa-mem2',
    threads=4,
    minimap2_profile='short',
    is_paired=False,
    prefix="",
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

        # Use prefix for filenames
        fastq_filename = f"{prefix}_{gene}.fq" if prefix else f"{gene}.fq"
        fastq_file = os.path.join(fastq_dir, fastq_filename)

        bam_filename_base = f"{prefix}_{gene}" if prefix else gene
        tmp_ref = os.path.join(output_dir, f"{bam_filename_base}_ref.fa")
        index_base = os.path.join(output_dir, f"{bam_filename_base}_index")
        sam_path = os.path.join(output_dir, f"{bam_filename_base}.sam")
        bam_path = os.path.join(output_dir, f"{bam_filename_base}.bam")
        sorted_bam = os.path.join(output_dir, f"{bam_filename_base}.sorted.bam")

        # Verify files exist
        if not os.path.exists(fastq_file):
            print(f"Error: FASTQ file {fastq_file} does not exist", file=sys.stderr)
            continue
        # Always write a fresh per-gene reference FASTA (avoid stale content)
        with open(tmp_ref, "w") as f:
            f.write(f">{gene}\n{ref_db[gene]}\n")

        # Build index and align based on aligner choice
        if aligner == 'bowtie2':
            # Build bowtie2 index
            build_cmd = f"bowtie2-build {tmp_ref} {index_base}"
            print(f"Building bowtie2 index: {build_cmd}", file=sys.stderr)
            subprocess.run(build_cmd, shell=True, check=True)

            # Align with bowtie2 using initial stage parameters
            # Note: Always use -U (unpaired) since FASTQ may contain mix of complete pairs and single mates
            align_cmd = f"bowtie2 --local --score-min G,40,40 -p {threads} -x {index_base} -U {fastq_file} -S {sam_path}"
            print(f"Running alignment: {align_cmd}", file=sys.stderr)
            subprocess.run(align_cmd, shell=True, check=True)

        elif aligner == 'bwa-mem2':
            # Build bwa-mem2 index
            build_cmd = f"bwa-mem2 index -p {index_base} {tmp_ref}"
            print(f"Building bwa-mem2 index: {build_cmd}", file=sys.stderr)
            subprocess.run(build_cmd, shell=True, check=True)

            # Align with bwa-mem2
            # Note: Do not use -p flag since FASTQ may contain mix of complete pairs and single mates
            align_cmd = f"bwa-mem2 mem -A 2 -B 8 -T 240 -t {threads} {index_base} {fastq_file} > {sam_path}"
            print(f"Running alignment: {align_cmd}", file=sys.stderr)
            subprocess.run(align_cmd, shell=True, check=True)

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

            # Add stringent score threshold only for short reads
            score_threshold = "-s 240" if preset == "sr" else ""

            # Build minimap2 index
            index_mmi = f"{index_base}.mmi"
            build_cmd = f"minimap2 -x {preset} -d {index_mmi} {tmp_ref}"
            print(f"Building minimap2 ({minimap2_profile}) index: {build_cmd}", file=sys.stderr)
            subprocess.run(build_cmd, shell=True, check=True)

            # Align with minimap2 (paired or single-end)
            align_cmd = f"minimap2 -ax {preset} --MD {score_threshold} -t {threads} {index_mmi} {fastq_file} > {sam_path}"
            print(f"Running alignment: {align_cmd}", file=sys.stderr)
            subprocess.run(align_cmd, shell=True, check=True)

        else:
            print(f"Error: Unsupported aligner '{aligner}'", file=sys.stderr)
            continue

        # Convert to BAM
        subprocess.run(f"samtools view -b {sam_path} > {bam_path}", shell=True, check=True)
        subprocess.run(f"samtools sort -o {sorted_bam} {bam_path}", shell=True, check=True)
        subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)

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
    aligner='bwa-mem2',
    threads=4,
    minimap2_profile='short',
    prefix="",
):
    # Ensure the output directories exist
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(bam_dir, exist_ok=True)

    # Step 1: Process files to create gene-specific FASTQs
    processed_genes = process_files(tsv_file, r1_file, r2_file, fastq_dir, prefix)

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
        prefix=prefix,
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create gene-specific FASTQ and BAM files from read mappings (supports both paired-end and single-end).')
    parser.add_argument('--tsv', required=True, help='Input TSV file mapping read names to genes')
    parser.add_argument('--r1', required=True, help='Forward reads FASTQ file (or single-end reads)')
    parser.add_argument('--r2', required=False, default=None, help='Reverse reads FASTQ file (optional, for paired-end)')
    parser.add_argument('--ref', required=True, dest='ref_fasta', help='Reference sequences FASTA file')
    parser.add_argument('--fastq-dir', default='mapped_fastq', help='Output directory for gene-specific FASTQ files')
    parser.add_argument('--bam-dir', default='bam_files', help='Output directory for BAM files')
    parser.add_argument(
        '--aligner',
        default='bwa-mem2',
        choices=['bwa-mem2', 'bowtie2', 'minimap2'],
        help='Aligner to use (default: BWA-MEM2 for short reads; Bowtie2 is an alternative short-read aligner)'
    )
    parser.add_argument('--threads', type=int, default=4, help='Number of threads for alignment (default: 4)')
    parser.add_argument('--minimap2-profile', default='short',
                        choices=['short', 'pacbio-hifi', 'pacbio-clr', 'ont-q20', 'ont-standard'],
                        help='Minimap2 profile (default: short)')
    parser.add_argument('--prefix', default='', help='Prefix to add to all output filenames (default: none)')

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
        minimap2_profile=args.minimap2_profile,
        prefix=args.prefix
    )
