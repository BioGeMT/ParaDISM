import argparse
import os
import subprocess
from subprocess import DEVNULL
from collections import defaultdict
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.sam import AlignmentIterator
from .msa_processing import map_seqcoords_to_alncoords


def load_msa(msa_fasta_path: str):
    """
    Returns:
        msa: MultipleSeqAlignment object (for base lookup)
        seq_to_aln: list of lists - seq_to_aln[gene_idx][seq_pos] = msa_col (0-based)
        gene_names: list of gene names
    """
    msa = AlignIO.read(msa_fasta_path, 'fasta')
    gene_names = [record.id for record in msa]
    seq_to_aln = map_seqcoords_to_alncoords(msa)

    return msa, seq_to_aln, gene_names


def process_read_simple(alignment, msa, seq_to_aln, gene_names, gene_idx_map):
    """
    Process a single read alignment and check c1/c2 conditions for each gene.
    Optimized version that uses alignment.aligned to skip gaps and processes c1/c2 in single pass.

    Returns: str
    The name of the mapped gene or NONE if no gene could be uniquely assigned.  
    """
    ref_gene = alignment.target.id
    try:
        ref_idx = gene_idx_map[ref_gene]
    except KeyError:
        raise KeyError('Reference gene name from read alignment not in the gene names list from MSA')
    
    # Calculate c1 and c2 for each gene
    c1_gene_id = -1
    c2_pass = True
    
    # Get aligned positions from BioPython alignment (skips gaps efficiently)
    target_intervals = alignment.aligned[0]
    query_intervals = alignment.aligned[1]
    query_seq = str(alignment.query.seq)
    
    # Single pass: check c1 and c2 together
    for target_interval, query_interval in zip(target_intervals, query_intervals):
        t_start, t_end = target_interval
        q_start, q_end = query_interval
        
        # Only process segments where both sequences advance equally (SNPs/matches, skip indels)
        if (t_end - t_start) != (q_end - q_start):
            continue
        
        for offset in range(t_end - t_start):
            ref_pos = t_start + offset  # 0-based ref position
            query_pos = q_start + offset
            
            if ref_pos >= len(seq_to_aln[ref_idx]):
                continue
            
            msa_col_id = seq_to_aln[ref_idx][ref_pos]
            read_base = query_seq[query_pos].upper()
            
            # Skip ambiguous nucleotides
            if read_base == 'N':
                continue
            
            # Check which genes match at this position
            matching_genes = []
            for gene_idx, gene_name in enumerate(gene_names):
                gene_base = str(msa[gene_idx, msa_col_id]).upper()
                if gene_base != '-' and read_base == gene_base:
                    matching_genes.append(gene_idx)
            
            if len(matching_genes) == 1:  # Unique matching gene
                if c1_gene_id == -1:  # No match seen yet - first match
                    c1_gene_id = matching_genes[0]
                elif c1_gene_id != matching_genes[0]:  # conflicting matches
                    c1_gene_id = -1  # reset - no match
                    break  # c1 failed - stop checking
            
            # Check c2 if we have a candidate target gene
            if c1_gene_id != -1:
                target_base = str(msa[c1_gene_id, msa_col_id]).upper()
                if read_base != target_base and matching_genes:
                    c2_pass = False
                    break  # no need to check further
        
        # Early termination checks
        if c1_gene_id == -1:
            break  # Already failed c1, no need to continue
        if not c2_pass:
            break  # Already failed c2, no need to continue
    
    if c1_gene_id != -1 and c2_pass:
        return gene_names[c1_gene_id]
    else:
        return "NONE"


def process_sam_to_dict_simple(sam_path, msa, seq_to_aln, gene_names):
    """
    Process SAM file and assign reads to genes based on c1/c2 conditions.
    
    Returns:
        dict: {read_name: gene_assignment} where gene_assignment is gene name or "NONE"
    """
    gene_idx_map = {g: i for i, g in enumerate(gene_names)}
    assignments = {}
    for alignment in AlignmentIterator(sam_path):
        qname = alignment.query.id
        assigned_gene = process_read_simple(alignment,
                                            msa,
                                            seq_to_aln,
                                            gene_names,
                                            gene_idx_map)
        if qname in assignments:
            # we've already seen the mate of this read
            if assignments[qname] == "NONE":
                # one mate could not be uniquely mapped,
                # but this one can - we go with this one
                assignments[qname] = assigned_gene
            elif assignments[qname] != assigned_gene:
                # two mates have conflicting assignments,
                # we don't assign
                assignments[qname] = "NONE"
            # else: two mates have identical assignments,
            # we don't need to do anything
        else:
            assignments[qname] = assigned_gene
    return assignments


def process_sam_to_dict_simple_parallel(sam_path, msa, seq_to_aln, gene_names, n_jobs):
    """
    Process SAM file and assign reads to genes based on c1/c2 conditions.
    
    Returns:
        dict: {read_name: gene_assignment} where gene_assignment is gene name or "NONE"
    """
    from joblib import Parallel, delayed # we import here in case someone doesn't have joblib
    def worker_generator(sam_path, msa, seq_to_aln, gene_names):
        for alignment in AlignmentIterator(sam_path):
            qname = alignment.query.id
            # delayed() may be unnecessary here, try to run without it: 
            assigned_gene = delayed(process_read_simple)(alignment,
                                            msa,
                                            seq_to_aln,
                                            gene_names)
            yield (qname, assigned_gene)
            
    assignments = {}
    workers = worker_generator(sam_path, msa, seq_to_aln, gene_names)
    parallel = Parallel(n_jobs = n_jobs, return_as = 'generator_unordered')
    for qname, assignment in parallel(workers):
        if qname in assignments:
            # we've already seen the mate of this read
            if assignments[qname] == "NONE":
                # one mate could not be uniquely mapped,
                # but this one can - we go with this one
                assignments[qname] = assigned_gene
            elif assignments[qname] != assigned_gene:
                # two mates have conflicting assignments,
                # we don't assign
                assignments[qname] = "NONE"
            # else: two mates have identical assignments,
            # we don't need to do anything
        else:
            assignments[qname] = assigned_gene
    return assignments



def _base_read_id(read_id: str) -> str:
    """Normalize read ID to match SAM read IDs."""
    if read_id.endswith("/1") or read_id.endswith("/2"):
        return read_id[:-2]
    return read_id


def write_fastq_outputs(assignments: dict[str, str], r1_path: str, r2_path: str | None, 
                       fastq_dir: str, prefix: str = "") -> list[str]:
    """Write per-gene FASTQ files directly from assignments dict."""
    os.makedirs(fastq_dir, exist_ok=True)
    
    # Load FASTQ records
    r1_records = {}
    for rec in SeqIO.parse(r1_path, "fastq"):
        r1_records[_base_read_id(rec.id)] = rec

    r2_records = {}
    if r2_path:
        for rec in SeqIO.parse(r2_path, "fastq"):
            r2_records[_base_read_id(rec.id)] = rec

    # Organize reads by gene
    gene_collection = defaultdict(list)
    is_paired = r2_path is not None

    for read_name, gene in assignments.items():
        if gene == "NONE":
            continue

        # Add R1
        r1 = r1_records.get(read_name)
        if r1 is not None:
            if is_paired:
                gene_collection[gene].append(SeqRecord(
                    r1.seq,
                    id=f"{read_name}/1",
                    description="",
                    letter_annotations=r1.letter_annotations
                ))
            else:
                gene_collection[gene].append(SeqRecord(
                    r1.seq,
                    id=read_name,
                    description="",
                    letter_annotations=r1.letter_annotations
                ))

        # Add R2 if paired-end
        if is_paired:
            r2 = r2_records.get(read_name)
            if r2 is not None:
                gene_collection[gene].append(SeqRecord(
                    r2.seq,
                    id=f"{read_name}/2",
                    description="",
                    letter_annotations=r2.letter_annotations
                ))

    # Write FASTQ files
    processed_genes = []
    for gene, records in gene_collection.items():
        if not records:
            continue

        filename = f"{prefix}_{gene}.fq" if prefix else f"{gene}.fq"
        output_file = os.path.join(fastq_dir, filename)

        with open(output_file, "w") as out_f:
            SeqIO.write(records, out_f, "fastq")

        processed_genes.append(gene)

    return processed_genes


def create_bam_files(genes: list[str], ref_fasta: str, fastq_dir: str, output_dir: str,
                     aligner: str = 'bwa-mem2', threads: int = 4, minimap2_profile: str = 'short',
                     prefix: str = "") -> None:
    """Create BAM files for each gene."""
    ref_db = {}
    for record in SeqIO.parse(ref_fasta, "fasta"):
        gene_name = record.id.split()[0]
        ref_db[gene_name] = str(record.seq)

    os.makedirs(output_dir, exist_ok=True)

    for gene in genes:
        if gene not in ref_db:
            continue

        fastq_filename = f"{prefix}_{gene}.fq" if prefix else f"{gene}.fq"
        fastq_file = os.path.join(fastq_dir, fastq_filename)

        if not os.path.exists(fastq_file):
            continue

        bam_filename_base = f"{prefix}_{gene}" if prefix else gene
        tmp_ref = os.path.join(output_dir, f"{bam_filename_base}_ref.fa")
        index_base = os.path.join(output_dir, f"{bam_filename_base}_index")
        sam_path = os.path.join(output_dir, f"{bam_filename_base}.sam")
        bam_path = os.path.join(output_dir, f"{bam_filename_base}.bam")
        sorted_bam = os.path.join(output_dir, f"{bam_filename_base}.sorted.bam")

        # Write per-gene reference
        with open(tmp_ref, "w") as f:
            f.write(f">{gene}\n{ref_db[gene]}\n")

        # Build index and align (suppress verbose output)
        if aligner == 'bowtie2':
            subprocess.run(f"bowtie2-build {tmp_ref} {index_base}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
            subprocess.run(f"bowtie2 --local --score-min G,40,40 -p {threads} -x {index_base} -U {fastq_file} -S {sam_path}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
        elif aligner == 'bwa-mem2':
            subprocess.run(f"bwa-mem2 index -p {index_base} {tmp_ref}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
            bwa_min_score = 240
            awk_filter = f"awk '/^@/{{print;next}} $3==\"*\"{{print;next}} {{for(i=12;i<=NF;i++)if($i~/^AS:i:/){{split($i,a,\":\");if(a[3]>={bwa_min_score})print;next}}}}'"
            subprocess.run(f"bwa-mem2 mem -A 2 -B 8 -T {bwa_min_score} -t {threads} {index_base} {fastq_file} | {awk_filter} > {sam_path}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
        elif aligner == 'minimap2':
            preset_map = {
                'short': 'sr',
                'pacbio-hifi': 'map-hifi',
                'pacbio-clr': 'map-pb',
                'ont-q20': 'lr:hq',
                'ont-standard': 'map-ont',
            }
            preset = preset_map.get(minimap2_profile, 'sr')
            score_threshold = "-s 240" if preset == "sr" else ""
            index_mmi = f"{index_base}.mmi"
            subprocess.run(f"minimap2 -x {preset} -d {index_mmi} {tmp_ref}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
            subprocess.run(f"minimap2 -ax {preset} --MD {score_threshold} -t {threads} {index_mmi} {fastq_file} > {sam_path}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)

        # Convert to BAM (suppress verbose output)
        subprocess.run(f"samtools view -b {sam_path} > {bam_path}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
        subprocess.run(f"samtools sort -o {sorted_bam} {bam_path}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)
        subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True, stdout=DEVNULL, stderr=DEVNULL)

        # Cleanup
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


def main():
    parser = argparse.ArgumentParser(
        description='ParaDISM: Paralog-aware read assignment using SNPs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sam', required=True,
                        help='Input SAM/BAM file')
    parser.add_argument('--msa', required=True,
                        help='Input MSA FASTA file')
    parser.add_argument('--r1', required=True,
                        help='R1 FASTQ file (or single-end reads)')
    parser.add_argument('--r2', required=False, default=None,
                        help='R2 FASTQ file (optional, for paired-end)')
    parser.add_argument('--ref', required=True,
                        help='Reference FASTA file')
    parser.add_argument('--fastq-dir', required=True,
                        help='Output directory for gene-specific FASTQ files')
    parser.add_argument('--bam-dir', required=True,
                        help='Output directory for BAM files')
    parser.add_argument('--aligner', default='bwa-mem2',
                        choices=['bwa-mem2', 'bowtie2', 'minimap2'],
                        help='Aligner to use')
    parser.add_argument('--threads', type=int, default=4,
                        help='Number of threads')
    parser.add_argument('--minimap2-profile', default='short',
                        choices=['short', 'pacbio-hifi', 'pacbio-clr', 'ont-q20', 'ont-standard'],
                        help='Minimap2 profile')
    parser.add_argument('--prefix', default='',
                        help='Prefix for output files')

    args = parser.parse_args()

    msa, seq_to_aln, gene_names = load_msa(args.msa)
    all_chars = set(''.join(str(alnseqrec.seq) for alnseqrec in msa))
    assert all(char.isupper() or char == '-' for char in all_chars), 'MSA needs to be uppercase'

    assignments = process_sam_to_dict_simple(args.sam, msa, seq_to_aln, gene_names)
    
    # Write FASTQ files
    genes = write_fastq_outputs(
        assignments,
        args.r1,
        args.r2,
        args.fastq_dir,
        args.prefix
    )
    
    # Create BAM files
    if genes:
        create_bam_files(
            genes,
            args.ref,
            args.fastq_dir,
            args.bam_dir,
            args.aligner,
            args.threads,
            args.minimap2_profile,
            args.prefix
        )



if __name__ == '__main__':
    main()
