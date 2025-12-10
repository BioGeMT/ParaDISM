import argparse
from typing import Dict, List
from collections import defaultdict
from Bio import AlignIO
from Bio.Align.sam import AlignmentIterator
from msa_processing import map_seqcoords_to_alncoords


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


def process_read(alignment, msa, seq_to_aln, gene_names):
    """
    Process a single read alignment and check c1/c2 conditions for each gene.

    Returns:
        c1_dict: {gene: bool} - True if read has unique pos for this gene
        c2_dict: {gene: bool} - True if read has no contradictions for this gene
    """
    ref_gene = alignment.target.id
    gene_idx_map = {g: i for i, g in enumerate(gene_names)}

    if ref_gene not in gene_idx_map:
        return {g: False for g in gene_names}, {g: True for g in gene_names}

    ref_idx = gene_idx_map[ref_gene]

    all_msa_cols = set()
    matching = defaultdict(set)  # gene -> set of msa_cols where read matches

    # Get aligned positions from BioPython alignment
    # alignment.aligned returns (target_intervals, query_intervals)
    target_intervals = alignment.aligned[0]
    query_intervals = alignment.aligned[1]
    query_seq = str(alignment.query.seq)

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

            msa_col = seq_to_aln[ref_idx][ref_pos]
            read_base = query_seq[query_pos].upper()

            all_msa_cols.add(msa_col)

            # Check which genes match at this position
            for gene in gene_names:
                gene_idx = gene_idx_map[gene]
                gene_base = str(msa[gene_idx, msa_col]).upper()
                if gene_base != '-' and read_base == gene_base:
                    matching[gene].add(msa_col)

    # Calculate c1 and c2 for each gene
    c1_dict = {}
    c2_dict = {}

    for gene in gene_names:
        # c1: exists a position where read matches this gene only
        c1 = any(
            msa_col in matching[gene] and
            all(msa_col not in matching[g] for g in gene_names if g != gene)
            for msa_col in all_msa_cols
        )

        # c2: for all positions, read matches this gene OR matches no other gene
        c2 = all(
            msa_col in matching[gene] or
            all(msa_col not in matching[g] for g in gene_names if g != gene)
            for msa_col in all_msa_cols
        )

        c1_dict[gene] = c1
        c2_dict[gene] = c2

    return c1_dict, c2_dict

def process_read_simple(alignment, msa, seq_to_aln, gene_names):
    """
    Process a single read alignment and check c1/c2 conditions for each gene.

    Returns:
        c1_dict: {gene: bool} - True if read has unique pos for this gene
        c2_dict: {gene: bool} - True if read has no contradictions for this gene
    """
    ref_gene = alignment.target.id
    gene_idx_map = {g: i for i, g in enumerate(gene_names)}

    if ref_gene not in gene_idx_map:
        raise ValueError('ref_gene not in gene_idx_map')

    ref_idx = gene_idx_map[ref_gene]

    # Calculate c1 and c2 for each gene
    c1_gene_id = -1
    c2_pass = True
    
    for aln_col_id in range(alignment.shape[1]):
        query_pos, target_pos = alignment.indices[:, aln_col_id]
        if query_pos == -1: continue
        if target_pos == -1: continue
        msa_col_id = seq_to_aln[ref_idx][target_pos]
        msa_column = msa[:, msa_col_id]
        query_nt = alignment[0, aln_col_id].upper()
        matching_genes = [gene_id for gene_id, nt in enumerate(msa_column) if nt == query_nt]
        if len(matching_genes) == 1: # Unique matching gene
            if c1_gene_id == -1: # No match seen yet - first match
                c1_gene_id == matching_genes[0]
            elif c1_gene_id != matching_genes[0]: # conflicting matches
                c1_gene_id = -1 # reset - no match
                break  # c1 failed - stop checking

        if c1_gene_id != -1: # c1 passed, check c2 
            for aln_col_id in range(alignment.shape[1]):
                query_pos, target_pos = alignment.indices[:, aln_col_id]
                if query_pos == -1: continue
                if target_pos == -1: continue
                msa_col_id = seq_to_aln[ref_idx][target_pos]
                msa_column = msa[:, msa_col_id]
                # Take query (read) nt and the c1_target nt
                # where c1_target is the gene that passed c1 
                query_nt = alignment[0, aln_col_id].upper()
                target_nt = msa_column[c1_gene_id]
                if query_nt != target_nt: # mismatch
                    if query_nt in msa_column: # match to other gene
                        c2_pass = False
                        break # no need to check further
        if c1_gene_id != -1 and c2_pass:
            return gene_names[c1_gene_id]
        else:
            return None



def process_sam(sam_path, msa, seq_to_aln, gene_names, output_path):
    """
    Process SAM file and assign reads to genes based on c1/c2 conditions.
    """
    # Track c1/c2 per read pair per gene
    qname_to_c1 = {gene: defaultdict(bool) for gene in gene_names}
    qname_to_c2 = {gene: defaultdict(lambda: True) for gene in gene_names}
    all_qnames = set()

    reads_processed = 0
    for alignment in AlignmentIterator(sam_path):
        qname = alignment.query.id
        all_qnames.add(qname)

        c1_dict, c2_dict = process_read(alignment, msa, seq_to_aln, gene_names)

        for gene in gene_names:
            if c1_dict[gene]:
                qname_to_c1[gene][qname] = True
            if not c2_dict[gene]:
                qname_to_c2[gene][qname] = False

        reads_processed += 1

    # Assign reads to genes
    with open(output_path, 'w') as out_f:
        out_f.write('Read_Name\tAssignment\n')

        for qname in sorted(all_qnames):
            passing_genes = []
            for gene in gene_names:
                c1 = qname_to_c1[gene][qname]
                c2 = qname_to_c2[gene][qname]
                if c1 and c2:
                    passing_genes.append(gene)

            if len(passing_genes) == 1:
                assignment = passing_genes[0]
            else:
                assignment = "NONE"

            out_f.write(f'{qname}\t{assignment}\n')


def main():
    parser = argparse.ArgumentParser(
        description='ParaDISM: Paralog-aware read assignment using SNPs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sam', required=True,
                        help='Input SAM/BAM file')
    parser.add_argument('--msa', required=True,
                        help='Input MSA FASTA file')
    parser.add_argument('--output', required=True,
                        help='Output file for unique mapping results')
    parser.add_argument('--single-end', action='store_true',
                        help='Treat input as single-end reads')

    args = parser.parse_args()

    msa, seq_to_aln, gene_names = load_msa(args.msa)
    all_chars = set(''.join(str(alnseqrec.seq) for alnseqrec in msa))
    assert all(char.isupper() or char == '-' for char in all_chars), 'MSA needs to be uppercase'

    process_sam(args.sam, msa, seq_to_aln, gene_names, args.output)


if __name__ == '__main__':
    main()
