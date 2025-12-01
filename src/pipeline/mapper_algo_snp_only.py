import argparse
import gc
import sys
from typing import Dict, List, Tuple
from collections import defaultdict

def load_msa_mapping(msa_filepath: str) -> Tuple[Dict, Dict, List]:
    """Load MSA mapping data, excluding gap handling since we only process SNPs"""
    sequence_positions = {}
    sequence_bases = {}
    gene_names = []

    with open(msa_filepath, 'r') as f:
        header = next(f).strip().split('\t')
        
        for col in header:
            if col.endswith('_Position') and col != 'MSA_Position':
                gene_name = col.replace('_Position', '')
                gene_names.append(gene_name)
                sequence_positions[gene_name] = {}

        for line in f:
            fields = line.strip().split('\t')
            msa_pos = int(fields[0])
            sequence_bases[msa_pos] = {}

            for gene in gene_names:
                pos_idx = header.index(f'{gene}_Position')
                base_idx = header.index(f'{gene}_Base')

                pos = fields[pos_idx]
                base = fields[base_idx]
                sequence_bases[msa_pos][gene] = base

                if pos != '-':
                    sequence_positions[gene][int(pos)] = msa_pos

    return sequence_positions, sequence_bases, gene_names

def process_read_data(read_name: str, read_rows: List[List], sequence_bases: Dict, gene_names: List) -> List[Dict]:
    """Process read data for SNP-only analysis"""
    processed_rows = []
    
    for scenario in read_rows:
        scenario_rows = []
        for row in scenario:
            read_name, read_pos, msa_pos, read_base = row
            row_dict = {
                'Read_Name': read_name,
                'Read_Position': read_pos,
                'MSA_Position': msa_pos,
                'Read_Base': read_base
            }
            for gene in gene_names:
                row_dict[f'{gene}_Base'] = sequence_bases[msa_pos][gene]
            
            scenario_rows.append(row_dict)
        processed_rows.append(scenario_rows)
    
    return processed_rows

def find_unique_mapping_snp_only(read_scenarios: Dict[str, List], genes: List[str]) -> Tuple[str, bool, str]:
    """Find unique mapping for SNP-only analysis; returns assignment, strand coverage, and which strand(s) mapped."""
    plus_data = read_scenarios.get('plus', [])
    minus_data = read_scenarios.get('minus', [])

    # Check if we have no data at all
    if not plus_data and not minus_data:
        return "NONE", False, ""

    # Determine if we have both strands (paired-end) or just one (single-end)
    has_both_strands = bool(plus_data and minus_data)

    def check_conditions(scenarios: List, gene: str) -> Tuple[bool, bool]:
        if not scenarios:
            return False, False

        all_positions = set()
        matching_positions = defaultdict(set)

        for scenario in scenarios:
            for row in scenario:
                msa_pos = row['MSA_Position']
                all_positions.add(msa_pos)
                for g in genes:
                    if row['Read_Base'].upper() == row[f'{g}_Base'].upper():
                        matching_positions[g].add(msa_pos)

        # Check condition 1: unique match
        c1_pass = any(
            pos in matching_positions[gene] and
            all(pos not in matching_positions[g] for g in genes if g != gene)
            for pos in all_positions
        )

        # Check condition 2: no contradictions
        c2_pass = all(
            pos in matching_positions[gene] or
            all(pos not in matching_positions[g] for g in genes if g != gene)
            for pos in all_positions
        )

        return c1_pass, c2_pass

    mapped_genes = []

    if has_both_strands:
        # Paired-end logic: require one strand with c1+c2, other strand with c2
        for gene in genes:
            plus_c1, plus_c2 = check_conditions(plus_data, gene)
            minus_c1, minus_c2 = check_conditions(minus_data, gene)

            if ((plus_c1 and plus_c2 and minus_c2) or
                (minus_c1 and minus_c2 and plus_c2)):
                mapped_genes.append(gene)
    else:
        # Single-end logic: require the available strand to pass both c1 and c2
        strand_data = plus_data if plus_data else minus_data
        for gene in genes:
            c1, c2 = check_conditions(strand_data, gene)
            if c1 and c2:
                mapped_genes.append(gene)

    assignment = mapped_genes[0] if len(mapped_genes) == 1 else "NONE"

    # Determine which strand(s) mapped
    strand_info = ""
    if assignment != "NONE":
        if has_both_strands:
            strand_info = "both"
        elif plus_data:
            strand_info = "plus"
        elif minus_data:
            strand_info = "minus"

    return assignment, has_both_strands, strand_info

def process_read_mappings_snp_only(read_map_filepath: str, sequence_positions: Dict,
                                   sequence_bases: Dict, gene_names: List, output_file: str,
                                   batch_size: int = 10000, paired_mode: bool = True) -> None:
    """Process read mappings for SNP-only analysis (no insertions)"""

    # First pass: count total unique read pairs
    total_reads = 0
    seen_bases = set()
    with open(read_map_filepath, 'r') as f:
        next(f)  # Skip header
        for line in f:
            read_name = line.split('\t')[0].rstrip('+-')
            if read_name not in seen_bases:
                seen_bases.add(read_name)
                total_reads += 1

    with open(output_file, 'w') as out_f:
        out_f.write('Read_Name\tUniquely_Mapped\n')

        current_read_base = None
        reads_processed = 0
        
        match_buffer = []
        current_read = None
        current_ref_gene = None
        read_scenarios = defaultdict(lambda: {'plus': [], 'minus': []})
        
        def process_buffer():
            nonlocal match_buffer, current_read, current_ref_gene, read_scenarios
            
            if not current_read or not match_buffer:
                return

            # For SNP-only, we only have one scenario per read (no insertion mapping)
            base_read_name = current_read.rstrip('+-')
            strand = 'plus' if current_read.endswith('+') else 'minus'
            
            processed = process_read_data(
                current_read,
                [match_buffer],  # Single scenario
                sequence_bases,
                gene_names
            )
            read_scenarios[base_read_name][strand].extend(processed)

            match_buffer.clear()

        def process_current_read():
            nonlocal current_read_base, reads_processed, read_scenarios
            
            if current_read_base and read_scenarios:
                unique_gene, has_both, strand_info = find_unique_mapping_snp_only(
                    read_scenarios[current_read_base],
                    gene_names
                )
                label = unique_gene
                if paired_mode and unique_gene != "NONE" and not has_both:
                    # Indicate which specific strand mapped
                    label = f"{unique_gene} ({strand_info})"
                out_f.write(f'{current_read_base}\t{label}\n')
                del read_scenarios[current_read_base]
                reads_processed += 1

                if reads_processed % 100 == 0:
                    sys.stdout.write(f"\rPaired reads processed: {reads_processed:,}/{total_reads:,}")
                    sys.stdout.flush()

                if reads_processed % batch_size == 0:
                    gc.collect()

        with open(read_map_filepath, 'r') as input_file:
            next(input_file)  # Skip header

            for line in input_file:
                fields = line.strip().split('\t')
                read_name = fields[0]
                read_pos = fields[1]
                ref_pos = fields[2]
                read_base = fields[3]
                alignment_type = fields[5]
                ref_gene = fields[6]

                read_base_name = read_name.rstrip('+-')
                if read_base_name != current_read_base:
                    process_buffer()
                    process_current_read()
                    current_read_base = read_base_name

                if read_name != current_read:
                    if current_read is not None:
                        process_buffer()
                    current_read = read_name
                    current_ref_gene = ref_gene

                # Only process matches and mismatches (skip insertions and deletions)
                if alignment_type in ['match', 'mismatch']:
                    if ref_pos != '-':
                        ref_pos_int = int(ref_pos)
                        if ref_pos_int in sequence_positions[ref_gene]:
                            msa_pos = sequence_positions[ref_gene][ref_pos_int]
                            row_data = (read_name, read_pos, msa_pos, read_base)
                            match_buffer.append(row_data)

            # Process final read
            if current_read is not None:
                process_buffer()
                process_current_read()

            sys.stdout.write(f"\rPaired reads processed: {reads_processed:,}/{total_reads:,}\n")
            sys.stdout.flush()

def main():
    parser = argparse.ArgumentParser(
        description='SNP-only mapper: Find unique mappings based on matches/mismatches only.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--read_map', required=True,
                      help='Input read mapping TSV file')
    parser.add_argument('--msa', required=True,
                      help='Input MSA mapping TSV file')
    parser.add_argument('--output_file', required=True,
                      help='Output file for unique mapping results')
    parser.add_argument('--batch_size', type=int, default=10000,
                      help='Number of reads to process before memory cleanup')
    parser.add_argument('--single-end', action='store_true',
                      help='Treat input as single-end (disables mate annotations)')

    args = parser.parse_args()

    sequence_positions, sequence_bases, gene_names = load_msa_mapping(args.msa)
    
    process_read_mappings_snp_only(
        args.read_map, sequence_positions, sequence_bases,
        gene_names, args.output_file, args.batch_size,
        paired_mode=not args.single_end
    )

if __name__ == '__main__':
    import time
    start_time = time.perf_counter()

    main()

    end_time = time.perf_counter()
    print(f"SNP-only mapper execution time: {end_time - start_time:.6f} seconds", file=sys.stderr)
