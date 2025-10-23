import argparse
import gc
import sys
from typing import Dict, List, Tuple
from collections import defaultdict

def load_msa_mapping(msa_filepath: str) -> Tuple[Dict, Dict, Dict, List]:
    sequence_positions = {}
    sequence_bases = {}
    gap_bases_dict = {}
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

            bases_at_position = set()
            for gene in gene_names:
                pos_idx = header.index(f'{gene}_Position')
                base_idx = header.index(f'{gene}_Base')

                pos = fields[pos_idx]
                base = fields[base_idx]
                sequence_bases[msa_pos][gene] = base

                if pos != '-':
                    sequence_positions[gene][int(pos)] = msa_pos
                    if base != '-':
                        bases_at_position.add(base.upper())

            if any(sequence_bases[msa_pos][gene] == '-' for gene in gene_names):
                gap_bases_dict[msa_pos] = bases_at_position

    return sequence_positions, sequence_bases, gap_bases_dict, gene_names

def map_insertions_dp(insertions: List[Tuple], valid_gaps: List[Tuple], max_scenarios: int = 100) -> List[List]:
    memo = {}

    def dp(i: int, j: int) -> Tuple[int, List]:
        if i == len(insertions) or j == len(valid_gaps):
            return (0, [[]])

        if (i, j) in memo:
            return memo[(i, j)]

        count1, mappings1 = dp(i+1, j)
        count2, mappings2 = dp(i, j+1)
        
        insertion_base = insertions[i][1].upper()
        gap_valid_bases = valid_gaps[j][1]
        if insertion_base in gap_valid_bases:
            count3, mappings3 = dp(i+1, j+1)
            count3 += 1
            new_mappings3 = [[(insertions[i][0], valid_gaps[j][0])] + m for m in mappings3]
        else:
            count3, new_mappings3 = (0, [])

        max_count = max(count1, count2, count3)
        mappings = []
        if count1 == max_count:
            mappings.extend(mappings1)
        if count2 == max_count:
            mappings.extend(mappings2)
        if count3 == max_count:
            mappings.extend(new_mappings3)

        if len(mappings) > max_scenarios:
            mappings = mappings[:max_scenarios]

        memo[(i, j)] = (max_count, mappings)
        return memo[(i, j)]

    max_mapped, all_mappings = dp(0, 0)
    
    unique_mappings = []
    seen = set()
    for m in all_mappings[:max_scenarios]:
        m_sorted = tuple(sorted(m))
        if m_sorted not in seen:
            seen.add(m_sorted)
            unique_mappings.append(m)
            if len(unique_mappings) >= max_scenarios:
                break

    memo.clear()
    return unique_mappings

def map_insertions(insertions: List[Tuple], start_msa_pos: int, end_msa_pos: int, 
                  sequence_bases: Dict, gap_bases_dict: Dict, ref_gene: str, 
                  max_scenarios: int = 100) -> List[List]:
    valid_gaps = []
    for pos in range(start_msa_pos + 1, end_msa_pos):
        if (pos in sequence_bases and
            sequence_bases[pos][ref_gene] == '-' and
            pos in gap_bases_dict):
            valid_gaps.append((pos, gap_bases_dict[pos]))

    if not valid_gaps or not insertions:
        return [[]]

    return map_insertions_dp(insertions, valid_gaps, max_scenarios)

def process_read_scenarios(read_name: str, scenario_rows: List[List], sequence_bases: Dict, gene_names: List) -> List[Dict]:
    processed_scenarios = []
    
    for scenario in scenario_rows:
        seen = set()
        unique_rows = []
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
            
            key = (read_name, read_pos, msa_pos, read_base)
            if key not in seen:
                seen.add(key)
                unique_rows.append(row_dict)
        
        processed_scenarios.append(unique_rows)
    
    return processed_scenarios

def find_unique_mapping_in_scenarios(read_scenarios: Dict[str, List], genes: List[str]) -> str:
    """Decide unique mapping; supports paired-end and single-end."""
    plus_data = read_scenarios.get('plus', [])
    minus_data = read_scenarios.get('minus', [])

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
        # Paired-end logic: one strand must pass c1+c2, the other c2
        for gene in genes:
            plus_c1, plus_c2 = check_conditions(plus_data, gene)
            minus_c1, minus_c2 = check_conditions(minus_data, gene)

            if ((plus_c1 and plus_c2 and minus_c2) or 
                (minus_c1 and minus_c2 and plus_c2)):
                mapped_genes.append(gene)
    else:
        # Single-end logic: the available strand must pass both c1 and c2
        strand_data = plus_data if plus_data else minus_data
        for gene in genes:
            c1, c2 = check_conditions(strand_data, gene)
            if c1 and c2:
                mapped_genes.append(gene)
    
    return mapped_genes[0] if len(mapped_genes) == 1 else "NONE"

def process_read_mappings(read_map_filepath: str, sequence_positions: Dict,
                         sequence_bases: Dict, gap_bases_dict: Dict,
                         gene_names: List, output_file: str,
                         batch_size: int = 10000, max_scenarios: int = 100) -> None:
    
    with open(output_file, 'w') as out_f:
        out_f.write('Read_Name\tUniquely_Mapped\n')
        
        current_read_base = None
        reads_processed = 0
        
        match_buffer = []
        insertion_buffer = []
        post_match_buffer = []
        last_ref_pos = None
        current_read = None
        current_ref_gene = None
        read_scenarios = defaultdict(lambda: {'plus': [], 'minus': []})
        
        def process_buffers():
            nonlocal match_buffer, insertion_buffer, post_match_buffer, last_ref_pos
            nonlocal current_read, current_ref_gene, read_scenarios
            
            if not current_read:
                return

            base_scenario = match_buffer.copy()
            scenarios_to_process = []

            if insertion_buffer:
                if last_ref_pos is not None and current_ref_gene is not None:
                    ref_pos_int = int(last_ref_pos)
                    current_msa_pos = sequence_positions[current_ref_gene][ref_pos_int]
                    all_positions = sorted(sequence_positions[current_ref_gene].keys())
                    next_ref_pos_candidates = [p for p in all_positions if p > ref_pos_int]
                    
                    if next_ref_pos_candidates:
                        next_ref_pos = min(next_ref_pos_candidates)
                        next_msa_pos = sequence_positions[current_ref_gene][next_ref_pos]
                        all_mapped_mappings = map_insertions(
                            insertion_buffer, current_msa_pos, next_msa_pos,
                            sequence_bases, gap_bases_dict, current_ref_gene,
                            max_scenarios
                        )
                        
                        for mapping in all_mapped_mappings:
                            scenario = base_scenario.copy()
                            for (read_pos, msa_pos) in mapping:
                                read_base = next((base for (pos, base) in insertion_buffer if pos == read_pos), '-')
                                scenario.append((current_read, read_pos, msa_pos, read_base))
                            scenario.extend(post_match_buffer)
                            scenarios_to_process.append(scenario)
            else:
                base_scenario.extend(post_match_buffer)
                scenarios_to_process.append(base_scenario)

            base_read_name = current_read.rstrip('+-')
            strand = 'plus' if current_read.endswith('+') else 'minus'
            
            processed = process_read_scenarios(
                current_read,
                scenarios_to_process[:max_scenarios],
                sequence_bases,
                gene_names
            )
            read_scenarios[base_read_name][strand].extend(processed)

            match_buffer.clear()
            insertion_buffer.clear()
            post_match_buffer.clear()
            last_ref_pos = None

        def process_current_read():
            nonlocal current_read_base, reads_processed, read_scenarios
            
            if current_read_base and read_scenarios:
                unique_gene = find_unique_mapping_in_scenarios(
                    read_scenarios[current_read_base],
                    gene_names
                )
                out_f.write(f'{current_read_base}\t{unique_gene}\n')
                del read_scenarios[current_read_base]
                reads_processed += 1
                
                if reads_processed % 100 == 0:
                    sys.stdout.write(f"\rReads processed: {reads_processed}")
                    sys.stdout.flush()
                
                if reads_processed % batch_size == 0:
                    gc.collect()

        with open(read_map_filepath, 'r') as input_file:
            next(input_file)
            
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
                    process_buffers()
                    process_current_read()
                    current_read_base = read_base_name

                if read_name != current_read:
                    if current_read is not None:
                        process_buffers()
                    current_read = read_name
                    current_ref_gene = ref_gene

                if alignment_type == 'deletion':
                    continue
                elif alignment_type in ['match', 'mismatch']:
                    if ref_pos != '-':
                        ref_pos_int = int(ref_pos)
                        if ref_pos_int in sequence_positions[ref_gene]:
                            msa_pos = sequence_positions[ref_gene][ref_pos_int]
                            row_data = (read_name, read_pos, msa_pos, read_base)

                            if not insertion_buffer:
                                match_buffer.append(row_data)
                                last_ref_pos = ref_pos
                            else:
                                post_match_buffer.append(row_data)

                elif alignment_type == 'insertion':
                    insertion_buffer.append((read_pos, read_base))

            # Process final read
            if current_read is not None:
                process_buffers()
                process_current_read()
                
            sys.stdout.write(f"\rReads processed: {reads_processed}\n")
            sys.stdout.flush()

def main():
    parser = argparse.ArgumentParser(
        description='Map read insertions to MSA gaps and find unique mappings.',
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
    parser.add_argument('--max_scenarios', type=int, default=100,
                      help='Maximum number of mapping scenarios per read')

    args = parser.parse_args()

    sequence_positions, sequence_bases, gap_bases_dict, gene_names = load_msa_mapping(args.msa)
    
    process_read_mappings(
        args.read_map, sequence_positions, sequence_bases,
        gap_bases_dict, gene_names, args.output_file,
        args.batch_size, args.max_scenarios
    )

if __name__ == '__main__':
    import time
    start_time = time.perf_counter()
    
    main()
    
    end_time = time.perf_counter()
    print(f"Pipeline execution time: {end_time - start_time:.6f} seconds")
