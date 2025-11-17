#!/usr/bin/env python3
"""Generate mutated PKD1/PKD1P references plus simulated read pairs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import numpy.random as rd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.stats import geom


NUCL = ["A", "C", "T", "G"]
GENES = ["PKD1"] + [f"PKD1P{i}" for i in range(1, 7)]

# Genome-level parameters (match notebook defaults)
SNP_PROB = 0.005
INDEL_PROB = 0.0005
INDEL_LEN_MEAN = 2

# Sequencing defaults
DEFAULT_ERROR = [0.01]
READ_LENGTH = 150
FRAGMENT_LENGTH_MEAN = 350
DEFAULT_READS = 100_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Simulate mutated paralog references and paired-end reads.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--reference", default="ref.fa", help="Input reference FASTA")
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for all generated files",
    )
    parser.add_argument(
        "--num-reads",
        type=int,
        default=DEFAULT_READS,
        help="Number of read pairs to simulate per error rate",
    )
    parser.add_argument(
        "--error-rate",
        type=float,
        action="append",
        help="Sequencing error rate (Phred substitution). Repeat to supply multiple values.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for numpy/scipy RNG",
    )
    parser.add_argument(
        "--snp-prob",
        type=float,
        default=SNP_PROB,
        help="Probability of a SNP at each reference base",
    )
    parser.add_argument(
        "--indel-prob",
        type=float,
        default=INDEL_PROB,
        help="Probability of an indel event at each reference base",
    )
    return parser.parse_args()


def seed_rng(seed: int) -> None:
    rd.seed(seed)
    np.random.seed(seed)


def load_references(path: Path) -> dict[str, SeqRecord]:
    records = {record.id: record for record in SeqIO.parse(str(path), "fasta")}
    missing = [gene for gene in GENES if gene not in records]
    if missing:
        raise ValueError(f"Reference FASTA missing genes: {', '.join(missing)}")
    return records


def mutate_sequences(
    refseq_dict: dict[str, SeqRecord],
    snp_prob: float,
    indel_prob: float,
) -> tuple[list[SeqRecord], list[tuple[str, int, str, str, str]]]:
    mutation_log: list[tuple[str, int, str, str, str]] = []
    mutated_records: list[SeqRecord] = []

    for gene in GENES:
        ref_seq = str(refseq_dict[gene].seq)
        mut_seq = list(ref_seq)

        # SNPs
        for pos, ref_base in enumerate(ref_seq):
            if rd.random() < snp_prob:
                alt_bases = [b for b in NUCL if b != ref_base]
                alt_base = rd.choice(alt_bases)
                mutation_log.append((gene, pos, "SNP", ref_base, alt_base))
                mut_seq[pos] = alt_base

        # Indels
        pos = 0
        while pos < len(ref_seq):
            if rd.random() < indel_prob:
                is_insertion = bool(rd.choice([0, 1]))
                indel_len = max(1, int(geom.rvs(p=1 / INDEL_LEN_MEAN)))

                if is_insertion:
                    ins_seq = "".join(rd.choice(NUCL, size=indel_len))
                    mutation_log.append((gene, pos, "INS", ref_seq[pos], ins_seq))

                    for idx, base in enumerate(ins_seq):
                        mut_seq.insert(pos + 1 + idx, base)
                    pos += indel_len
                else:
                    indel_len = min(indel_len, len(ref_seq) - pos)
                    if indel_len > 0:
                        del_seq = ref_seq[pos : pos + indel_len]
                        mutation_log.append((gene, pos, "DEL", del_seq, ""))
                        del mut_seq[pos : pos + indel_len]
            pos += 1

        mut_record = SeqRecord(Seq("".join(mut_seq)), id=f"{gene}_mutated", description="")
        mutated_records.append(mut_record)

    return mutated_records, mutation_log


def write_reference_outputs(
    refseq_dict: dict[str, SeqRecord],
    mutated_records: list[SeqRecord],
    mutation_log: list[tuple[str, int, str, str, str]],
    output_dir: Path,
) -> None:
    ref_lengths = {gene: len(refseq_dict[gene]) for gene in GENES}

    fasta_path = output_dir / "all_genes_sequences.fasta"
    with fasta_path.open("w") as handle:
        for idx, gene in enumerate(GENES):
            handle.write(f">{gene}_reference\n{refseq_dict[gene].seq}\n")
            handle.write(f">{gene}_mutated\n{mutated_records[idx].seq}\n")

    sorted_mutations = sorted(mutation_log, key=lambda entry: (entry[0], entry[1]))

    vcf_path = output_dir / "all_genes_mutations.vcf"
    with vcf_path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##fileDate=20250315\n")
        handle.write("##source=ReadSimulator\n")
        handle.write("##reference=ref.fa\n")
        for gene in GENES:
            handle.write(f"##contig=<ID={gene},length={ref_lengths[gene]}>\n")
        handle.write("##phasing=none\n")
        handle.write('##INFO=<ID=TYPE,Number=A,Type=String,Description="Allele type (SNP/INS/DEL)">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tINFO\n")

        for gene, pos, mut_type, ref_seq, alt_seq in sorted_mutations:
            ref_sequence = str(refseq_dict[gene].seq)
            if mut_type == "SNP":
                handle.write(f"{gene}\t{pos+1}\t.\t{ref_seq}\t{alt_seq}\tTYPE=SNP\n")
            elif mut_type == "DEL" and pos > 0:
                preceding = ref_sequence[pos - 1]
                handle.write(f"{gene}\t{pos}\t.\t{preceding}{ref_seq}\t{preceding}\tTYPE=DEL\n")
            elif mut_type == "INS":
                handle.write(f"{gene}\t{pos+1}\t.\t{ref_seq}\t{ref_seq}{alt_seq}\tTYPE=INS\n")

    tsv_path = output_dir / "all_genes_mutations.tsv"
    with tsv_path.open("w") as handle:
        handle.write("GENE\tPOS\tREF\tMUT\tTYPE\n")
        gene_mut_map: dict[str, dict[int, tuple[str, str, str]]] = {gene: {} for gene in GENES}
        gene_ins_map: dict[str, dict[int, str]] = {gene: {} for gene in GENES}

        for gene, pos, mut_type, ref_seq, alt_seq in sorted_mutations:
            if mut_type == "SNP":
                gene_mut_map[gene][pos] = ("SNP", ref_seq, alt_seq)
            elif mut_type == "DEL":
                for idx, base in enumerate(ref_seq):
                    gene_mut_map[gene][pos + idx] = ("DEL", base, "-")
            elif mut_type == "INS":
                gene_ins_map[gene][pos] = alt_seq

        for gene in GENES:
            ref_sequence = str(refseq_dict[gene].seq)
            for current_pos, ref_base in enumerate(ref_sequence):
                if current_pos in gene_mut_map[gene]:
                    mut_type, ref_val, alt_val = gene_mut_map[gene][current_pos]
                    handle.write(f"{gene}\t{current_pos+1}\t{ref_val}\t{alt_val}\t{mut_type}\n")
                else:
                    handle.write(f"{gene}\t{current_pos+1}\t{ref_base}\t{ref_base}\tMATCH\n")

                if current_pos in gene_ins_map[gene]:
                    for base in gene_ins_map[gene][current_pos]:
                        handle.write(f"{gene}\t\t-\t{base}\tINS\n")


def simulate_reads(
    mutated_records: list[SeqRecord],
    num_reads: int,
    error_rates: list[float],
    output_dir: Path,
) -> None:
    mut_sequences = [str(record.seq) for record in mutated_records]
    mut_lengths = [len(seq) for seq in mut_sequences]

    for error in error_rates:
        print(f"Simulating reads with error rate: {error}")
        r1_records: list[SeqRecord] = []
        r2_records: list[SeqRecord] = []

        selected_genes = rd.choice(len(GENES), size=num_reads)
        for read_id, gene_idx in enumerate(selected_genes):
            fragment_extra = geom.rvs(p=1 / (FRAGMENT_LENGTH_MEAN - READ_LENGTH))
            frag_len = READ_LENGTH + fragment_extra
            max_start = mut_lengths[gene_idx] - frag_len + 1
            start = rd.choice(max_start)
            gene_name = GENES[gene_idx]
            desc = f"{gene_name}; {start}-{start + frag_len}"

            forward = list(mut_sequences[gene_idx][start : start + READ_LENGTH])
            mutate_mask = rd.choice(2, size=READ_LENGTH, p=[1 - error, error])
            forward = [rd.choice(NUCL) if mut else base for base, mut in zip(forward, mutate_mask)]
            r1 = SeqRecord(Seq("".join(forward)), id=f"Read{read_id}", description=desc)
            r1.letter_annotations["phred_quality"] = [42] * READ_LENGTH
            r1_records.append(r1)

            end_slice = mut_sequences[gene_idx][start + frag_len - READ_LENGTH : start + frag_len]
            reverse = Seq(end_slice).reverse_complement()
            mutate_mask = rd.choice(2, size=READ_LENGTH, p=[1 - error, error])
            reverse_seq = [
                rd.choice(NUCL) if mut else base for base, mut in zip(str(reverse), mutate_mask)
            ]
            r2 = SeqRecord(Seq("".join(reverse_seq)), id=f"Read{read_id}", description=desc)
            r2.letter_annotations["phred_quality"] = [42] * READ_LENGTH
            r2_records.append(r2)

        suffix = f"{error:.3f}".replace(".", "_")
        r1_path = output_dir / f"simulated_r1_err_{suffix}.fq"
        r2_path = output_dir / f"simulated_r2_err_{suffix}.fq"
        with r1_path.open("w") as handle:
            SeqIO.write(r1_records, handle, "fastq")
        with r2_path.open("w") as handle:
            SeqIO.write(r2_records, handle, "fastq")

        print(f"Saved reads to {r1_path.name} and {r2_path.name}")


def main() -> None:
    args = parse_args()
    error_rates = args.error_rate if args.error_rate else DEFAULT_ERROR
    seed_rng(args.seed)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    refseq_dict = load_references(Path(args.reference))
    mutated_records, mutation_log = mutate_sequences(
        refseq_dict,
        snp_prob=args.snp_prob,
        indel_prob=args.indel_prob,
    )

    write_reference_outputs(refseq_dict, mutated_records, mutation_log, output_dir)
    print("Exported reference FASTA/VCF/TSV files.")

    simulate_reads(mutated_records, args.num_reads, error_rates, output_dir)


if __name__ == "__main__":
    main()
