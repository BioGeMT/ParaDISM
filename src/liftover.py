#!/usr/bin/env python3
"""
Liftover tool for converting ParaDISM VCF/BED coordinates from gene-local
to chromosomal (GRCh38) coordinates.

ParaDISM produces variants called against a custom multi-gene reference
where each gene is a separate contig. This tool converts those gene-local
positions to standard chromosomal coordinates so they can be compared with
results from standard mapping/calling pipelines.

Usage:
    python src/liftover.py --vcf input.vcf --positions positions.txt --output lifted.vcf
    python src/liftover.py --bed input.bed --positions positions.txt --output lifted.bed
"""

import argparse
import gzip
import re
import sys
from pathlib import Path


_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# GRCh38 chromosome lengths (for VCF contig headers)
_CHR_LENGTHS = {"chr16": 90338345, "16": 90338345}

# The reference sequences in ParaDISM have this many bp of flanking
# sequence on each side beyond the gene's chromosomal span.
_FLANKING_BP = 600

_GENE_INFO_ID = "PARADISM_GENE"


def parse_positions_file(path):
    """Parse a gene positions file.

    Expected format (one gene per line):
        GENE (Human Gene) ENSID CHR:START-END:STRAND
        e.g. PKD1 (Human Gene) ENSG00000008710 16:2088708-2135898:-1

    Returns:
        dict: {gene_name: (chrom_name, start, end, strand)}
              where strand is +1 or -1, start/end are ints (1-based inclusive).
    """
    genes = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            gene = parts[0]
            # Last field is CHR:START-END:STRAND
            coord_str = parts[-1]
            m = re.match(r"(\w+):(\d+)-(\d+):(-?1)", coord_str)
            if not m:
                print(f"Warning: could not parse coordinates for {gene}: {coord_str}", file=sys.stderr)
                continue
            chrom_num = m.group(1)
            chrom = f"chr{chrom_num}"
            start = int(m.group(2))
            end = int(m.group(3))
            strand = int(m.group(4))
            genes[gene] = (chrom, start, end, strand)
    return genes


def local_to_chrom(local_pos, gene_start, gene_end, strand):
    """Convert a 1-based gene-local position to a chromosomal position.

    The gene reference has _FLANKING_BP bp of flanking on each side.
    For + strand genes, position 1 = gene_start - _FLANKING_BP on the chromosome.
    For - strand genes, the sequence is reverse-complemented, so
    position 1 = gene_end + _FLANKING_BP on the chromosome.
    """
    if strand == 1:
        return gene_start - _FLANKING_BP - 1 + local_pos
    else:
        return gene_end + _FLANKING_BP + 1 - local_pos


def local_to_chrom_bed(local_start, local_end, gene_start, gene_end, strand):
    """Convert 0-based half-open BED coordinates to chromosomal coordinates.

    Returns (chrom_start, chrom_end) in 0-based half-open format.
    """
    if strand == 1:
        # + strand: straightforward shift
        offset = gene_start - _FLANKING_BP - 1  # 0-based offset
        return (offset + local_start, offset + local_end)
    else:
        # - strand: reverse-complement, so coordinates invert
        # local_start (0-based) corresponds to higher chr position
        # local_end (0-based exclusive) corresponds to lower chr position
        offset = gene_end + _FLANKING_BP  # 1-based end of rev-comp seq
        chrom_end = offset - local_start
        chrom_start = offset - local_end
        return (chrom_start, chrom_end)


def complement_seq(seq):
    """Complement a DNA sequence (does NOT reverse)."""
    return seq.translate(_COMPLEMENT)


def _upsert_info_tag(info_field, key, value):
    """Insert or replace an INFO key=value tag in a VCF INFO column."""
    tag = f"{key}={value}"
    if info_field in ("", "."):
        return tag
    items = [item for item in info_field.split(";") if item and not item.startswith(f"{key}=")]
    items.append(tag)
    return ";".join(items)


def liftover_vcf(vcf_path, positions, output_path):
    """Liftover a VCF from gene-local to chromosomal coordinates."""
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    open_mode = "rt" if str(vcf_path).endswith(".gz") else "r"

    header_lines = []
    variant_lines = []
    contigs_seen = set()
    target_chrom = None

    has_gene_info_header = False
    has_liftover_source_header = False

    with opener(vcf_path, open_mode) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                if line.startswith(f"##INFO=<ID={_GENE_INFO_ID},"):
                    has_gene_info_header = True
                if line == "##liftover_source=ParaDISM":
                    has_liftover_source_header = True
                # Replace contig headers
                if line.startswith("##contig="):
                    m = re.match(r"##contig=<ID=(\w+)", line)
                    if m:
                        gene = m.group(1)
                        if gene in positions:
                            chrom = positions[gene][0]
                            target_chrom = chrom
                            if chrom not in contigs_seen:
                                contigs_seen.add(chrom)
                                chrom_len = _CHR_LENGTHS.get(chrom, "")
                                if chrom_len:
                                    header_lines.append(f"##contig=<ID={chrom},length={chrom_len}>")
                                else:
                                    header_lines.append(f"##contig=<ID={chrom}>")
                        else:
                            header_lines.append(line)
                else:
                    header_lines.append(line)
            elif line.startswith("#CHROM"):
                # Add liftover provenance before column header
                if not has_gene_info_header:
                    header_lines.append(
                        f"##INFO=<ID={_GENE_INFO_ID},Number=1,Type=String,"
                        "Description=\"Original ParaDISM gene/pseudogene contig before liftover\">"
                    )
                if not has_liftover_source_header:
                    header_lines.append("##liftover_source=ParaDISM")
                header_lines.append(line)
            else:
                # Variant line
                fields = line.split("\t")
                gene = fields[0]
                if gene not in positions:
                    print(f"Warning: gene '{gene}' not in positions file, skipping variant", file=sys.stderr)
                    continue
                chrom, gene_start, gene_end, strand = positions[gene]
                local_pos = int(fields[1])
                chr_pos = local_to_chrom(local_pos, gene_start, gene_end, strand)
                fields[0] = chrom
                fields[1] = str(chr_pos)
                fields[7] = _upsert_info_tag(fields[7], _GENE_INFO_ID, gene)
                if strand == -1:
                    # Complement REF and ALT
                    fields[3] = complement_seq(fields[3])
                    # ALT can have multiple alleles separated by commas
                    alts = fields[4].split(",")
                    fields[4] = ",".join(complement_seq(a) for a in alts)
                variant_lines.append(fields)

    # Sort variants by chromosome position
    variant_lines.sort(key=lambda f: int(f[1]))

    with open(output_path, "w") as out:
        for line in header_lines:
            out.write(line + "\n")
        for fields in variant_lines:
            out.write("\t".join(fields) + "\n")

    n_variants = len(variant_lines)
    print(f"Lifted {n_variants} variants → {output_path}")


def liftover_bed(bed_path, positions, output_path):
    """Liftover a BED file from gene-local to chromosomal coordinates."""
    lifted = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split("\t")
            gene = fields[0]
            if gene not in positions:
                print(f"Warning: gene '{gene}' not in positions file, skipping: {line}", file=sys.stderr)
                continue
            chrom, gene_start, gene_end, strand = positions[gene]
            local_start = int(fields[1])
            local_end = int(fields[2])
            chr_start, chr_end = local_to_chrom_bed(local_start, local_end, gene_start, gene_end, strand)
            fields[0] = chrom
            fields[1] = str(chr_start)
            fields[2] = str(chr_end)
            # Update strand field if present (column 6 in standard BED)
            if len(fields) >= 6:
                fields[5] = "-" if strand == -1 else "+"
            lifted.append(fields)

    # Sort by chromosome position
    lifted.sort(key=lambda f: int(f[1]))

    with open(output_path, "w") as out:
        for fields in lifted:
            out.write("\t".join(fields) + "\n")

    print(f"Lifted {len(lifted)} regions → {output_path}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Liftover ParaDISM VCF/BED from gene-local to chromosomal coordinates",
    )
    parser.add_argument("--positions", required=True, help="Gene positions file (e.g., PKD1_b38_pseudogene_positions.txt)")
    parser.add_argument("--output", "-o", required=True, help="Output file path")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--vcf", help="Input VCF file to liftover")
    group.add_argument("--bed", help="Input BED file to liftover")

    return parser


def run_liftover(args):
    """Run liftover from parsed arguments (used by both standalone and paradism.py)."""
    positions = parse_positions_file(args.positions)
    if not positions:
        print("Error: no gene positions parsed from file", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(positions)} gene positions:")
    for gene, (chrom, start, end, strand) in sorted(positions.items()):
        strand_str = "+" if strand == 1 else "-"
        print(f"  {gene}: {chrom}:{start}-{end} ({strand_str})")

    if args.vcf:
        liftover_vcf(args.vcf, positions, args.output)
    elif args.bed:
        liftover_bed(args.bed, positions, args.output)


def main():
    parser = build_parser()
    args = parser.parse_args()
    run_liftover(args)


if __name__ == "__main__":
    main()
