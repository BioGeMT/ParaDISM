"""
Liftover implementation for converting ParaDISM VCF/BED coordinates from
reference contigs to chromosomal coordinates.

ParaDISM produces variants called against a custom multi-gene reference
where each gene is a separate contig. The positions file gives the exact
chromosomal interval each contig was cut from, so contig positions can be
converted to chromosomal positions and compared with results from standard
mapping/calling pipelines.
"""

import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

_INTERVAL_PATTERN = re.compile(r"^(.+):(\d+)-(\d+):(-?1)$")

_SEQUENCE_ALLELE_PATTERN = re.compile(r"^[ACGTNacgtn]+$")

_GENE_INFO_ID = "PARADISM_GENE"


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def parse_positions_file(path):
    """Parse a positions file.

    Expected format (one contig per line):
        CONTIG ... CHR:START-END:STRAND
        e.g. PKD1 chr16:2088108-2136498:-1

    START-END is the exact 1-based inclusive interval the contig was cut from,
    including any flanking sequence. STRAND -1 means the contig is the reverse
    complement of that interval.

    Returns:
        dict: {contig: (chrom, start, end, strand)} in file order.
    """
    intervals = {}
    with open(path) as fh:
        for line_number, line in enumerate(fh, 1):
            fields = line.split()
            if not fields:
                continue
            match = _INTERVAL_PATTERN.match(fields[-1]) if len(fields) >= 2 else None
            if not match:
                sys.exit(f"Error: {path}:{line_number}: expected 'CONTIG CHR:START-END:STRAND', got: {line.strip()}")
            contig = fields[0]
            chrom, start, end, strand = match.group(1), int(match.group(2)), int(match.group(3)), int(match.group(4))
            if start < 1 or end < start:
                sys.exit(f"Error: {path}:{line_number}: invalid interval {fields[-1]}")
            if contig in intervals:
                sys.exit(f"Error: {path}:{line_number}: {contig} is listed more than once")
            intervals[contig] = (chrom, start, end, strand)
    if not intervals:
        sys.exit(f"Error: no contig intervals in {path}")
    return intervals


def load_oriented_contigs(ref_path, intervals):
    """Return each contig's sequence in chromosome orientation.

    Stops if a contig is missing from the reference or its length differs
    from its interval.
    """
    sequences = {record.id: str(record.seq).upper() for record in SeqIO.parse(ref_path, "fasta")}
    oriented = {}
    for contig, (chrom, start, end, strand) in intervals.items():
        if contig not in sequences:
            sys.exit(f"Error: {contig} from the positions file is not a contig in {ref_path}")
        sequence = sequences[contig]
        interval_length = end - start + 1
        if len(sequence) != interval_length:
            sys.exit(
                f"Error: {contig} is {len(sequence)} bp in {ref_path}, but its interval "
                f"{chrom}:{start}-{end} is {interval_length} bp. The positions file must give "
                "the exact interval of each contig, including any flanking sequence."
            )
        oriented[contig] = sequence if strand == 1 else reverse_complement(sequence)
    return oriented


def check_genome(genome_path, intervals, oriented):
    """Check each contig against the genome and return the chromosome lengths."""
    import pysam

    genome = pysam.FastaFile(genome_path)
    chrom_lengths = {}
    for contig, (chrom, start, end, strand) in intervals.items():
        if chrom not in genome.references:
            sys.exit(f"Error: {chrom} (interval of {contig}) is not a sequence in {genome_path}")
        genomic = genome.fetch(chrom, start - 1, end).upper()
        if genomic != oriented[contig]:
            differences = abs(len(genomic) - len(oriented[contig])) + sum(
                a != b for a, b in zip(genomic, oriented[contig])
            )
            sys.exit(
                f"Error: {contig} differs from {genome_path} {chrom}:{start}-{end} "
                f"(strand {strand}) at {differences} positions"
            )
        chrom_lengths[chrom] = genome.get_reference_length(chrom)
    return chrom_lengths


def _upsert_info_tag(info_field, key, value):
    """Insert or replace an INFO key=value tag in a VCF INFO column."""
    tag = f"{key}={value}"
    if info_field in ("", "."):
        return tag
    items = [item for item in info_field.split(";") if item and not item.startswith(f"{key}=")]
    items.append(tag)
    return ";".join(items)


def _orient_variant(fields, contig_length, strand, vcf_path):
    """Convert a VCF record from contig to chromosome orientation in place."""
    pos = int(fields[1])
    ref = fields[3]
    last = pos + len(ref) - 1
    if pos < 1 or last > contig_length:
        sys.exit(f"Error: {vcf_path}: {fields[0]}:{pos} {ref} lies outside the contig")
    if strand == 1:
        return
    alts = fields[4].split(",")
    for allele in [ref] + alts:
        if allele not in ("*", ".") and not _SEQUENCE_ALLELE_PATTERN.match(allele):
            sys.exit(f"Error: {vcf_path}: {fields[0]}:{pos} has allele {allele}, which cannot be reverse complemented")
    fields[1] = str(contig_length - last + 1)
    fields[3] = reverse_complement(ref)
    fields[4] = ",".join(allele if allele in ("*", ".") else reverse_complement(allele) for allele in alts)


def liftover_vcf(vcf_path, intervals, oriented, chrom_lengths, output_path):
    """Liftover a VCF from contig to chromosomal coordinates."""
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        sys.exit("Error: bcftools is required for VCF liftover")

    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    open_mode = "rt" if str(vcf_path).endswith(".gz") else "r"

    input_header = []
    records = []
    with opener(vcf_path, open_mode) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##contig="):
                continue
            if line.startswith("#"):
                input_header.append(line)
                continue
            fields = line.split("\t")
            contig = fields[0]
            if contig not in intervals:
                sys.exit(f"Error: {vcf_path} has variants on {contig}, which is not in the positions file")
            _orient_variant(fields, len(oriented[contig]), intervals[contig][3], vcf_path)
            records.append(fields)

    contig_order = {contig: index for index, contig in enumerate(intervals)}
    records.sort(key=lambda f: (contig_order[f[0]], int(f[1])))

    # Reverse-complemented indels are anchored on their right base; bcftools norm
    # re-anchors and left-aligns them against the contigs in chromosome orientation.
    with tempfile.TemporaryDirectory() as tmp:
        oriented_fasta = Path(tmp) / "oriented_contigs.fa"
        oriented_vcf = Path(tmp) / "oriented.vcf"
        normalized_vcf = Path(tmp) / "normalized.vcf"
        with open(oriented_fasta, "w") as out:
            for contig, sequence in oriented.items():
                out.write(f">{contig}\n{sequence}\n")
        with open(oriented_vcf, "w") as out:
            for line in input_header:
                if line.startswith("#CHROM"):
                    for contig, sequence in oriented.items():
                        out.write(f"##contig=<ID={contig},length={len(sequence)}>\n")
                out.write(line + "\n")
            for fields in records:
                out.write("\t".join(fields) + "\n")
        result = subprocess.run(
            [bcftools, "norm", "--no-version", "--check-ref", "e", "-f", str(oriented_fasta),
             "-O", "v", "-o", str(normalized_vcf), str(oriented_vcf)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            sys.exit(
                "Error: bcftools norm failed (positions and alleles below are on the contig "
                f"in chromosome orientation):\n{result.stderr}"
            )
        with open(normalized_vcf) as fh:
            normalized_lines = [line.rstrip("\n") for line in fh]

    header_lines = []
    variant_lines = []
    for line in normalized_lines:
        if line.startswith("##contig=") or line.startswith("#CHROM"):
            continue
        if line.startswith("#"):
            header_lines.append(line)
            continue
        fields = line.split("\t")
        contig = fields[0]
        chrom, start, _, _ = intervals[contig]
        fields[0] = chrom
        fields[1] = str(start + int(fields[1]) - 1)
        fields[7] = _upsert_info_tag(fields[7], _GENE_INFO_ID, contig)
        variant_lines.append(fields)

    chrom_order = {}
    for chrom, _, _, _ in intervals.values():
        chrom_order.setdefault(chrom, len(chrom_order))
    for chrom in chrom_order:
        if chrom in chrom_lengths:
            header_lines.append(f"##contig=<ID={chrom},length={chrom_lengths[chrom]}>")
        else:
            header_lines.append(f"##contig=<ID={chrom}>")
    if not any(line.startswith(f"##INFO=<ID={_GENE_INFO_ID},") for line in header_lines):
        header_lines.append(
            f"##INFO=<ID={_GENE_INFO_ID},Number=1,Type=String,"
            "Description=\"Original ParaDISM gene/pseudogene contig before liftover\">"
        )
    if "##liftover_source=ParaDISM" not in header_lines:
        header_lines.append("##liftover_source=ParaDISM")
    header_lines.append(next(line for line in input_header if line.startswith("#CHROM")))

    variant_lines.sort(key=lambda f: (chrom_order[f[0]], int(f[1])))

    with open(output_path, "w") as out:
        for line in header_lines:
            out.write(line + "\n")
        for fields in variant_lines:
            out.write("\t".join(fields) + "\n")

    print(f"Lifted {len(variant_lines)} variants → {output_path}")


def liftover_bed(bed_path, intervals, oriented, output_path):
    """Liftover a BED file from contig to chromosomal coordinates."""
    chrom_order = {}
    for chrom, _, _, _ in intervals.values():
        chrom_order.setdefault(chrom, len(chrom_order))

    lifted = []
    with open(bed_path) as fh:
        for line_number, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) > 6:
                sys.exit(f"Error: {bed_path}:{line_number}: only BED3 to BED6 are supported")
            contig = fields[0]
            if contig not in intervals:
                sys.exit(f"Error: {bed_path}:{line_number}: {contig} is not in the positions file")
            chrom, start, end, strand = intervals[contig]
            local_start = int(fields[1])
            local_end = int(fields[2])
            if not 0 <= local_start <= local_end <= len(oriented[contig]):
                sys.exit(f"Error: {bed_path}:{line_number}: interval lies outside {contig}")
            if strand == 1:
                chr_start, chr_end = start - 1 + local_start, start - 1 + local_end
            else:
                chr_start, chr_end = end - local_end, end - local_start
            fields[0] = chrom
            fields[1] = str(chr_start)
            fields[2] = str(chr_end)
            if len(fields) == 6 and strand == -1:
                fields[5] = {"+": "-", "-": "+"}.get(fields[5], fields[5])
            lifted.append(fields)

    lifted.sort(key=lambda f: (chrom_order[f[0]], int(f[1])))

    with open(output_path, "w") as out:
        for fields in lifted:
            out.write("\t".join(fields) + "\n")

    print(f"Lifted {len(lifted)} regions → {output_path}")


def run_liftover(args):
    """Run liftover from parsed paradism.py arguments."""
    intervals = parse_positions_file(args.positions)
    oriented = load_oriented_contigs(args.ref, intervals)
    chrom_lengths = check_genome(args.genome, intervals, oriented) if args.genome else {}

    print(f"Loaded {len(intervals)} contig intervals:")
    for contig, (chrom, start, end, strand) in intervals.items():
        strand_str = "+" if strand == 1 else "-"
        print(f"  {contig}: {chrom}:{start}-{end} ({strand_str})")
    if args.genome:
        print(f"Contig sequences match {args.genome}")

    if args.vcf:
        liftover_vcf(args.vcf, intervals, oriented, chrom_lengths, args.output)
    elif args.bed:
        liftover_bed(args.bed, intervals, oriented, args.output)


if __name__ == "__main__":
    raise SystemExit("Run liftover through: python paradism.py liftover ...")
