#!/usr/bin/env python3
"""
Convert GIAB HG002 truth VCF (GRCh38 chr16 coordinates) into ParaDISM gene-contig
coordinates (PKD1 / PKD1P1-6) so it can be directly compared to per-gene call VCFs.

Output is a minimal, indexable VCF with CHROM in {PKD1,PKD1P1..PKD1P6} and POS as
1-based contig coordinates matching giab_benchmark/ref.fa.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent

# Coordinate mapping: contig -> (chr16_start, chr16_end) on GRCh38 (1-based, inclusive).
# These intervals are defined to match the sequences in giab_benchmark/ref.fa and can
# include flanking sequence beyond the annotated gene models.
GENE_COORDS: dict[str, tuple[int, int]] = {
    "PKD1": (2088707, 2135898),
    "PKD1P1": (16310133, 16334190),
    "PKD1P2": (16356223, 16377507),
    "PKD1P3": (14911550, 14935708),
    "PKD1P4": (18334399, 18352476),
    "PKD1P5": (18374520, 18402014),
    "PKD1P6": (15125138, 15154873),
}


@dataclass(frozen=True)
class VariantRow:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: str
    flt: str
    gt: str


def _load_fasta_sequences(fasta_path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                sequences[current] = []
                continue
            if current is None:
                raise ValueError(f"FASTA parse error in {fasta_path}: sequence before header")
            sequences[current].append(line.strip())
    return {k: "".join(v).upper() for k, v in sequences.items()}


def _load_benchmark_intervals_chr16(bed_path: Path) -> tuple[list[int], list[tuple[int, int]]] | None:
    """
    Load and merge chr16 intervals from a GIAB benchmark BED.

    BED coordinates are 0-based, end-exclusive.
    Returns (starts, merged_intervals) where starts[i] == merged_intervals[i][0].
    """
    if not bed_path.exists():
        return None

    intervals: list[tuple[int, int]] = []
    with open(bed_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom, start_s, end_s = parts[0], parts[1], parts[2]
            if chrom != "chr16":
                continue
            try:
                start = int(start_s)
                end = int(end_s)
            except ValueError:
                continue
            if end <= start:
                continue
            intervals.append((start, end))

    if not intervals:
        return None

    intervals.sort()
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start > cur_end:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
            continue
        cur_end = max(cur_end, end)
    merged.append((cur_start, cur_end))

    starts = [s for s, _ in merged]
    return starts, merged


def _is_benchmarkable_pos(pos_1based: int, intervals: tuple[list[int], list[tuple[int, int]]] | None) -> bool:
    if intervals is None:
        return True
    if pos_1based <= 0:
        return False

    starts, merged = intervals
    pos0 = pos_1based - 1
    idx = bisect_right(starts, pos0) - 1
    if idx < 0:
        return False
    start0, end0 = merged[idx]
    return start0 <= pos0 < end0


def _find_gene_for_chr16_pos(pos_1based: int) -> tuple[str, int] | None:
    for gene, (start, end) in GENE_COORDS.items():
        if start <= pos_1based <= end:
            return gene, pos_1based - start + 1
    return None


def _bgzip_and_index(vcf_path: Path, out_vcf_gz: Path) -> None:
    """
    bgzip-compress and bcftools-index the VCF.
    Prefers system tools; falls back to `conda run -n paradism_env`.
    """
    bgzip = shutil.which("bgzip")
    bcftools = shutil.which("bcftools")

    if not (bgzip and bcftools):
        # Avoid writing binary data through `conda run` (it can print to stdout and corrupt output).
        conda = shutil.which("conda")
        if not conda:
            raise RuntimeError(
                "bgzip/bcftools not found on PATH and conda not found; cannot create .vcf.gz + index"
            )

        bgzip_p = subprocess.run(
            [conda, "run", "-n", "paradism_env", "which", "bgzip"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        bcftools_p = subprocess.run(
            [conda, "run", "-n", "paradism_env", "which", "bcftools"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()

        bgzip = bgzip_p or None
        bcftools = bcftools_p or None

    if not (bgzip and bcftools):
        raise RuntimeError("Could not locate bgzip/bcftools executables")

    with open(out_vcf_gz, "wb") as out_f:
        subprocess.run([bgzip, "-c", str(vcf_path)], check=True, stdout=out_f)
    subprocess.run([bcftools, "index", str(out_vcf_gz)], check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert GIAB HG002 truth VCF (chr16) to PKD1/PKD1P* gene-contig coordinates.",
    )
    parser.add_argument(
        "--truth-vcf",
        type=Path,
        default=(
            SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact_benchmarkable.vcf.gz"
            if (SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact_benchmarkable.vcf.gz").exists()
            else SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact.vcf.gz"
        ),
        help="Input truth VCF (.vcf.gz). Default: benchmarkable truth if present, else exact truth.",
    )
    parser.add_argument(
        "--benchmark-bed",
        type=Path,
        default=SCRIPT_DIR / "giab_hg002_vcf/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        help="GIAB benchmark BED for filtering to benchmarkable regions (chr16 only).",
    )
    parser.add_argument(
        "--ref-fasta",
        type=Path,
        default=SCRIPT_DIR / "ref.fa",
        help="Reference FASTA used by ParaDISM per-gene contigs (giab_benchmark/ref.fa).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output VCF.GZ path. Default: <truth_vcf stem>_gene_coords.vcf.gz in giab_hg002_vcf/.",
    )
    args = parser.parse_args()

    truth_vcf: Path = args.truth_vcf
    if not truth_vcf.exists():
        print(f"ERROR: truth VCF not found: {truth_vcf}", file=sys.stderr)
        return 2

    benchmark_intervals = _load_benchmark_intervals_chr16(args.benchmark_bed)
    ref_seqs = _load_fasta_sequences(args.ref_fasta)
    contig_order = [name for name in ref_seqs.keys() if name in GENE_COORDS]
    if set(GENE_COORDS.keys()) - set(ref_seqs.keys()):
        missing = sorted(set(GENE_COORDS.keys()) - set(ref_seqs.keys()))
        print(f"ERROR: ref fasta missing contigs: {', '.join(missing)}", file=sys.stderr)
        return 2

    if args.output is None:
        out_dir = SCRIPT_DIR / "giab_hg002_vcf"
        out_dir.mkdir(exist_ok=True)
        stem = truth_vcf.name
        if stem.endswith(".vcf.gz"):
            stem = stem[: -len(".vcf.gz")]
        out_vcf_gz = out_dir / f"{stem}_gene_coords.vcf.gz"
    else:
        out_vcf_gz = args.output

    rows: list[VariantRow] = []
    skipped_not_chr16 = 0
    skipped_not_benchmarkable = 0
    skipped_outside_gene = 0
    skipped_non_simple_snp = 0
    ref_mismatch = 0

    with gzip.open(truth_vcf, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom = parts[0]
            if chrom != "chr16":
                skipped_not_chr16 += 1
                continue
            pos = int(parts[1])
            ref = parts[3].upper()
            alt = parts[4].upper()

            if "," in alt or len(ref) != 1 or len(alt) != 1:
                skipped_non_simple_snp += 1
                continue

            if not _is_benchmarkable_pos(pos, benchmark_intervals):
                skipped_not_benchmarkable += 1
                continue

            mapped = _find_gene_for_chr16_pos(pos)
            if mapped is None:
                skipped_outside_gene += 1
                continue
            gene, gene_pos = mapped

            ref_base = ref_seqs[gene][gene_pos - 1]
            if ref_base != ref:
                ref_mismatch += 1
                continue

            qual = parts[5]
            flt = parts[6]
            fmt = parts[8].split(":")
            sample = parts[9].split(":")
            gt = "."
            if "GT" in fmt:
                idx = fmt.index("GT")
                if idx < len(sample):
                    gt = sample[idx]

            rows.append(VariantRow(gene, gene_pos, ref, alt, qual, flt, gt))

    def _sort_key(r: VariantRow) -> tuple[int, int]:
        return (contig_order.index(r.chrom), r.pos)

    rows.sort(key=_sort_key)

    # Write uncompressed VCF, then bgzip + index for tabix/CSI compatibility.
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_vcf = Path(tmpdir) / "truth_gene_coords.vcf"
        with open(tmp_vcf, "w") as out:
            out.write("##fileformat=VCFv4.2\n")
            out.write("##source=make_truth_gene_coords_vcf.py\n")
            out.write(f"##reference={args.ref_fasta}\n")
            for contig in contig_order:
                out.write(f"##contig=<ID={contig},length={len(ref_seqs[contig])}>\n")
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")

            for r in rows:
                out.write(
                    f"{r.chrom}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t{r.qual}\t{r.flt}\t.\tGT\t{r.gt}\n"
                )

        _bgzip_and_index(tmp_vcf, out_vcf_gz)

    print(f"Wrote: {out_vcf_gz}")
    print(f"Variants: {len(rows)}")
    if benchmark_intervals is None:
        print("Note: benchmark BED not found; did not restrict to benchmarkable regions.")
    if skipped_outside_gene or skipped_not_benchmarkable or skipped_non_simple_snp or ref_mismatch:
        print(
            "Skipped:"
            f" not_chr16={skipped_not_chr16},"
            f" not_benchmarkable={skipped_not_benchmarkable},"
            f" outside_gene={skipped_outside_gene},"
            f" non_simple_snp={skipped_non_simple_snp},"
            f" ref_mismatch={ref_mismatch}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
