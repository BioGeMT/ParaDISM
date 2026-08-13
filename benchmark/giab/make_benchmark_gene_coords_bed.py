#!/usr/bin/env python3
"""
Convert GIAB chr16 benchmark BED intervals into ParaDISM gene-contig BED intervals.

Input BED is expected to be GRCh38 coordinates (0-based, end-exclusive).
Output BED is written in PKD1/PKD1P* contig coordinates (0-based, end-exclusive),
matching benchmark/references/pkd1_panel.fa.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from pkd1_panel_coordinates import EVALUATION_REGIONS, PANEL_COORDINATES


def load_and_merge_chr16_bed_intervals(path: Path) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    with open(path, "r") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if chrom != "chr16":
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue
            intervals.append((start, end))

    if not intervals:
        return []

    intervals.sort()
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start > cur_end:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    merged.append((cur_start, cur_end))
    return merged


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = intervals[0]
    for start, end in intervals[1:]:
        if start > cur_end:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    merged.append((cur_start, cur_end))
    return merged


def map_chr16_to_gene_coords(
    chr16_intervals: list[tuple[int, int]],
) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {
        gene: [] for gene in PANEL_COORDINATES
    }
    for gene, coordinate in PANEL_COORDINATES.items():
        evaluation_start0, evaluation_end0 = EVALUATION_REGIONS[gene]
        mapped: list[tuple[int, int]] = []
        for start0, end0 in chr16_intervals:
            overlap_start0 = max(start0, evaluation_start0)
            overlap_end0 = min(end0, evaluation_end0)
            interval = coordinate.genomic_to_contig_interval(
                overlap_start0, overlap_end0
            )
            if interval is None:
                continue
            mapped.append(interval)

        out[gene] = merge_intervals(mapped)
    return out


def total_bp(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in intervals)


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description="Create benchmarkable BED in PKD1/PKD1P* gene coordinates.",
    )
    parser.add_argument(
        "--benchmark-bed",
        type=Path,
        default=script_dir / "giab_hg002_vcf/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        help="Input GIAB benchmark BED (chr16 on GRCh38).",
    )
    parser.add_argument(
        "--output-bed",
        type=Path,
        default=script_dir / "giab_hg002_vcf/HG002_PKD1_genes_benchmarkable_gene_coords.bed",
        help="Output BED in PKD1/PKD1P* contig coordinates.",
    )
    args = parser.parse_args()

    if not args.benchmark_bed.exists():
        print(f"ERROR: benchmark BED not found: {args.benchmark_bed}")
        return 2

    chr16_intervals = load_and_merge_chr16_bed_intervals(args.benchmark_bed)
    if not chr16_intervals:
        print("ERROR: no chr16 intervals found in benchmark BED.")
        return 2

    gene_intervals = map_chr16_to_gene_coords(chr16_intervals)

    args.output_bed.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_bed, "w") as out:
        for gene in PANEL_COORDINATES:
            for start0, end0 in gene_intervals[gene]:
                out.write(f"{gene}\t{start0}\t{end0}\n")

    genes_with_intervals = sum(
        1 for gene in PANEL_COORDINATES if gene_intervals[gene]
    )
    total_intervals = sum(len(v) for v in gene_intervals.values())
    total_bases = sum(total_bp(v) for v in gene_intervals.values())

    print(f"Wrote: {args.output_bed}")
    print(
        "Genes with benchmarkable intervals:"
        f" {genes_with_intervals}/{len(PANEL_COORDINATES)}"
    )
    print(f"Total intervals: {total_intervals}")
    print(f"Total benchmarkable bases (gene coords): {total_bases}")
    for gene in PANEL_COORDINATES:
        intervals = gene_intervals[gene]
        if not intervals:
            continue
        print(f"  {gene}: intervals={len(intervals)} bases={total_bp(intervals)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
