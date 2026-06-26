#!/usr/bin/env python3
"""
Build a multi-sequence ParaDISM reference FASTA from Ensembl REST.

Use gene symbols for ordinary cases:

    python benchmark/references/make_reference.py \
      --gene HBA1 --gene HBA2 \
      --out benchmark/references/hba_pair.fa

Use explicit regions when an exact locus is needed:

    python benchmark/references/make_reference.py \
      --region PKD1=16:2088707-2135898:1 \
      --out benchmark/references/custom.fa
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


ENSEMBL_REST = "https://rest.ensembl.org"


@dataclass(frozen=True)
class Region:
    label: str
    seq_region: str
    start: int
    end: int
    strand: int


def request_text(url: str, accept: str) -> str:
    request = Request(url, headers={"Accept": accept, "User-Agent": "ParaDISM-reference-helper"})
    try:
        with urlopen(request, timeout=60) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {exc.code} for {url}: {detail}") from exc
    except URLError as exc:
        raise RuntimeError(f"Failed to reach Ensembl REST for {url}: {exc}") from exc


def lookup_gene(symbol: str, species: str, flank: int) -> Region:
    url = f"{ENSEMBL_REST}/lookup/symbol/{quote(species)}/{quote(symbol)}?{urlencode({'expand': 0})}"
    payload = json.loads(request_text(url, "application/json"))

    start = int(payload["start"]) - flank
    end = int(payload["end"]) + flank
    if start < 1:
        start = 1

    return Region(
        label=symbol,
        seq_region=str(payload["seq_region_name"]),
        start=start,
        end=end,
        strand=int(payload.get("strand", 1)),
    )


def parse_region(raw: str) -> Region:
    if "=" not in raw:
        raise ValueError(f"Region must use LABEL=SEQ_REGION:START-END[:STRAND], got: {raw}")
    label, region_text = raw.split("=", 1)
    parts = region_text.split(":")
    if len(parts) not in {2, 3}:
        raise ValueError(f"Region must use LABEL=SEQ_REGION:START-END[:STRAND], got: {raw}")

    seq_region = parts[0]
    start_text, end_text = parts[1].split("-", 1)
    strand = int(parts[2]) if len(parts) == 3 else 1
    if strand not in {-1, 1}:
        raise ValueError(f"Region strand must be 1 or -1, got: {raw}")

    start = int(start_text.replace(",", ""))
    end = int(end_text.replace(",", ""))
    if start < 1 or end < start:
        raise ValueError(f"Invalid region coordinates: {raw}")

    return Region(label=label, seq_region=seq_region, start=start, end=end, strand=strand)


def fetch_sequence(region: Region, species: str, assembly: str | None) -> str:
    region_spec = f"{region.seq_region}:{region.start}..{region.end}:{region.strand}"
    query = {}
    if assembly:
        query["coord_system_version"] = assembly

    suffix = f"?{urlencode(query)}" if query else ""
    url = f"{ENSEMBL_REST}/sequence/region/{quote(species)}/{quote(region_spec)}{suffix}"
    fasta = request_text(url, "text/x-fasta")
    sequence_lines = [line.strip() for line in fasta.splitlines() if line and not line.startswith(">")]
    sequence = "".join(sequence_lines).upper()
    if not sequence:
        raise RuntimeError(f"Ensembl returned an empty sequence for {region.label}: {region_spec}")
    return sequence


def write_fasta(regions: list[Region], sequences: list[str], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        for region, sequence in zip(regions, sequences, strict=True):
            handle.write(f">{region.label}\n")
            for i in range(0, len(sequence), 80):
                handle.write(sequence[i : i + 80] + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a ParaDISM multi-contig reference FASTA from Ensembl.",
    )
    parser.add_argument("--gene", action="append", default=[], help="Gene symbol to fetch, repeatable.")
    parser.add_argument(
        "--region",
        action="append",
        default=[],
        help="Explicit region as LABEL=SEQ_REGION:START-END[:STRAND], repeatable.",
    )
    parser.add_argument("--species", default="homo_sapiens", help="Ensembl species name.")
    parser.add_argument("--assembly", default="GRCh38", help="Coordinate-system version, or empty string.")
    parser.add_argument("--flank", type=int, default=0, help="Bases to add around gene-symbol intervals.")
    parser.add_argument("--sleep", type=float, default=0.2, help="Delay between Ensembl requests.")
    parser.add_argument("--out", type=Path, required=True, help="Output FASTA path.")
    args = parser.parse_args()

    if not args.gene and not args.region:
        parser.error("Provide at least one --gene or --region.")
    if args.flank < 0:
        parser.error("--flank must be >= 0.")

    try:
        regions = [lookup_gene(gene, args.species, args.flank) for gene in args.gene]
        regions.extend(parse_region(region) for region in args.region)
        assembly = args.assembly or None

        sequences: list[str] = []
        for index, region in enumerate(regions):
            if index:
                time.sleep(args.sleep)
            sequence = fetch_sequence(region, args.species, assembly)
            sequences.append(sequence)
            print(
                f"{region.label}: {region.seq_region}:{region.start}-{region.end}:{region.strand} "
                f"({len(sequence)} bp)"
            )

        write_fasta(regions, sequences, args.out)
        print(f"Wrote: {args.out}")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
