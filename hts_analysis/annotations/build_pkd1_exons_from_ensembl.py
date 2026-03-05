#!/usr/bin/env python3
"""
Build PKD1 exon BED files from Ensembl canonical transcript ENST00000262304 (GRCh38),
then validate coordinates by matching exon sequences against local PKD1 contigs.

Outputs:
- PKD1_exons_ref.fa.annotation_ensembl.bed
- PKD1_exons_1ref.fa.annotation_ensembl.bed
- PKD1_exons_annotation_validation.tsv
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ENSEMBL_REST = "https://rest.ensembl.org"
DEFAULT_TRANSCRIPT = "ENST00000262304"
DEFAULT_ASSEMBLY = "GRCh38"
DEFAULT_CONTIG = "PKD1"


@dataclass(frozen=True)
class Exon:
    tx_rank: int
    exon_id: str
    chrom: str
    start: int  # 1-based inclusive, genome
    end: int  # 1-based inclusive, genome
    strand: int
    seq: str  # genome forward (+) sequence


def run_curl_json(url: str) -> dict:
    out = subprocess.check_output(
        ["curl", "-fsSL", "-H", "Accept: application/json", url],
        text=True,
    )
    return json.loads(out)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_contig_sequence(fasta_path: Path, contig: str) -> str:
    chunks: list[str] = []
    current = ""
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                continue
            if current == contig:
                chunks.append(line)
    if not chunks:
        raise ValueError(f"Contig '{contig}' not found in {fasta_path}")
    return "".join(chunks).upper()


def all_match_positions(haystack: str, needle: str) -> list[int]:
    pos: list[int] = []
    start = haystack.find(needle)
    while start != -1:
        pos.append(start)
        start = haystack.find(needle, start + 1)
    return pos


def load_exons_from_ensembl(transcript_id: str, assembly: str) -> tuple[dict, list[Exon]]:
    lookup_url = f"{ENSEMBL_REST}/lookup/id/{transcript_id}?expand=1"
    lookup = run_curl_json(lookup_url)
    if lookup.get("assembly_name") != assembly:
        raise RuntimeError(
            f"Ensembl returned assembly {lookup.get('assembly_name')} (expected {assembly})"
        )
    chrom = str(lookup["seq_region_name"])
    strand = int(lookup["strand"])

    raw_exons = lookup.get("Exon", [])
    if not raw_exons:
        raise RuntimeError(f"No exons found for transcript {transcript_id}")

    # Ensembl returns exons in transcript order.
    exons: list[Exon] = []
    for i, ex in enumerate(raw_exons, start=1):
        start = int(ex["start"])
        end = int(ex["end"])
        seq_url = (
            f"{ENSEMBL_REST}/sequence/region/human/"
            f"{chrom}:{start}..{end}:1?coord_system_version={assembly}"
        )
        seq_json = run_curl_json(seq_url)
        seq = seq_json["seq"].upper()
        expected_len = end - start + 1
        if len(seq) != expected_len:
            raise RuntimeError(
                f"Length mismatch for exon {ex['id']}: seq={len(seq)} expected={expected_len}"
            )
        exons.append(
            Exon(
                tx_rank=i,
                exon_id=ex["id"],
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                seq=seq,
            )
        )
    return lookup, exons


def detect_orientation_and_anchor(contig_seq: str, exons: list[Exon]) -> tuple[str, int]:
    # Forward orientation (contig bases match genome forward strand)
    forward_anchors: list[int] = []
    forward_ok = True
    for ex in exons:
        positions = all_match_positions(contig_seq, ex.seq)
        if len(positions) != 1:
            forward_ok = False
            break
        forward_anchors.append(ex.start - positions[0])

    if forward_ok and len(set(forward_anchors)) == 1:
        return "forward", forward_anchors[0]

    # Reverse orientation (contig bases match reverse complement of genome forward strand)
    reverse_anchors: list[int] = []
    reverse_ok = True
    for ex in exons:
        seq_rc = revcomp(ex.seq)
        positions = all_match_positions(contig_seq, seq_rc)
        if len(positions) != 1:
            reverse_ok = False
            break
        reverse_anchors.append(ex.end + positions[0])

    if reverse_ok and len(set(reverse_anchors)) == 1:
        return "reverse", reverse_anchors[0]

    raise RuntimeError(
        "Could not determine a single consistent orientation/anchor from unique exon matches."
    )


def project_exons_to_contig(
    contig_seq: str,
    exons: list[Exon],
    orientation: str,
) -> list[tuple[int, int, Exon]]:
    projected: list[tuple[int, int, Exon]] = []
    for ex in exons:
        if orientation == "forward":
            target = ex.seq
        else:
            target = revcomp(ex.seq)

        matches = all_match_positions(contig_seq, target)
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected unique match for {ex.exon_id}; found {len(matches)} matches"
            )
        start0 = matches[0]
        end0 = start0 + len(target)
        projected.append((start0, end0, ex))
    projected.sort(key=lambda row: (row[0], row[1]))
    return projected


def write_bed(
    out_bed: Path,
    contig: str,
    strand_symbol: str,
    projected: list[tuple[int, int, Exon]],
) -> None:
    with open(out_bed, "w", encoding="utf-8") as handle:
        for start0, end0, ex in projected:
            name = f"tx_exon_{ex.tx_rank:02d}_{ex.exon_id}"
            handle.write(f"{contig}\t{start0}\t{end0}\t{name}\t0\t{strand_symbol}\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build annotation-derived PKD1 exon BEDs and validate by sequence matching."
    )
    parser.add_argument(
        "--transcript-id",
        default=DEFAULT_TRANSCRIPT,
        help=f"Ensembl transcript ID (default: {DEFAULT_TRANSCRIPT})",
    )
    parser.add_argument(
        "--assembly",
        default=DEFAULT_ASSEMBLY,
        help=f"Genome assembly (default: {DEFAULT_ASSEMBLY})",
    )
    parser.add_argument(
        "--contig",
        default=DEFAULT_CONTIG,
        help=f"Target contig in local FASTAs (default: {DEFAULT_CONTIG})",
    )
    parser.add_argument(
        "--reference",
        action="append",
        dest="references",
        default=[],
        help="Reference FASTA path. Repeat for multiple (default: ref.fa and 1ref.fa).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Output directory (default: script directory).",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    refs = [Path(p) for p in args.references] if args.references else [
        script_dir.parent.parent / "ref.fa",
        script_dir.parent.parent / "1ref.fa",
    ]

    for ref in refs:
        if not ref.exists():
            print(f"ERROR: reference not found: {ref}", file=sys.stderr)
            return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)
    lookup, exons = load_exons_from_ensembl(args.transcript_id, args.assembly)

    lookup_cache = args.out_dir / f"{args.transcript_id}_{args.assembly}_lookup.json"
    with open(lookup_cache, "w", encoding="utf-8") as handle:
        json.dump(lookup, handle, indent=2, sort_keys=True)

    report_path = args.out_dir / "PKD1_exons_annotation_validation.tsv"
    with open(report_path, "w", encoding="utf-8") as report:
        report.write(
            "reference\tcontig\tcontig_len\torientation\tanchor_constant\texons\tstatus\n"
        )

        for ref_path in refs:
            contig_seq = load_contig_sequence(ref_path, args.contig)
            orientation, anchor = detect_orientation_and_anchor(contig_seq, exons)
            projected = project_exons_to_contig(contig_seq, exons, orientation)

            strand_symbol = "-" if int(lookup["strand"]) < 0 else "+"
            out_name = f"PKD1_exons_{ref_path.name}.annotation_ensembl.bed"
            out_bed = args.out_dir / out_name
            write_bed(out_bed, args.contig, strand_symbol, projected)

            report.write(
                f"{ref_path}\t{args.contig}\t{len(contig_seq)}\t{orientation}\t"
                f"{anchor}\t{len(projected)}\tOK\n"
            )

    print(f"Wrote lookup cache: {lookup_cache}")
    print(f"Wrote validation report: {report_path}")
    print("Wrote BED files:")
    for ref_path in refs:
        print(f"  - {args.out_dir / f'PKD1_exons_{ref_path.name}.annotation_ensembl.bed'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
