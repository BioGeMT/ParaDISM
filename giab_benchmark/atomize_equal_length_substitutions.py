#!/usr/bin/env python3
"""
Atomize equal-length substitutions (including MNPs) into primitive SNP records.

Input:  VCF from stdin or file path
Output: VCF to stdout

Behavior:
- For biallelic records where REF and ALT reduce to an equal-length substitution
  after trimming common prefix/suffix, emit one SNP record per differing base.
- Non-equal-length variants (indels/complex) are passed through unchanged.
- Header lines are preserved and a source header is added.

This script is intended to run before SNP-only filtering.
"""

from __future__ import annotations

import argparse
import gzip
import io
import sys
from typing import Iterable, TextIO


def open_text(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def emit_source_header(line: str) -> str:
    return "##source=atomize_equal_length_substitutions.py\n" + line


def atomize_record(fields: list[str]) -> list[list[str]]:
    if len(fields) < 8:
        return []

    ref = fields[3].upper()
    alt = fields[4].upper()

    if "," in alt:
        return [fields]
    if alt == "." or "<" in alt or ">" in alt or alt == "*":
        return [fields]

    start = 0
    while start < len(ref) and start < len(alt) and ref[start] == alt[start]:
        start += 1

    end_ref = len(ref)
    end_alt = len(alt)
    while end_ref > start and end_alt > start and ref[end_ref - 1] == alt[end_alt - 1]:
        end_ref -= 1
        end_alt -= 1

    core_ref = ref[start:end_ref]
    core_alt = alt[start:end_alt]

    if len(core_ref) == 0 and len(core_alt) == 0:
        return []

    if len(core_ref) != len(core_alt) or len(core_ref) == 0:
        return [fields]

    pos0 = int(fields[1])
    out_rows: list[list[str]] = []
    for idx, (ref_base, alt_base) in enumerate(zip(core_ref, core_alt)):
        if ref_base == alt_base:
            continue
        new_fields = fields.copy()
        new_fields[1] = str(pos0 + start + idx)
        new_fields[3] = ref_base
        new_fields[4] = alt_base
        out_rows.append(new_fields)

    return out_rows


def process(handle: Iterable[str]) -> None:
    source_header_written = False
    for line in handle:
        if line.startswith("##"):
            sys.stdout.write(line)
            continue
        if line.startswith("#CHROM"):
            if not source_header_written:
                sys.stdout.write("##source=atomize_equal_length_substitutions.py\n")
                source_header_written = True
            sys.stdout.write(line)
            continue
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        for out_fields in atomize_record(fields):
            sys.stdout.write("\t".join(out_fields) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Atomize equal-length substitutions into SNPs.")
    parser.add_argument("vcf", nargs="?", default="-", help="Input VCF path or '-' for stdin")
    args = parser.parse_args()

    with open_text(args.vcf) as handle:
        process(handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

