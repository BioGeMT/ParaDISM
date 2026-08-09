#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path


EXPECTED_ROWS = (
    ("PKD1", 5, 5, 0, 0),
    ("PKD1P1", 3, 0, 3, 0),
    ("PKD1P2", 2, 0, 2, 0),
    ("PKD1P3", 3, 0, 3, 0),
    ("PKD1P4", 2, 0, 2, 0),
    ("PKD1P5", 3, 2, 1, 0),
    ("PKD1P6", 2, 2, 0, 0),
)


def read_fastq_ids(path):
    with path.open("r", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()
            if not sequence or not separator or not quality:
                raise ValueError(f"Incomplete FASTQ record in {path}")
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError(f"Invalid FASTQ record in {path}")

            read_id = header[1:].split()[0]
            if read_id.endswith("/1") or read_id.endswith("/2"):
                read_id = read_id[:-2]
            yield read_id


def add_assignments(assignments, read_ids, assignment):
    for read_id in read_ids:
        previous = assignments.setdefault(read_id, assignment)
        if previous != assignment:
            raise ValueError(
                f"Read {read_id} appears in both {previous} and {assignment} outputs"
            )


def collect_assignments(assigned_dir, unassigned_read1, prefix):
    assignments = {}
    filename_prefix = f"{prefix}_"

    for fastq_path in sorted(assigned_dir.glob(f"{filename_prefix}*.fq")):
        gene = fastq_path.stem.removeprefix(filename_prefix)
        add_assignments(assignments, read_fastq_ids(fastq_path), gene)

    add_assignments(assignments, read_fastq_ids(unassigned_read1), "NONE")
    return assignments


def summarize_assignments(truth_by_read, assignments):
    missing = set(truth_by_read) - set(assignments)
    unexpected = set(assignments) - set(truth_by_read)
    if missing or unexpected:
        raise ValueError(
            "Output read IDs do not match the input: "
            f"{len(missing)} missing, {len(unexpected)} unexpected"
        )

    rows = []
    for source in sorted(set(truth_by_read.values())):
        source_reads = [
            read_id
            for read_id, true_source in truth_by_read.items()
            if true_source == source
        ]
        observed = Counter(assignments[read_id] for read_id in source_reads)
        correct = observed[source]
        unassigned = observed["NONE"]
        incorrect = len(source_reads) - correct - unassigned
        rows.append((source, len(source_reads), correct, unassigned, incorrect))
    return tuple(rows)


def write_summary(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("true_source", "input_pairs", "correct", "unassigned", "incorrect")
        )
        writer.writerows(rows)


def print_summary(rows):
    print("\nDemo assignment summary (paired-end fragments)")
    print(
        f"{'True source':<12}{'Input':>8}{'Correct':>10}"
        f"{'Unassigned':>13}{'Incorrect':>11}"
    )
    for source, total, correct, unassigned, incorrect in rows:
        print(f"{source:<12}{total:>8}{correct:>10}{unassigned:>13}{incorrect:>11}")

    total = sum(row[1] for row in rows)
    correct = sum(row[2] for row in rows)
    unassigned = sum(row[3] for row in rows)
    incorrect = sum(row[4] for row in rows)
    assigned = correct + incorrect
    precision = 100 * correct / assigned
    recall = 100 * correct / total

    print(f"{'Total':<12}{total:>8}{correct:>10}{unassigned:>13}{incorrect:>11}")
    print(f"Assignment precision: {precision:.1f}% ({correct}/{assigned})")
    print(f"Assignment recall:    {recall:.1f}% ({correct}/{total})")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate the deterministic PKD1 demo assignments"
    )
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--assigned-dir", type=Path, required=True)
    parser.add_argument("--unassigned-read1", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--summary", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    truth_by_read = {
        read_id: read_id.split("_", 1)[0]
        for read_id in read_fastq_ids(args.read1)
    }
    assignments = collect_assignments(
        args.assigned_dir,
        args.unassigned_read1,
        args.prefix,
    )
    rows = summarize_assignments(truth_by_read, assignments)
    print_summary(rows)
    write_summary(args.summary, rows)

    if rows != EXPECTED_ROWS:
        raise SystemExit(
            "ERROR: Demo assignments differ from the documented expected result"
        )
    print("Expected demo result verified.")
    print(f"Summary table: {args.summary}")


if __name__ == "__main__":
    main()
