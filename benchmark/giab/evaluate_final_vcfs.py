#!/usr/bin/env python3
"""Evaluate final GIAB SNP VCFs and write CSV/JSON metrics only."""

import argparse
import csv
import json
import subprocess
from bisect import bisect_right
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
DEFAULT_DATASET_DIR = SCRIPT_DIR / "giab_hg002_output_bowtie2_G60_min5_qfilters"
DEFAULT_OUT_DIR = SCRIPT_DIR / "vcf_out"
PARADISM_VCF_SUBPATH = (
    "variant_calling/paradism_raw/"
    "variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
)
BASE_VCF_SUBPATH = (
    "variant_calling/basealigner_raw/"
    "variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
)
DEFAULT_TRUTH_VCF = SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact_benchmarkable.vcf.gz"
DEFAULT_TRUTH_FALLBACK = SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact.vcf.gz"
DEFAULT_BENCHMARK_BED = (
    SCRIPT_DIR / "giab_hg002_vcf/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
)

# GRCh38 1-based inclusive coordinates matching benchmark/references/pkd1_panel.fa.
GENE_COORDS = {
    "PKD1": (2088707, 2135898),
    "PKD1P1": (16310133, 16334190),
    "PKD1P2": (16356223, 16377507),
    "PKD1P3": (14911550, 14935708),
    "PKD1P4": (18334399, 18352476),
    "PKD1P5": (18374520, 18402014),
    "PKD1P6": (15125138, 15154873),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate final ParaDISM/base-aligner VCFs against GIAB truth."
    )
    parser.add_argument("--dataset-dir", default=str(DEFAULT_DATASET_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--truth-vcf", default=None)
    parser.add_argument("--benchmark-bed", default=str(DEFAULT_BENCHMARK_BED))
    return parser.parse_args()


def load_benchmark_intervals_chr16(benchmark_bed: Path) -> tuple[list[int], list[tuple[int, int]]] | None:
    if not benchmark_bed.exists():
        return None

    intervals = []
    with open(benchmark_bed) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != "chr16":
                continue
            start = int(parts[1])
            end = int(parts[2])
            if end > start:
                intervals.append((start, end))

    if not intervals:
        return None

    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start > prev_end:
            merged.append((start, end))
        else:
            merged[-1] = (prev_start, max(prev_end, end))

    return [start for start, _ in merged], merged


def is_benchmarkable_chr16_pos(
    pos_1based: int,
    benchmark_intervals: tuple[list[int], list[tuple[int, int]]] | None,
) -> bool:
    if benchmark_intervals is None:
        return True
    starts, intervals = benchmark_intervals
    pos0 = pos_1based - 1
    idx = bisect_right(starts, pos0) - 1
    if idx < 0:
        return False
    start0, end0 = intervals[idx]
    return start0 <= pos0 < end0


def query_vcf(vcf_path: Path) -> list[tuple[str, int, str, str]]:
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        chrom, pos, ref, alt = line.split("\t")
        rows.append((chrom, int(pos), ref, alt))
    return rows


def load_truth_variants(
    truth_vcf: Path,
    benchmark_intervals: tuple[list[int], list[tuple[int, int]]] | None,
) -> set[tuple[str, int, str, str]]:
    variants = set()
    for chrom, pos, ref, alt in query_vcf(truth_vcf):
        if chrom != "chr16":
            continue
        if not is_benchmarkable_chr16_pos(pos, benchmark_intervals):
            continue
        for gene, (start, end) in GENE_COORDS.items():
            if start <= pos <= end:
                variants.add((gene, pos - start + 1, ref, alt))
                break
    return variants


def load_called_variants(
    vcf_path: Path,
    benchmark_intervals: tuple[list[int], list[tuple[int, int]]] | None,
) -> set[tuple[str, int, str, str]]:
    variants = set()
    for chrom, pos, ref, alt in query_vcf(vcf_path):
        if chrom not in GENE_COORDS:
            continue
        gene_start, _ = GENE_COORDS[chrom]
        genome_pos = gene_start + pos - 1
        if is_benchmarkable_chr16_pos(genome_pos, benchmark_intervals):
            variants.add((chrom, pos, ref, alt))
    return variants


def metric_counts(
    truth: set[tuple[str, int, str, str]],
    called: set[tuple[str, int, str, str]],
) -> dict:
    tp = len(truth & called)
    fp = len(called - truth)
    fn = len(truth - called)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def write_overall(out_dir: Path, truth, paradism_called, base_called) -> None:
    metrics = {
        "truth_variants": len(truth),
        "methods": {
            "ParaDISM": metric_counts(truth, paradism_called),
            "BaseAligner": metric_counts(truth, base_called),
        },
    }

    with open(out_dir / "variant_calling_metrics.json", "w") as handle:
        json.dump(metrics, handle, indent=2)

    with open(out_dir / "variant_calling_metrics.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "TP", "FP", "FN", "precision", "recall", "f1"])
        for method, values in metrics["methods"].items():
            writer.writerow(
                [
                    method,
                    values["TP"],
                    values["FP"],
                    values["FN"],
                    round(values["precision"], 6),
                    round(values["recall"], 6),
                    round(values["f1"], 6),
                ]
            )

    with open(out_dir / "confusion_matrices_overall.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "TP", "FP", "FN"])
        for method, values in metrics["methods"].items():
            writer.writerow([method, values["TP"], values["FP"], values["FN"]])


def write_per_gene(out_dir: Path, truth, paradism_called, base_called) -> None:
    per_gene_dir = out_dir / "per_gene_metrics"
    per_gene_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for gene in GENE_COORDS:
        gene_truth = {variant for variant in truth if variant[0] == gene}
        for method, called in (
            ("ParaDISM", paradism_called),
            ("BaseAligner", base_called),
        ):
            values = metric_counts(gene_truth, {variant for variant in called if variant[0] == gene})
            rows.append(
                {
                    "method": method,
                    "gene": gene,
                    "truth_count": len(gene_truth),
                    **values,
                }
            )

    with open(per_gene_dir / "per_gene_metrics.csv", "w", newline="") as handle:
        fieldnames = ["method", "gene", "truth_count", "TP", "FP", "FN", "precision", "recall", "f1"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "precision": round(row["precision"], 6),
                    "recall": round(row["recall"], 6),
                    "f1": round(row["f1"], 6),
                }
            )

    with open(per_gene_dir / "per_gene_metrics.json", "w") as handle:
        json.dump(rows, handle, indent=2)


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    truth_vcf = Path(args.truth_vcf) if args.truth_vcf else DEFAULT_TRUTH_VCF
    if not truth_vcf.exists():
        truth_vcf = DEFAULT_TRUTH_FALLBACK

    benchmark_intervals = load_benchmark_intervals_chr16(Path(args.benchmark_bed))
    truth = load_truth_variants(truth_vcf, benchmark_intervals)
    paradism_called = load_called_variants(dataset_dir / PARADISM_VCF_SUBPATH, benchmark_intervals)
    base_called = load_called_variants(dataset_dir / BASE_VCF_SUBPATH, benchmark_intervals)

    write_overall(out_dir, truth, paradism_called, base_called)
    write_per_gene(out_dir, truth, paradism_called, base_called)

    print(f"Truth variants: {len(truth)}")
    print(f"Saved metrics to: {out_dir}")


if __name__ == "__main__":
    main()
