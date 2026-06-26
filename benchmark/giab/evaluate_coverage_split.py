#!/usr/bin/env python3
"""Write GIAB coverage-split confusion counts without generating figures."""

import argparse
import csv
import json
import statistics
import subprocess
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
DEFAULT_DATASET_DIR = SCRIPT_DIR / "giab_hg002_output_bowtie2_G60_min5_qfilters"
DEFAULT_OUT_DIR = SCRIPT_DIR / "vcf_out"
DEFAULT_TRUTH_VCF = (
    SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact_benchmarkable_gene_coords.vcf.gz"
)
DEFAULT_TRUTH_FALLBACK = SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact_gene_coords.vcf.gz"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate final GIAB VCF confusion counts split by coverage."
    )
    parser.add_argument("--dataset-dir", default=str(DEFAULT_DATASET_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--truth-vcf", default=None)
    parser.add_argument(
        "--n",
        type=int,
        default=None,
        help="Coverage split threshold. Defaults to rounded pooled median truth depth.",
    )
    return parser.parse_args()


def load_variants(vcf_path: Path) -> set[tuple[str, int, str, str]]:
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    variants = set()
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        chrom, pos, ref, alt = line.split("\t")
        variants.add((chrom, int(pos), ref, alt))
    return variants


def get_depth_at_position(bam_path: Path, chrom: str, pos: int) -> int:
    cmd = ["samtools", "depth", "-aa", "-r", f"{chrom}:{pos}-{pos}", str(bam_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    line = result.stdout.strip()
    if not line:
        return 0
    return int(line.split("\t")[2])


def build_depth_lookup(
    positions: set[tuple[str, int]],
    method: str,
    paradism_bam_dir: Path,
    base_bam: Path,
    dataset_prefix: str,
) -> dict[tuple[str, int], int]:
    lookup = {}
    for chrom, pos in sorted(positions):
        if method == "ParaDISM":
            bam_path = paradism_bam_dir / f"{dataset_prefix}_{chrom}.sorted.bam"
        else:
            bam_path = base_bam
        if not bam_path.exists():
            raise FileNotFoundError(f"BAM not found for {method}: {bam_path}")
        lookup[(chrom, pos)] = get_depth_at_position(bam_path, chrom, pos)
    return lookup


def choose_n(
    truth_variants: set[tuple[str, int, str, str]],
    paradism_bam_dir: Path,
    base_bam: Path,
    dataset_prefix: str,
) -> tuple[int, dict]:
    positions = {(chrom, pos) for chrom, pos, _, _ in truth_variants}
    paradism_depth = build_depth_lookup(positions, "ParaDISM", paradism_bam_dir, base_bam, dataset_prefix)
    base_depth = build_depth_lookup(positions, "BaseAligner", paradism_bam_dir, base_bam, dataset_prefix)
    pooled = list(paradism_depth.values()) + list(base_depth.values())
    if not pooled:
        return 0, {"method": "rounded pooled median truth depth", "pooled_min": 0, "pooled_median": 0, "pooled_max": 0}
    return int(round(statistics.median(pooled))), {
        "method": "rounded pooled median truth depth",
        "pooled_min": min(pooled),
        "pooled_median": statistics.median(pooled),
        "pooled_max": max(pooled),
    }


def confusion_counts(truth, called) -> dict[str, int]:
    return {
        "TP": len(truth & called),
        "FP": len(called - truth),
        "FN": len(truth - called),
    }


def split_by_coverage(truth, called, depth_lookup, n: int) -> dict[str, dict[str, int]]:
    split = {"le_n": {"TP": 0, "FP": 0, "FN": 0}, "gt_n": {"TP": 0, "FP": 0, "FN": 0}}

    for chrom, pos, ref, alt in truth & called:
        bucket = "le_n" if depth_lookup.get((chrom, pos), 0) <= n else "gt_n"
        split[bucket]["TP"] += 1
    for chrom, pos, ref, alt in called - truth:
        bucket = "le_n" if depth_lookup.get((chrom, pos), 0) <= n else "gt_n"
        split[bucket]["FP"] += 1
    for chrom, pos, ref, alt in truth - called:
        bucket = "le_n" if depth_lookup.get((chrom, pos), 0) <= n else "gt_n"
        split[bucket]["FN"] += 1

    return split


def main() -> None:
    args = parse_args()
    dataset_dir = Path(args.dataset_dir)
    dataset_prefix = dataset_dir.name
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    truth_vcf = Path(args.truth_vcf) if args.truth_vcf else DEFAULT_TRUTH_VCF
    if not truth_vcf.exists():
        truth_vcf = DEFAULT_TRUTH_FALLBACK

    truth = load_variants(truth_vcf)
    paradism_vcf = dataset_dir / "variant_calling/paradism_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
    base_vcf = dataset_dir / "variant_calling/basealigner_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
    paradism_called = load_variants(paradism_vcf)
    base_called = load_variants(base_vcf)

    paradism_bam_dir = dataset_dir / f"final_outputs/{dataset_prefix}_bam"
    base_bam = dataset_dir / "variant_calling/basealigner_raw/mapped_reads.sorted.bam"

    if args.n is None:
        n, threshold_details = choose_n(truth, paradism_bam_dir, base_bam, dataset_prefix)
    else:
        n = args.n
        threshold_details = {"method": "user-specified", "value": n}

    paradism_positions = {(chrom, pos) for chrom, pos, _, _ in (truth | paradism_called)}
    base_positions = {(chrom, pos) for chrom, pos, _, _ in (truth | base_called)}
    paradism_depth = build_depth_lookup(
        paradism_positions, "ParaDISM", paradism_bam_dir, base_bam, dataset_prefix
    )
    base_depth = build_depth_lookup(
        base_positions, "BaseAligner", paradism_bam_dir, base_bam, dataset_prefix
    )

    overall = {
        "ParaDISM": confusion_counts(truth, paradism_called),
        "BaseAligner": confusion_counts(truth, base_called),
    }
    split = {
        "ParaDISM": split_by_coverage(truth, paradism_called, paradism_depth, n),
        "BaseAligner": split_by_coverage(truth, base_called, base_depth, n),
    }

    with open(out_dir / "coverage_split_confusion_metrics.json", "w") as handle:
        json.dump(
            {
                "threshold_n_reads": n,
                "threshold_details": threshold_details,
                "overall": overall,
                "coverage_split": split,
            },
            handle,
            indent=2,
        )

    with open(out_dir / "confusion_matrices_by_coverage.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["coverage_bin", "method", "TP", "FP", "FN"])
        for bin_key, label in (("le_n", f"<= {n}"), ("gt_n", f"> {n}")):
            writer.writerow(
                [
                    label,
                    "ParaDISM",
                    split["ParaDISM"][bin_key]["TP"],
                    split["ParaDISM"][bin_key]["FP"],
                    split["ParaDISM"][bin_key]["FN"],
                ]
            )
            writer.writerow(
                [
                    label,
                    "BaseAligner",
                    split["BaseAligner"][bin_key]["TP"],
                    split["BaseAligner"][bin_key]["FP"],
                    split["BaseAligner"][bin_key]["FN"],
                ]
            )

    print(f"Truth variants: {len(truth)}")
    print(f"Coverage split threshold n={n} ({threshold_details['method']})")
    print(f"Saved coverage-split metrics to: {out_dir}")


if __name__ == "__main__":
    main()
