#!/usr/bin/env python3
"""
Create coverage-split confusion matrices for final filtered SNP calls.
"""

import argparse
import csv
import json
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / "vcf_out"

DATASET_DIR = SCRIPT_DIR / "giab_hg002_output_bowtie2_G60_min5_qfilters_genome10x_balanced_v1"
DATASET_PREFIX = DATASET_DIR.name

TRUTH_VCF = SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact_benchmarkable_gene_coords.vcf.gz"
if not TRUTH_VCF.exists():
    TRUTH_VCF = SCRIPT_DIR / "giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact_gene_coords.vcf.gz"

PARADISM_VCF = DATASET_DIR / "variant_calling/paradism_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
BASE_VCF = DATASET_DIR / "variant_calling/basealigner_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"

PARADISM_BAM_DIR = DATASET_DIR / f"final_outputs/{DATASET_PREFIX}_bam"
BASE_BAM = DATASET_DIR / "variant_calling/basealigner_raw/mapped_reads.sorted.bam"

FILTER_LABEL = "AO>=3, QUAL>=10, SAF>=1, SAR>=1 + benchmark BED"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create confusion matrices split by coverage <=n vs >n."
    )
    parser.add_argument(
        "--dataset-dir",
        default=str(DATASET_DIR),
        help=f"Run directory containing variant_calling outputs (default: {DATASET_DIR})",
    )
    parser.add_argument(
        "--out-dir",
        default=str(OUTPUT_DIR),
        help=f"Output directory (default: {OUTPUT_DIR})",
    )
    parser.add_argument(
        "--n",
        type=int,
        default=None,
        help="Coverage split threshold n. If omitted, chooses rounded pooled median truth depth.",
    )
    return parser.parse_args()


def load_variants(vcf_path: Path) -> set[tuple[str, int, str, str]]:
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    variants: set[tuple[str, int, str, str]] = set()
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


def get_paradism_bam(contig: str) -> Path:
    return PARADISM_BAM_DIR / f"{DATASET_PREFIX}_{contig}.sorted.bam"


def build_depth_lookup(method: str, positions: set[tuple[str, int]]) -> dict[tuple[str, int], int]:
    lookup: dict[tuple[str, int], int] = {}
    for chrom, pos in sorted(positions):
        if method == "ParaDISM":
            bam_path = get_paradism_bam(chrom)
        else:
            bam_path = BASE_BAM
        if not bam_path.exists():
            raise FileNotFoundError(f"BAM not found for {method}: {bam_path}")
        lookup[(chrom, pos)] = get_depth_at_position(bam_path, chrom, pos)
    return lookup


def choose_n_from_truth(truth_variants: set[tuple[str, int, str, str]]) -> tuple[int, dict]:
    truth_positions = {(chrom, pos) for chrom, pos, _, _ in truth_variants}
    paradism_depth = build_depth_lookup("ParaDISM", truth_positions)
    base_depth = build_depth_lookup("BaseAligner", truth_positions)
    pooled = list(paradism_depth.values()) + list(base_depth.values())
    if not pooled:
        return 0, {"method": "rounded pooled median truth depth", "pooled_min": 0, "pooled_median": 0.0, "pooled_max": 0}
    n = int(round(float(np.median(pooled))))
    return n, {
        "method": "rounded pooled median truth depth",
        "pooled_min": int(min(pooled)),
        "pooled_median": float(np.median(pooled)),
        "pooled_max": int(max(pooled)),
    }


def confusion_counts(
    truth_variants: set[tuple[str, int, str, str]],
    called_variants: set[tuple[str, int, str, str]],
) -> dict[str, int]:
    return {
        "TP": len(truth_variants & called_variants),
        "FP": len(called_variants - truth_variants),
        "FN": len(truth_variants - called_variants),
    }


def split_by_coverage(
    truth_variants: set[tuple[str, int, str, str]],
    called_variants: set[tuple[str, int, str, str]],
    depth_lookup: dict[tuple[str, int], int],
    n: int,
) -> dict[str, dict[str, int]]:
    split = {"le_n": {"TP": 0, "FP": 0, "FN": 0}, "gt_n": {"TP": 0, "FP": 0, "FN": 0}}

    for chrom, pos, ref, alt in truth_variants & called_variants:
        bucket = "le_n" if depth_lookup.get((chrom, pos), 0) <= n else "gt_n"
        split[bucket]["TP"] += 1
    for chrom, pos, ref, alt in called_variants - truth_variants:
        bucket = "le_n" if depth_lookup.get((chrom, pos), 0) <= n else "gt_n"
        split[bucket]["FP"] += 1
    for chrom, pos, ref, alt in truth_variants - called_variants:
        bucket = "le_n" if depth_lookup.get((chrom, pos), 0) <= n else "gt_n"
        split[bucket]["FN"] += 1
    return split


def plot_confusion_matrix(title: str, paradism: dict[str, int], base: dict[str, int], output_png: Path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for idx, (method, counts) in enumerate([("ParaDISM", paradism), ("Base Aligner", base)]):
        tp = counts["TP"]
        fp = counts["FP"]
        fn = counts["FN"]
        matrix = np.array([[tp, fp], [fn, 0]])
        labels = [["TP", "FP"], ["FN", "—"]]

        ax = axes[idx]
        ax.imshow(matrix, cmap="Blues", vmin=0, vmax=max(tp + fp + fn, 1))
        for i in range(2):
            for j in range(2):
                if labels[i][j] == "—":
                    text, color = "—", "gray"
                else:
                    text = f"{labels[i][j]}\n{matrix[i, j]}"
                    color = "white" if matrix[i, j] > (tp + fp + fn) / 2 else "black"
                ax.text(j, i, text, ha="center", va="center", fontsize=14, color=color)
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_title(method, fontsize=12, fontweight="bold")

    plt.suptitle(f"{title}\n({FILTER_LABEL})", fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    args = parse_args()
    global DATASET_DIR, DATASET_PREFIX, PARADISM_VCF, BASE_VCF, PARADISM_BAM_DIR, BASE_BAM
    DATASET_DIR = Path(args.dataset_dir)
    DATASET_PREFIX = DATASET_DIR.name
    PARADISM_VCF = DATASET_DIR / "variant_calling/paradism_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
    BASE_VCF = DATASET_DIR / "variant_calling/basealigner_raw/variants_simple_snps_acgt_benchmarkable_regions_final.vcf.gz"
    PARADISM_BAM_DIR = DATASET_DIR / f"final_outputs/{DATASET_PREFIX}_bam"
    BASE_BAM = DATASET_DIR / "variant_calling/basealigner_raw/mapped_reads.sorted.bam"
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    truth = load_variants(TRUTH_VCF)
    paradism_called = load_variants(PARADISM_VCF)
    base_called = load_variants(BASE_VCF)

    overall = {
        "ParaDISM": confusion_counts(truth, paradism_called),
        "BaseAligner": confusion_counts(truth, base_called),
    }

    if args.n is None:
        n, n_info = choose_n_from_truth(truth)
    else:
        n = int(args.n)
        n_info = {"method": "user-specified", "value": n}

    paradism_positions = {(chrom, pos) for chrom, pos, _, _ in (truth | paradism_called)}
    base_positions = {(chrom, pos) for chrom, pos, _, _ in (truth | base_called)}
    paradism_depth = build_depth_lookup("ParaDISM", paradism_positions)
    base_depth = build_depth_lookup("BaseAligner", base_positions)

    split = {
        "ParaDISM": split_by_coverage(truth, paradism_called, paradism_depth, n),
        "BaseAligner": split_by_coverage(truth, base_called, base_depth, n),
    }

    print(f"Truth variants: {len(truth)}")
    print(f"Coverage split threshold n={n} ({n_info['method']})")
    print(f"Overall ParaDISM: TP={overall['ParaDISM']['TP']} FP={overall['ParaDISM']['FP']} FN={overall['ParaDISM']['FN']}")
    print(f"Overall Base:     TP={overall['BaseAligner']['TP']} FP={overall['BaseAligner']['FP']} FN={overall['BaseAligner']['FN']}")
    print(f"<= {n} reads ParaDISM: TP={split['ParaDISM']['le_n']['TP']} FP={split['ParaDISM']['le_n']['FP']} FN={split['ParaDISM']['le_n']['FN']}")
    print(f"<= {n} reads Base:     TP={split['BaseAligner']['le_n']['TP']} FP={split['BaseAligner']['le_n']['FP']} FN={split['BaseAligner']['le_n']['FN']}")
    print(f"> {n} reads ParaDISM: TP={split['ParaDISM']['gt_n']['TP']} FP={split['ParaDISM']['gt_n']['FP']} FN={split['ParaDISM']['gt_n']['FN']}")
    print(f"> {n} reads Base:     TP={split['BaseAligner']['gt_n']['TP']} FP={split['BaseAligner']['gt_n']['FP']} FN={split['BaseAligner']['gt_n']['FN']}")

    json_path = out_dir / "coverage_split_confusion_metrics.json"
    with open(json_path, "w") as f:
        json.dump(
            {
                "threshold_n_reads": n,
                "threshold_details": n_info,
                "overall": overall,
                "coverage_split": split,
            },
            f,
            indent=2,
        )
    print(f"Saved {json_path}")

    csv_path = out_dir / "confusion_matrices_by_coverage.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["coverage_bin", "method", "TP", "FP", "FN"])
        writer.writerow([f"<= {n}", "ParaDISM", split["ParaDISM"]["le_n"]["TP"], split["ParaDISM"]["le_n"]["FP"], split["ParaDISM"]["le_n"]["FN"]])
        writer.writerow([f"<= {n}", "BaseAligner", split["BaseAligner"]["le_n"]["TP"], split["BaseAligner"]["le_n"]["FP"], split["BaseAligner"]["le_n"]["FN"]])
        writer.writerow([f"> {n}", "ParaDISM", split["ParaDISM"]["gt_n"]["TP"], split["ParaDISM"]["gt_n"]["FP"], split["ParaDISM"]["gt_n"]["FN"]])
        writer.writerow([f"> {n}", "BaseAligner", split["BaseAligner"]["gt_n"]["TP"], split["BaseAligner"]["gt_n"]["FP"], split["BaseAligner"]["gt_n"]["FN"]])
    print(f"Saved {csv_path}")

    plot_confusion_matrix(
        f"Confusion Matrices - G60 (coverage <= {n})",
        split["ParaDISM"]["le_n"],
        split["BaseAligner"]["le_n"],
        out_dir / f"confusion_matrix_G60_cov_le_{n}.png",
    )
    plot_confusion_matrix(
        f"Confusion Matrices - G60 (coverage > {n})",
        split["ParaDISM"]["gt_n"],
        split["BaseAligner"]["gt_n"],
        out_dir / f"confusion_matrix_G60_cov_gt_{n}.png",
    )
    print(f"Saved {out_dir / f'confusion_matrix_G60_cov_le_{n}.png'}")
    print(f"Saved {out_dir / f'confusion_matrix_G60_cov_gt_{n}.png'}")


if __name__ == "__main__":
    main()
