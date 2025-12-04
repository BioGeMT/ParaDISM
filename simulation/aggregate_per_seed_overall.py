#!/usr/bin/env python3
"""
Aggregate per-seed overall metrics for ParaDISM (mapper) and direct aligner.
Reads per-seed summary CSVs from read_mapping_analysis and writes a combined CSV.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support, confusion_matrix


def calculate_weighted_specificity(y_true, y_pred, labels):
    """Mirror ParaDISM specificity: weighted by class support."""
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    specs = []
    supports = []
    for i in range(len(labels)):
        tn = np.sum(cm) - np.sum(cm[i, :]) - np.sum(cm[:, i]) + cm[i, i]
        fp = np.sum(cm[:, i]) - cm[i, i]
        spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
        specs.append(spec)
        supports.append(np.sum(cm[i, :]))
    total = np.sum(supports)
    return (np.sum(np.array(specs) * np.array(supports)) / total) if total > 0 else 0.0


def extract_overall(summary_path: Path, aligner: str) -> dict | None:
    if not summary_path.exists():
        return None
    df = pd.read_csv(summary_path)
    if df.empty:
        return None

    y_true = df["Ground_Truth"]
    if "Mapper_Prediction" in df.columns:
        y_mapper = df["Mapper_Prediction"]
    else:
        y_mapper = df["ParaDISM_Prediction"]

    direct_col = None
    for col in df.columns:
        if col.endswith("_Prediction") and col not in ("Mapper_Prediction", "ParaDISM_Prediction"):
            direct_col = col
            break
    if direct_col is None:
        direct_col = df.columns[3]
    y_direct = df[direct_col]

    labels = sorted(set(y_true) | set(y_mapper) | set(y_direct))

    prec_map, rec_map, _, _ = precision_recall_fscore_support(
        y_true, y_mapper, labels=labels, average="weighted", zero_division=0
    )
    spec_map = calculate_weighted_specificity(y_true, y_mapper, labels)

    prec_dir, rec_dir, _, _ = precision_recall_fscore_support(
        y_true, y_direct, labels=labels, average="weighted", zero_division=0
    )
    spec_dir = calculate_weighted_specificity(y_true, y_direct, labels)

    return {
        "ParaDISM_Precision": prec_map,
        "ParaDISM_Recall": rec_map,
        "ParaDISM_Specificity": spec_map,
        "Reads": len(df),
        f"{aligner.upper()}_Precision": prec_dir,
        f"{aligner.upper()}_Recall": rec_dir,
        f"{aligner.upper()}_Specificity": spec_dir,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-seed overall metrics across aligners."
    )
    parser.add_argument("--sim-output-base", required=True, help="Simulation output base")
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-end", type=int, default=1000)
    parser.add_argument(
        "--aligners",
        nargs="+",
        default=["bwa-mem2", "bowtie2", "minimap2"],
        help="Aligners to include",
    )
    parser.add_argument(
        "--output",
        default="per_seed_overall_metrics.csv",
        help="Output CSV path",
    )
    args = parser.parse_args()

    rows = []
    for seed in range(args.seed_start, args.seed_end + 1):
        for aligner in args.aligners:
            summary = (
                Path(args.sim_output_base)
                / f"seed_{seed}"
                / aligner
                / "read_mapping_analysis"
                / f"seed_{seed}_{aligner}_summary.csv"
            )
            metrics = extract_overall(summary, aligner)
            if metrics is None:
                continue
            row = {"Seed": seed, "Aligner": aligner}
            row.update(metrics)
            rows.append(row)

    if not rows:
        print("No summaries found; nothing to write.", file=sys.stderr)
        return

    df = pd.DataFrame(rows)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"Wrote {len(df)} rows to {out_path}")


if __name__ == "__main__":
    main()
