#!/usr/bin/env python3
"""Aggregate read-mapping CSV metrics across simulation seeds."""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support


def parse_read_mapping_summary_per_gene(csv_path: Path) -> dict:
    """Parse one per-seed read-mapping summary CSV into per-gene metrics."""
    df = pd.read_csv(csv_path)

    ground_truth = df["Ground_Truth"]
    if "Mapper_Prediction" in df.columns:
        mapper_pred = df["Mapper_Prediction"]
    elif "ParaDISM_Prediction" in df.columns:
        mapper_pred = df["ParaDISM_Prediction"]
    else:
        raise KeyError("Could not find Mapper_Prediction or ParaDISM_Prediction column")

    direct_pred_col = None
    for col in df.columns:
        if col.endswith("_Prediction") and col not in {
            "Mapper_Prediction",
            "ParaDISM_Prediction",
        }:
            direct_pred_col = col
            break
    direct_pred = df[direct_pred_col] if direct_pred_col else df.iloc[:, 3]

    labels = sorted(set(ground_truth) | set(mapper_pred) | set(direct_pred))

    return {
        "mapper": calculate_per_gene_metrics(ground_truth, mapper_pred, labels),
        "direct": calculate_per_gene_metrics(ground_truth, direct_pred, labels),
    }


def calculate_per_gene_metrics(y_true, y_pred, labels: list[str]) -> dict:
    """Calculate precision, recall, and specificity per gene plus weighted overall."""
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    metrics = {}

    for i, label in enumerate(labels):
        if label == "NONE":
            continue

        tp = cm[i, i]
        fp = np.sum(cm[:, i]) - tp
        fn = np.sum(cm[i, :]) - tp
        tn = np.sum(cm) - tp - fp - fn

        metrics[label] = {
            "precision": tp / (tp + fp) if (tp + fp) > 0 else 0.0,
            "recall": tp / (tp + fn) if (tp + fn) > 0 else 0.0,
            "specificity": tn / (tn + fp) if (tn + fp) > 0 else 0.0,
        }

    precision, recall, _, _ = precision_recall_fscore_support(
        y_true, y_pred, labels=labels, average="weighted", zero_division=0
    )

    specificities = []
    supports = []
    for i, _label in enumerate(labels):
        tn = np.sum(cm) - np.sum(cm[i, :]) - np.sum(cm[:, i]) + cm[i, i]
        fp = np.sum(cm[:, i]) - cm[i, i]
        specificities.append(tn / (tn + fp) if (tn + fp) > 0 else 0.0)
        supports.append(np.sum(cm[i, :]))

    total_support = np.sum(supports)
    specificity = (
        np.sum(np.array(specificities) * np.array(supports)) / total_support
        if total_support > 0
        else 0.0
    )

    metrics["Overall"] = {
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
    }
    return metrics


def aggregate_read_mapping_results(
    sim_output_base: Path,
    aligners: list[str],
    seed_start: int,
    seed_end: int,
) -> dict:
    collected = {
        aligner: {
            "mapper": defaultdict(lambda: defaultdict(list)),
            "direct": defaultdict(lambda: defaultdict(list)),
        }
        for aligner in aligners
    }

    for seed in range(seed_start, seed_end + 1):
        seed_dir = sim_output_base / f"seed_{seed}"
        for aligner in aligners:
            csv_path = (
                seed_dir
                / aligner
                / "read_mapping_analysis"
                / f"seed_{seed}_{aligner}_summary.csv"
            )
            if not csv_path.exists():
                continue

            try:
                metrics = parse_read_mapping_summary_per_gene(csv_path)
            except Exception as exc:
                print(f"Warning: failed to parse {csv_path}: {exc}", file=sys.stderr)
                continue

            for approach in ("mapper", "direct"):
                for gene, gene_metrics in metrics[approach].items():
                    for metric_name, value in gene_metrics.items():
                        collected[aligner][approach][gene][metric_name].append(value)

    aggregated = {}
    for aligner in aligners:
        aggregated[aligner] = {"mapper": {}, "direct": {}}
        for approach in ("mapper", "direct"):
            for gene, metrics in collected[aligner][approach].items():
                if not metrics["precision"]:
                    continue
                aggregated[aligner][approach][gene] = {
                    "precision_mean": np.mean(metrics["precision"]),
                    "precision_std": np.std(metrics["precision"]),
                    "recall_mean": np.mean(metrics["recall"]),
                    "recall_std": np.std(metrics["recall"]),
                    "specificity_mean": np.mean(metrics["specificity"]),
                    "specificity_std": np.std(metrics["specificity"]),
                    "n_seeds": len(metrics["precision"]),
                }

    return aggregated


def write_read_mapping_csv(aggregated: dict, aligners: list[str], output_dir: Path) -> Path | None:
    rows = []
    for aligner in aligners:
        if aligner not in aggregated:
            continue
        for approach in ("mapper", "direct"):
            for gene, metrics in aggregated[aligner][approach].items():
                rows.append(
                    {
                        "Aligner": aligner.upper(),
                        "Approach": approach.capitalize(),
                        "Gene": gene,
                        "Precision_Mean": metrics["precision_mean"],
                        "Precision_Std": metrics["precision_std"],
                        "Recall_Mean": metrics["recall_mean"],
                        "Recall_Std": metrics["recall_std"],
                        "Specificity_Mean": metrics["specificity_mean"],
                        "Specificity_Std": metrics["specificity_std"],
                        "N_Seeds": metrics["n_seeds"],
                    }
                )

    if not rows:
        return None

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "read_mapping_aggregated_summary.csv"
    pd.DataFrame(rows).to_csv(output_path, index=False)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate simulation read-mapping metrics across seeds."
    )
    parser.add_argument(
        "--sim-output-base",
        default="sim_output",
        help="Base directory for simulation outputs",
    )
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-end", type=int, default=1000)
    parser.add_argument(
        "--aligners",
        nargs="+",
        default=["bowtie2", "bwa-mem2", "minimap2"],
        help="Aligners to aggregate",
    )
    parser.add_argument(
        "--output-dir",
        default="aggregated_results",
        help="Output directory for aggregated CSV results",
    )
    args = parser.parse_args()

    aggregated = aggregate_read_mapping_results(
        Path(args.sim_output_base), args.aligners, args.seed_start, args.seed_end
    )
    output_path = write_read_mapping_csv(aggregated, args.aligners, Path(args.output_dir))

    if output_path is None:
        print("No read-mapping summaries found; nothing to aggregate.", file=sys.stderr)
        return

    n_aligners = len([aligner for aligner in args.aligners if aligner in aggregated])
    print(f"Saved: {output_path}")
    print(f"Aggregated read-mapping metrics for {n_aligners} aligners.")


if __name__ == "__main__":
    main()
