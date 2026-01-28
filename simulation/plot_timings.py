#!/usr/bin/env python3
"""
Plot wall-clock timing per aligner with error bars.
Reads timing_data.csv (seed, aligner, wall_clock_seconds) and writes a bar plot and summary CSV.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Plot timing by aligner with error bars.")
    parser.add_argument(
        "--timing-csv",
        required=True,
        help="Path to timing_data.csv",
    )
    parser.add_argument(
        "--output-dir",
        default="aggregated_results",
        help="Directory to write timing plot/summary",
    )
    parser.add_argument(
        "--unit",
        choices=["seconds", "minutes"],
        default="minutes",
        help="Y-axis unit for the plot",
    )
    args = parser.parse_args()

    timing_path = Path(args.timing_csv)
    if not timing_path.exists():
        raise FileNotFoundError(f"Timing CSV not found: {timing_path}")

    df = pd.read_csv(timing_path)
    if df.empty:
        raise ValueError("Timing CSV is empty.")

    grouped = (
        df.groupby("aligner")["wall_clock_seconds"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    if grouped.empty:
        raise ValueError("No timing data after grouping.")

    scale = 60.0 if args.unit == "minutes" else 1.0
    grouped["mean_scaled"] = grouped["mean"] / scale
    grouped["std_scaled"] = grouped["std"].fillna(0.0) / scale

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save summary CSV
    summary_csv = output_dir / "timing_by_aligner.csv"
    grouped.rename(
        columns={"mean": "mean_seconds", "std": "std_seconds", "count": "n"},
        inplace=True,
    )
    grouped.to_csv(summary_csv, index=False)

    # Plot
    plt.style.use("seaborn-v0_8")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(
        grouped["aligner"],
        grouped["mean_scaled"],
        yerr=grouped["std_scaled"],
        capsize=6,
        color="#4c78a8",
        alpha=0.85,
    )
    ax.set_ylabel(f"Wall-clock time ({args.unit})")
    ax.set_xlabel("Aligner")
    ax.set_title("ParaDISM wall-clock time by aligner")
    ax.margins(y=0.1)
    plt.tight_layout()

    plot_path = output_dir / "timing_by_aligner.png"
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)

    print(f"Wrote timing summary: {summary_csv}")
    print(f"Wrote timing plot:    {plot_path}")


if __name__ == "__main__":
    main()
