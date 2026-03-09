#!/usr/bin/env python3
"""
Build a GIAB summary figure:
- Top row: 2 Precision/Recall bar plots (Full, 10x)
- Below each bar plot: 3 TP/FP/FN colored boxes aligned under the bars
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


METHOD_ORDER = ["ParaDISM", "BaseAligner"]
METRIC_ORDER = ["precision", "recall"]
COUNT_COLS = ["TP", "FP", "FN"]


def _load_metrics(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df


def _load_confusion(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df


def _get_metric(df: pd.DataFrame, method: str, metric: str) -> float:
    sub = df[df["method"] == method]
    if sub.empty:
        return 0.0
    return float(sub.iloc[0][metric])


def _get_counts(conf_df: pd.DataFrame, method: str) -> dict[str, int]:
    sub = conf_df[conf_df["method"] == method]
    if sub.empty:
        return {"TP": 0, "FP": 0, "FN": 0}
    r = sub.iloc[0]
    return {"TP": int(r["TP"]), "FP": int(r["FP"]), "FN": int(r["FN"])}


def draw_metric_panel(ax: plt.Axes, metrics_df: pd.DataFrame, title: str) -> None:
    x = np.arange(len(METRIC_ORDER), dtype=float)
    w = 0.34

    paradism_vals = [_get_metric(metrics_df, "ParaDISM", m) for m in METRIC_ORDER]
    base_vals = [_get_metric(metrics_df, "BaseAligner", m) for m in METRIC_ORDER]

    b1 = ax.bar(
        x - w / 2.0, paradism_vals, width=w,
        color="#D62728", edgecolor="black", linewidth=3.5,
    )
    b2 = ax.bar(
        x + w / 2.0, base_vals, width=w,
        color="#1F77B4", edgecolor="black", linewidth=3.5,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(["Precision", "Recall"], fontsize=34, fontweight="bold")
    ax.set_ylim(0.0, 1.08)
    ax.set_title(title, fontsize=38, fontweight="bold", pad=10)
    ax.grid(axis="y", alpha=0.28)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", length=0, labelsize=34)
    ax.tick_params(axis="y", labelsize=28)
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_color("#222222")

    for bars in (b1, b2):
        for bar in bars:
            v = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0, v + 0.015,
                f"{v:.3f}", ha="center", va="bottom",
                fontsize=26, fontweight="bold",
            )


def draw_count_boxes(ax: plt.Axes, para_counts: dict, base_counts: dict) -> None:
    """Draw 3 side-by-side TP/FP/FN boxes with ParaDISM (red) and Base (blue) values."""
    ax.set_xlim(-0.3, 4.2)
    ax.set_ylim(-0.1, 1.15)
    ax.axis("off")

    box_w = 1.25
    box_h = 0.75
    y0 = 0.25
    gap = 1.40

    for i, m in enumerate(COUNT_COLS):
        x0 = i * gap + 0.05
        box = FancyBboxPatch(
            (x0, y0), box_w, box_h,
            boxstyle="round,pad=0.03,rounding_size=0.04",
            facecolor="#f5f5f5",
            edgecolor="#222222",
            linewidth=3.0,
        )
        ax.add_patch(box)
        # ParaDISM value (red)
        ax.text(
            x0 + box_w / 2.0, y0 + 0.50,
            f"ParaDISM  {para_counts[m]}",
            ha="center", va="center",
            fontsize=30, fontweight="bold", color="#D62728",
        )
        # Base value (blue)
        ax.text(
            x0 + box_w / 2.0, y0 + 0.22,
            f"Base  {base_counts[m]}",
            ha="center", va="center",
            fontsize=30, fontweight="bold", color="#1F77B4",
        )
        # Label below box
        ax.text(
            x0 + box_w / 2.0, 0.08,
            m, ha="center", va="center",
            fontsize=30, fontweight="bold", color="#222222",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create GIAB Full vs 10x summary plot.")
    parser.add_argument("--full-dir", type=Path, default=Path("giab_benchmark/vcf_out_full"))
    parser.add_argument("--tenx-dir", type=Path, default=Path("giab_benchmark/vcf_out_10x"))
    parser.add_argument("--output", type=Path, default=Path("giab_benchmark/giab_full_10x_summary.svg"))
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    full_metrics = _load_metrics(args.full_dir / "balanced_comparison_metrics.csv")
    tenx_metrics = _load_metrics(args.tenx_dir / "balanced_comparison_metrics.csv")
    full_conf = _load_confusion(args.full_dir / "confusion_matrices_overall.csv")
    tenx_conf = _load_confusion(args.tenx_dir / "confusion_matrices_overall.csv")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })

    # Layout: 2 columns, bars on top, boxes below, compact
    fig = plt.figure(figsize=(26, 13))

    # Bar plots — upper portion
    # left bar: x=0.06..0.48, y=0.38..0.88
    # right bar: x=0.55..0.97, y=0.38..0.88
    ax_full = fig.add_axes([0.06, 0.36, 0.42, 0.54])
    ax_tenx = fig.add_axes([0.55, 0.36, 0.42, 0.54])

    # Box panels — directly below each bar plot, same x range
    ax_boxes_full = fig.add_axes([0.04, 0.03, 0.45, 0.28])
    ax_boxes_tenx = fig.add_axes([0.53, 0.03, 0.45, 0.28])

    draw_metric_panel(ax_full, full_metrics, "Full Coverage")
    draw_metric_panel(ax_tenx, tenx_metrics, "10x Coverage")

    # Legend spanning the top-left bar plot width
    legend_elements = [
        mpatches.Patch(facecolor="#D62728", edgecolor="black", linewidth=2.0, label="ParaDISM"),
        mpatches.Patch(facecolor="#1F77B4", edgecolor="black", linewidth=2.0, label="Base Aligner"),
    ]
    ax_full.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.16),
        frameon=True,
        ncol=2,
        fontsize=28,
        mode="expand",
        borderaxespad=0,
        edgecolor="#222222",
    )
    leg = ax_full.get_legend()
    leg.get_frame().set_linewidth(2.5)

    # Draw TP/FP/FN boxes
    para_full = _get_counts(full_conf, "ParaDISM")
    base_full = _get_counts(full_conf, "BaseAligner")
    para_tenx = _get_counts(tenx_conf, "ParaDISM")
    base_tenx = _get_counts(tenx_conf, "BaseAligner")

    draw_count_boxes(ax_boxes_full, para_full, base_full)
    draw_count_boxes(ax_boxes_tenx, para_tenx, base_tenx)

    # Panel labels
    panel_fs = 40
    panel_bbox = {"facecolor": "white", "edgecolor": "none", "alpha": 0.95, "pad": 0.15}
    for fx, fy, label in [(0.02, 0.96, "A"), (0.51, 0.96, "B")]:
        fig.text(fx, fy, label, fontsize=panel_fs, fontweight="bold",
                 ha="left", va="top", bbox=panel_bbox)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    png_path = args.output.with_suffix(".png")
    fig.savefig(png_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure: {args.output}")
    print(f"Saved figure: {png_path}")


if __name__ == "__main__":
    main()
