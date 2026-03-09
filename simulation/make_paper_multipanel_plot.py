#!/usr/bin/env python3
"""
Create a 2x2 summary figure from existing simulation CSV outputs.

Layout:
  [0,0] Three grouped bar plots (Overall, PKD1, Pseudogene Avg; red+blue in each metric group)
  [0,1] Averaged iteration progression plot
  [1,0] Averaged ParaDISM confusion matrix
  [1,1] Averaged Direct confusion matrix
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import lines as mlines
from matplotlib import patches as mpatches
from matplotlib.colors import LinearSegmentedColormap


ALIGNER_ORDER = ["bwa-mem2", "bowtie2", "minimap2"]
GENE_ORDER = ["PKD1", "PKD1P1", "PKD1P2", "PKD1P3", "PKD1P4", "PKD1P5", "PKD1P6", "NONE"]
METRIC_ORDER = ["Precision", "Recall", "Specificity"]
PSEUDOGENE_ORDER = ["PKD1P1", "PKD1P2", "PKD1P3", "PKD1P4", "PKD1P5", "PKD1P6"]
PARA_COLOR = "#B7A1E6"
DIRECT_COLOR = "#F6C08B"


def configure_style(font_scale: float) -> None:
    base = 20.0 * font_scale
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans Mono",
            "font.size": base,
            "axes.titlesize": base + 7,
            "axes.labelsize": base + 4,
            "xtick.labelsize": base + 1,
            "ytick.labelsize": base + 1,
            "legend.fontsize": base,
            "font.weight": "bold",
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def _combine_mean_std(mean_stack: np.ndarray, std_stack: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean_avg = mean_stack.mean(axis=0)
    var_avg = np.maximum((std_stack**2 + mean_stack**2).mean(axis=0) - mean_avg**2, 0.0)
    return mean_avg, np.sqrt(var_avg)


def _style_axis(ax: plt.Axes) -> None:
    ax.grid(axis="y", alpha=0.30)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(2.0)
    ax.spines["left"].set_color("#222222")
    ax.spines["bottom"].set_linewidth(2.0)
    ax.spines["bottom"].set_color("#222222")


def _aggregate_rows_for_metric(rows: pd.DataFrame, metric: str) -> tuple[float, float] | None:
    if rows.empty:
        return None
    vals = rows[f"{metric}_Mean"].to_numpy(dtype=float).reshape(-1, 1)
    stds = rows[f"{metric}_Std"].to_numpy(dtype=float).reshape(-1, 1)
    m, s = _combine_mean_std(vals, stds)
    return float(m[0]), float(s[0])


def build_avg_bar_panel(ax: plt.Axes, panel_b_all_aligners_csv: Path, gene: str | list[str], title: str) -> None:
    """
    Build one grouped bar panel with 3 metric groups:
      Precision, Recall, Specificity
    and two bars per group:
      ParaDISM (red), Direct (blue), each averaged across aligners.
    """
    df = pd.read_csv(panel_b_all_aligners_csv)
    if isinstance(gene, list):
        df = df[df["Gene"].isin(gene)].copy()
    else:
        df = df[df["Gene"] == gene].copy()

    para_df = df[df["Approach"] == "ParaDISM"].copy()
    direct_df = df[df["Approach"] == "Direct"].copy()

    para_means = []
    para_stds = []
    direct_means = []
    direct_stds = []

    for metric in METRIC_ORDER:
        p_vals = []
        p_stds = []
        d_vals = []
        d_stds = []
        for aligner in ALIGNER_ORDER:
            p_rows = para_df[para_df["Aligner"] == aligner]
            p_agg = _aggregate_rows_for_metric(p_rows, metric)
            if p_agg is not None:
                p_vals.append(p_agg[0])
                p_stds.append(p_agg[1])

            d_rows = direct_df[direct_df["Aligner"] == aligner]
            d_agg = _aggregate_rows_for_metric(d_rows, metric)
            if d_agg is not None:
                d_vals.append(d_agg[0])
                d_stds.append(d_agg[1])

        if p_vals:
            m, s = _combine_mean_std(np.array(p_vals).reshape(-1, 1), np.array(p_stds).reshape(-1, 1))
            para_means.append(float(m[0]))
            para_stds.append(float(s[0]))
        else:
            para_means.append(0.0)
            para_stds.append(0.0)

        if d_vals:
            m, s = _combine_mean_std(np.array(d_vals).reshape(-1, 1), np.array(d_stds).reshape(-1, 1))
            direct_means.append(float(m[0]))
            direct_stds.append(float(s[0]))
        else:
            direct_means.append(0.0)
            direct_stds.append(0.0)

    x = np.arange(len(METRIC_ORDER), dtype=float)
    # Wider bars to reduce spacing between metric groups while keeping
    # the gap between the two bars inside each group unchanged.
    w = 0.36
    inner_gap = 0.02

    para_bars = ax.bar(
        x - w / 2.0 - inner_gap,
        para_means,
        yerr=para_stds,
        width=w,
        color=PARA_COLOR,
        capsize=8,
        edgecolor="black",
        linewidth=7.0,
        error_kw={"elinewidth": 2.0, "ecolor": "#111111", "capthick": 2.0},
    )
    direct_bars = ax.bar(
        x + w / 2.0 + inner_gap,
        direct_means,
        yerr=direct_stds,
        width=w,
        color=DIRECT_COLOR,
        capsize=8,
        edgecolor="black",
        linewidth=7.0,
        error_kw={"elinewidth": 2.0, "ecolor": "#111111", "capthick": 2.0},
    )

    # Match bar-panel metric label size to iteration-panel text sizing.
    metric_label_size = max(30, int(plt.rcParams["font.size"] * 1.05))
    tick_size = max(34, int(plt.rcParams["font.size"] * 1.25))
    title_size = tick_size + 4
    ax.set_ylim(0.0, 1.08)
    ax.set_xticks(x)
    ax.set_xticklabels(["Precision", "Recall", "Specificity"], fontsize=metric_label_size, fontweight="bold")
    ax.set_title("")
    ax.tick_params(axis="x", length=0, labelsize=metric_label_size, pad=12)
    ax.tick_params(axis="y", labelsize=tick_size)
    for tick in ax.get_yticklabels():
        tick.set_fontweight("normal")
    _style_axis(ax)

    for bar, v, e in zip(para_bars, para_means, para_stds):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            v + e + 0.012,
            f"{v:.3f}",
            ha="center",
            va="bottom",
            fontsize=max(24, int(plt.rcParams["font.size"] * 0.95)),
            fontweight="normal",
        )
    for bar, v, e in zip(direct_bars, direct_means, direct_stds):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            v + e + 0.012,
            f"{v:.3f}",
            ha="center",
            va="bottom",
            fontsize=max(24, int(plt.rcParams["font.size"] * 0.95)),
            fontweight="normal",
        )


def _matrix_for_aligner(df: pd.DataFrame, aligner: str, approach_group: str) -> tuple[np.ndarray, np.ndarray]:
    if approach_group == "paradism":
        sub = df[(df["Aligner"] == aligner) & (df["Approach"] == "ParaDISM")]
    elif approach_group == "direct":
        sub = df[(df["Aligner"] == aligner) & (df["Approach"].astype(str).str.startswith("Direct"))]
    else:
        raise ValueError(f"Unknown approach_group: {approach_group}")

    if sub.empty:
        raise ValueError(f"No confusion rows for aligner={aligner}, approach={approach_group}")

    mean_m = (
        sub.pivot(index="Origin_Gene", columns="Mapped_Gene", values="Mean_Reads")
        .reindex(index=GENE_ORDER, columns=GENE_ORDER)
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    std_m = (
        sub.pivot(index="Origin_Gene", columns="Mapped_Gene", values="Std_Reads")
        .reindex(index=GENE_ORDER, columns=GENE_ORDER)
        .fillna(0.0)
        .to_numpy(dtype=float)
    )
    return mean_m, std_m


def load_confusion_averages(panel_a_all_aligners_csv: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(panel_a_all_aligners_csv)
    para_means = []
    para_stds = []
    direct_means = []
    direct_stds = []
    for aligner in ALIGNER_ORDER:
        m, s = _matrix_for_aligner(df, aligner, "paradism")
        para_means.append(m)
        para_stds.append(s)
        m, s = _matrix_for_aligner(df, aligner, "direct")
        direct_means.append(m)
        direct_stds.append(s)

    para_mean, para_std = _combine_mean_std(np.stack(para_means, axis=0), np.stack(para_stds, axis=0))
    direct_mean, direct_std = _combine_mean_std(np.stack(direct_means, axis=0), np.stack(direct_stds, axis=0))
    return para_mean, para_std, direct_mean, direct_std


def draw_confusion(ax: plt.Axes, mean_m: np.ndarray, std_m: np.ndarray, vmax: float, cmap=None) -> plt.Axes:
    xlabels = ["P", "P1", "P2", "P3", "P4", "P5", "P6", "None"]
    ylabels = [g if g != "NONE" else "None" for g in GENE_ORDER]
    if cmap is None:
        cmap = LinearSegmentedColormap.from_list(
            "orange_brown",
            ["#fffdf8", "#fef3df", "#f9d08f", "#f3a44a", "#ea7f1d", "#d8640d"],
        )
    # Use auto aspect so cells can be slightly wider than tall.
    im = ax.imshow(mean_m, cmap=cmap, vmin=0.0, vmax=vmax, aspect="auto")
    ax.set_facecolor("white")
    idx = np.arange(len(GENE_ORDER))
    ax.set_xticks(idx)
    ax.set_yticks(idx)
    tick_size = max(24, int(plt.rcParams["font.size"] * 0.90))
    ax.set_xticklabels(xlabels, rotation=0, ha="center", fontsize=tick_size)
    ax.set_yticklabels(ylabels, fontsize=tick_size)
    ax.tick_params(length=0)
    ax.tick_params(axis="y", pad=10)
    ax.set_xticks(np.arange(-0.5, len(GENE_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(GENE_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0, alpha=0.75)
    ax.tick_params(which="minor", bottom=False, left=False)

    txt_size = max(30, int(plt.rcParams["font.size"] * 1.15))
    for i in range(mean_m.shape[0]):
        for j in range(mean_m.shape[1]):
            ax.text(
                j,
                i,
                f"{mean_m[i, j]:.0f}\n±{std_m[i, j]:.0f}",
                ha="center",
                va="center",
                color="#111111",
                fontsize=txt_size,
                fontweight="normal",
                linespacing=0.85,
            )

    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_color("#222222")
    return im


def _load_iteration_frames(iteration_plots_dir: Path) -> dict[str, pd.DataFrame]:
    frames: dict[str, pd.DataFrame] = {}
    for aligner in ALIGNER_ORDER:
        p = iteration_plots_dir / f"iteration_progression_{aligner}.csv"
        if p.exists():
            frames[aligner] = pd.read_csv(p).sort_values("Iteration").reset_index(drop=True)
    if not frames:
        raise FileNotFoundError(f"No iteration progression CSVs found in {iteration_plots_dir}")
    return frames


def _common_iterations(frames: dict[str, pd.DataFrame]) -> np.ndarray:
    common: set[int] | None = None
    for df in frames.values():
        vals = set(df["Iteration"].astype(int).tolist())
        common = vals if common is None else (common & vals)
    if not common:
        raise ValueError("No common iterations across aligners")
    return np.array(sorted(common), dtype=int)


def _avg_metric(frames: dict[str, pd.DataFrame], iterations: np.ndarray, mean_col: str, std_col: str) -> tuple[np.ndarray, np.ndarray]:
    mean_rows = []
    std_rows = []
    for aligner in ALIGNER_ORDER:
        if aligner not in frames:
            continue
        sub = frames[aligner]
        sub = sub[sub["Iteration"].isin(iterations)].sort_values("Iteration")
        mean_rows.append(sub[mean_col].to_numpy(dtype=float))
        std_rows.append(sub[std_col].to_numpy(dtype=float))
    return _combine_mean_std(np.stack(mean_rows, axis=0), np.stack(std_rows, axis=0))


def draw_iteration_panels(axes: list[plt.Axes], iteration_plots_dir: Path, base_fontsize: float) -> None:
    frames = _load_iteration_frames(iteration_plots_dir)
    iterations = _common_iterations(frames)
    x = iterations.astype(float)

    metric_colors = {
        "Precision": (PARA_COLOR, DIRECT_COLOR),
        "Recall": (PARA_COLOR, DIRECT_COLOR),
        "Specificity": (PARA_COLOR, DIRECT_COLOR),
    }
    metric_style = {
        "Precision": {"linestyle": "-", "marker": "o"},
        "Recall": {"linestyle": "-", "marker": "s"},
        "Specificity": {"linestyle": "-", "marker": "^"},
    }
    for i, (ax, metric) in enumerate(zip(axes, METRIC_ORDER)):
        para_mean, para_std = _avg_metric(
            frames,
            iterations,
            mean_col=f"ParaDISM_{metric}_Mean",
            std_col=f"ParaDISM_{metric}_Std",
        )
        direct_mean, direct_std = _avg_metric(
            frames,
            iterations,
            mean_col=f"Direct_{metric}_Mean",
            std_col=f"Direct_{metric}_Std",
        )

        st = metric_style[metric]
        para_color, direct_color = metric_colors[metric]
        ax.plot(
            x,
            para_mean,
            color=para_color,
            linestyle="-",
            marker=st["marker"],
            markersize=20,
            markeredgewidth=2.2,
            linewidth=8.0,
            label=f"ParaDISM {metric}",
        )
        ax.fill_between(x, para_mean - para_std, para_mean + para_std, color=para_color, alpha=0.35)
        ax.plot(
            x,
            direct_mean,
            color=direct_color,
            linestyle="-",
            marker=st["marker"],
            markersize=20,
            markeredgewidth=2.2,
            linewidth=8.0,
            label=f"Direct {metric}",
        )
        ax.fill_between(x, direct_mean - direct_std, direct_mean + direct_std, color=direct_color, alpha=0.30)

        ax.set_xlim(float(x.min()) - 0.25, float(x.max()) + 0.25)
        if metric == "Precision":
            y0, y1 = 0.70, 1.00
        elif metric == "Recall":
            y0, y1 = 0.40, 1.00
        else:  # Specificity
            y0, y1 = 0.90, 1.00
        ax.set_ylim(y0, y1)
        # Keep only start/end y ticks (no intermediate values)
        ax.set_yticks([y0, y1])
        ax.grid(axis="both", alpha=0.30)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_linewidth(2.0)
            spine.set_color("#222222")

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        iter_label_fs = max(28, int(base_fontsize * 1.00))
        iter_tick_fs = max(24, int(base_fontsize * 0.95))
        ax.set_title("")
        ax.set_ylabel(metric, fontsize=iter_label_fs, fontweight="bold", rotation=90, labelpad=10)
        ax.tick_params(axis="both", labelsize=iter_tick_fs)
        for tick in ax.get_yticklabels() + ax.get_xticklabels():
            tick.set_fontweight("normal")
        if i < 2:
            ax.set_xticklabels([])
            ax.set_xlabel("")
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="x", length=0)
        else:
            ax.set_xlabel("Iteration", fontsize=iter_label_fs, fontweight="bold", labelpad=6)
        ax.yaxis.set_label_position("left")


def build_figure(
    panel_a_all_aligners_csv: Path,
    panel_b_all_aligners_csv: Path,
    iteration_plots_dir: Path,
    output_path: Path,
    dpi: int,
    font_scale: float,
) -> None:
    configure_style(font_scale=font_scale)
    fig = plt.figure(figsize=(36, 34))
    # Leave room at top for shared legend
    gs = fig.add_gridspec(2, 2, height_ratios=[0.80, 1.30], wspace=0.24, hspace=0.14, top=0.94, left=0.14, right=0.98)

    # Top-left [0,0]: Overall bar chart
    ax_bar_overall = fig.add_subplot(gs[0, 0])

    # Top-right [0,1]: Iteration progression (3 stacked subpanels)
    gs_iter = gs[0, 1].subgridspec(3, 1, hspace=0.26)
    ax_iter = [fig.add_subplot(gs_iter[i, 0]) for i in range(3)]

    # Bottom row: confusion matrices in a tight subgrid so they sit closer together.
    gs_cm = gs[1, :].subgridspec(1, 2, wspace=0.03)

    # Bottom-left [1,0]: ParaDISM confusion matrix
    ax_cm_para = fig.add_subplot(gs_cm[0, 0])
    ax_cm_para.set_box_aspect(0.95)

    # Bottom-right [1,1]: Direct confusion matrix
    ax_cm_direct = fig.add_subplot(gs_cm[0, 1])
    ax_cm_direct.set_box_aspect(0.95)

    # Panel A: Overall bar chart only
    build_avg_bar_panel(ax_bar_overall, panel_b_all_aligners_csv, gene="Overall", title="Overall")
    bar_xtick_size = (
        ax_bar_overall.get_xticklabels()[0].get_fontsize()
        if ax_bar_overall.get_xticklabels()
        else plt.rcParams["xtick.labelsize"]
    )

    # Shared legend at top of figure, spanning full width
    legend_elements = [
        mpatches.Patch(facecolor=PARA_COLOR, edgecolor="black", linewidth=7.0, label="ParaDISM"),
        mpatches.Patch(facecolor=DIRECT_COLOR, edgecolor="black", linewidth=7.0, label="Direct Alignment"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.56, 0.99),
        frameon=False,
        ncol=2,
        fontsize=max(24, int(bar_xtick_size * 0.9)),
        handlelength=1.25,
        handleheight=1.25,
        columnspacing=4.0,
    )

    # Panel B: Iteration progression — full quadrant width
    draw_iteration_panels(ax_iter, iteration_plots_dir, base_fontsize=bar_xtick_size)

    # Panel C & D: Confusion matrices
    para_mean, para_std, direct_mean, direct_std = load_confusion_averages(panel_a_all_aligners_csv)
    vmax = max(float(para_mean.max()), float(direct_mean.max()))
    cmap_purple = LinearSegmentedColormap.from_list(
        "purple_cmap", ["#fbf9ff", "#efe7fb", "#d9c7f3", "#b79de4", PARA_COLOR, "#7f67c2"],
    )
    cmap_orange = LinearSegmentedColormap.from_list(
        "orange_cmap", ["#fffaf4", "#fdebd8", "#f9d4b0", "#f3b77f", DIRECT_COLOR, "#d78f4e"],
    )
    im = draw_confusion(ax_cm_para, para_mean, para_std, vmax=vmax, cmap=cmap_purple)
    _ = draw_confusion(ax_cm_direct, direct_mean, direct_std, vmax=vmax, cmap=cmap_orange)
    # Show y labels on left matrix, hide on right matrix.
    ax_cm_para.tick_params(axis="y", labelleft=True)
    ax_cm_direct.tick_params(axis="y", labelleft=False)

    # Panel labels centered above each panel.
    panel_fs = bar_xtick_size + 10
    panel_axes = [
        ("A", ax_bar_overall),
        ("B", ax_iter[0]),
        ("C", ax_cm_para),
        ("D", ax_cm_direct),
    ]
    label_x_shift = -0.015
    label_y_base = 0.006
    label_extra_y = {"A": -0.008, "B": -0.014}
    for label, pax in panel_axes:
        pos = pax.get_position()
        cx = pos.x0 + pos.width / 2.0 + label_x_shift
        cy = pos.y1 + label_y_base + label_extra_y.get(label, 0.0)
        fig.text(
            cx,
            cy,
            label,
            fontsize=panel_fs,
            fontweight="bold",
            ha="center",
            va="bottom",
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight")
    png_path = output_path.with_suffix(".png")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved PNG: {png_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create 2x2 multi-panel plot from existing aggregated CSVs.")
    parser.add_argument(
        "--figure1-data-dir",
        type=Path,
        default=Path("simulation/simulations_outputs_1000/figure1_data"),
        help="Directory containing panel_A_all_aligners.csv and panel_B_all_aligners.csv",
    )
    parser.add_argument(
        "--iteration-plots-dir",
        type=Path,
        default=Path("simulation/simulations_outputs_1000/iteration_plots"),
        help="Directory containing iteration_progression_<aligner>.csv files",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("simulation/simulations_outputs_1000/figure1_data/paper_multipanel_figure.svg"),
        help="Output image path",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Output DPI")
    parser.add_argument("--font-scale", type=float, default=1.9, help="Global font scale multiplier")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    panel_a = args.figure1_data_dir / "panel_A_all_aligners.csv"
    panel_b = args.figure1_data_dir / "panel_B_all_aligners.csv"
    if not panel_a.exists():
        raise FileNotFoundError(f"Missing required CSV: {panel_a}")
    if not panel_b.exists():
        raise FileNotFoundError(f"Missing required CSV: {panel_b}")
    if not args.iteration_plots_dir.exists():
        raise FileNotFoundError(f"Missing iteration plots directory: {args.iteration_plots_dir}")

    build_figure(
        panel_a_all_aligners_csv=panel_a,
        panel_b_all_aligners_csv=panel_b,
        iteration_plots_dir=args.iteration_plots_dir,
        output_path=args.output,
        dpi=args.dpi,
        font_scale=args.font_scale,
    )
    print(f"Saved multi-panel figure to: {args.output}")


if __name__ == "__main__":
    main()
