#!/usr/bin/env python3
"""
Build CSV data for Figure 1 panels A and B.

Panel A: Confusion matrices (avg ± std over 1000 seeds) for ParaDISM (final iteration).
Panel B: Per-gene precision, recall, specificity (avg ± std over 1000 seeds).
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

GENES = ["PKD1", "PKD1P1", "PKD1P2", "PKD1P3", "PKD1P4", "PKD1P5", "PKD1P6"]
ALL_LABELS = GENES + ["NONE"]
LABEL_TO_IDX = {l: i for i, l in enumerate(ALL_LABELS)}
N = len(ALL_LABELS)


def build_cm_fast(csv_path):
    """Build confusion matrices from a per-read summary CSV using vectorized ops."""
    df = pd.read_csv(csv_path, usecols=lambda c: c in {
        "Ground_Truth", "ParaDISM_Prediction"
    } or c.endswith("_Prediction"))

    # Find direct aligner prediction column
    direct_col = None
    for col in df.columns:
        if col.endswith("_Prediction") and col != "ParaDISM_Prediction":
            direct_col = col
            break

    # Map labels to indices
    gt_idx = df["Ground_Truth"].map(LABEL_TO_IDX)
    para_idx = df["ParaDISM_Prediction"].map(LABEL_TO_IDX)

    # Drop rows with unknown labels
    valid = gt_idx.notna() & para_idx.notna()
    gt_idx = gt_idx[valid].astype(int).values
    para_idx = para_idx[valid].astype(int).values

    # Build ParaDISM CM
    cm = np.zeros((N, N), dtype=float)
    np.add.at(cm, (gt_idx, para_idx), 1)

    # Build direct aligner CM
    cm_d = None
    if direct_col:
        d_idx = df[direct_col].map(LABEL_TO_IDX)
        valid_d = d_idx.notna()
        # Use same valid mask intersected
        mask = valid & valid_d
        gt_d = df["Ground_Truth"].map(LABEL_TO_IDX)[mask].astype(int).values
        d_vals = d_idx[mask].astype(int).values
        cm_d = np.zeros((N, N), dtype=float)
        np.add.at(cm_d, (gt_d, d_vals), 1)

    return cm, cm_d


def metrics_from_cm(cm, n_genes=7):
    """Compute per-gene precision, recall, specificity from confusion matrix."""
    prec = np.zeros(n_genes)
    rec = np.zeros(n_genes)
    spec = np.zeros(n_genes)
    for i in range(n_genes):
        tp = cm[i, i]
        fp = cm[:, i].sum() - tp
        fn = cm[i, :].sum() - tp
        tn = cm.sum() - tp - fp - fn
        prec[i] = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec[i] = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        spec[i] = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    return prec, rec, spec


def main():
    parser = argparse.ArgumentParser(description="Build Figure 1 CSV data")
    parser.add_argument("--sim-output-base", default="simulations_outputs_1000")
    parser.add_argument("--aligner", default="bowtie2")
    parser.add_argument("--seed-start", type=int, default=1)
    parser.add_argument("--seed-end", type=int, default=1000)
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    base = Path(args.sim_output_base)
    out_dir = Path(args.output_dir) if args.output_dir else base / "figure1_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    aligner = args.aligner

    all_cm = []
    all_cm_d = []
    all_prec, all_rec, all_spec = [], [], []
    all_prec_d, all_rec_d, all_spec_d = [], [], []

    for seed in range(args.seed_start, args.seed_end + 1):
        seed_dir = base / f"seed_{seed}" / aligner / "read_mapping_analysis"
        # Prefer final iteration (iter10), fall back to base
        csv_path = seed_dir / f"seed_{seed}_{aligner}_iter10_summary.csv"
        if not csv_path.exists():
            csv_path = seed_dir / f"seed_{seed}_{aligner}_summary.csv"
        if not csv_path.exists():
            continue

        cm, cm_d = build_cm_fast(csv_path)
        all_cm.append(cm)

        p, r, s = metrics_from_cm(cm)
        all_prec.append(p)
        all_rec.append(r)
        all_spec.append(s)

        if cm_d is not None:
            all_cm_d.append(cm_d)
            pd_, rd_, sd_ = metrics_from_cm(cm_d)
            all_prec_d.append(pd_)
            all_rec_d.append(rd_)
            all_spec_d.append(sd_)

        if seed % 100 == 0:
            print(f"  Processed seed {seed}...")

    n_seeds = len(all_cm)
    print(f"Processed {n_seeds} seeds")

    # ===== Panel A: Confusion matrices =====
    cm_stack = np.stack(all_cm)
    cm_mean = cm_stack.mean(axis=0)
    cm_std = cm_stack.std(axis=0)

    rows_a = []
    for i, gt in enumerate(ALL_LABELS):
        for j, pred in enumerate(ALL_LABELS):
            rows_a.append({
                "Approach": "ParaDISM",
                "Origin_Gene": gt,
                "Mapped_Gene": pred,
                "Mean_Reads": round(cm_mean[i, j], 1),
                "Std_Reads": round(cm_std[i, j], 1),
            })

    if all_cm_d:
        cm_d_stack = np.stack(all_cm_d)
        cm_d_mean = cm_d_stack.mean(axis=0)
        cm_d_std = cm_d_stack.std(axis=0)
        for i, gt in enumerate(ALL_LABELS):
            for j, pred in enumerate(ALL_LABELS):
                rows_a.append({
                    "Approach": f"Direct ({aligner})",
                    "Origin_Gene": gt,
                    "Mapped_Gene": pred,
                    "Mean_Reads": round(cm_d_mean[i, j], 1),
                    "Std_Reads": round(cm_d_std[i, j], 1),
                })

    pd.DataFrame(rows_a).to_csv(out_dir / "panel_A_confusion_matrix.csv", index=False)

    # Compact matrix format
    with open(out_dir / "panel_A_confusion_matrix_compact.csv", "w") as f:
        f.write(f"ParaDISM Confusion Matrix (Mean +/- Std over {n_seeds} seeds)\n")
        f.write("Rows=Origin Cols=Mapped\n")
        f.write("," + ",".join(ALL_LABELS) + "\n")
        for i, gt in enumerate(ALL_LABELS):
            vals = [f"{cm_mean[i,j]:.1f} +/- {cm_std[i,j]:.1f}" for j in range(N)]
            f.write(gt + "," + ",".join(vals) + "\n")

        if all_cm_d:
            f.write(f"\nDirect {aligner} Confusion Matrix (Mean +/- Std over {n_seeds} seeds)\n")
            f.write("Rows=Origin Cols=Mapped\n")
            f.write("," + ",".join(ALL_LABELS) + "\n")
            for i, gt in enumerate(ALL_LABELS):
                vals = [f"{cm_d_mean[i,j]:.1f} +/- {cm_d_std[i,j]:.1f}" for j in range(N)]
                f.write(gt + "," + ",".join(vals) + "\n")

    print(f"Panel A: {out_dir / 'panel_A_confusion_matrix.csv'}")
    print(f"Panel A compact: {out_dir / 'panel_A_confusion_matrix_compact.csv'}")

    # ===== Panel B: Per-gene metrics =====
    prec_arr = np.stack(all_prec)
    rec_arr = np.stack(all_rec)
    spec_arr = np.stack(all_spec)

    support_mean = cm_mean[:len(GENES), :].sum(axis=1)

    rows_b = []
    for i, gene in enumerate(GENES):
        rows_b.append({
            "Approach": "ParaDISM", "Gene": gene,
            "Precision_Mean": round(prec_arr[:, i].mean(), 4),
            "Precision_Std": round(prec_arr[:, i].std(), 4),
            "Recall_Mean": round(rec_arr[:, i].mean(), 4),
            "Recall_Std": round(rec_arr[:, i].std(), 4),
            "Specificity_Mean": round(spec_arr[:, i].mean(), 4),
            "Specificity_Std": round(spec_arr[:, i].std(), 4),
        })
    rows_b.append({
        "Approach": "ParaDISM", "Gene": "Overall",
        "Precision_Mean": round(np.average(prec_arr.mean(axis=0), weights=support_mean), 4),
        "Precision_Std": round(np.average(prec_arr.std(axis=0), weights=support_mean), 4),
        "Recall_Mean": round(np.average(rec_arr.mean(axis=0), weights=support_mean), 4),
        "Recall_Std": round(np.average(rec_arr.std(axis=0), weights=support_mean), 4),
        "Specificity_Mean": round(np.average(spec_arr.mean(axis=0), weights=support_mean), 4),
        "Specificity_Std": round(np.average(spec_arr.std(axis=0), weights=support_mean), 4),
    })

    if all_prec_d:
        prec_d_arr = np.stack(all_prec_d)
        rec_d_arr = np.stack(all_rec_d)
        spec_d_arr = np.stack(all_spec_d)
        support_d = cm_d_mean[:len(GENES), :].sum(axis=1) if all_cm_d else support_mean

        for i, gene in enumerate(GENES):
            rows_b.append({
                "Approach": f"Direct ({aligner})", "Gene": gene,
                "Precision_Mean": round(prec_d_arr[:, i].mean(), 4),
                "Precision_Std": round(prec_d_arr[:, i].std(), 4),
                "Recall_Mean": round(rec_d_arr[:, i].mean(), 4),
                "Recall_Std": round(rec_d_arr[:, i].std(), 4),
                "Specificity_Mean": round(spec_d_arr[:, i].mean(), 4),
                "Specificity_Std": round(spec_d_arr[:, i].std(), 4),
            })
        rows_b.append({
            "Approach": f"Direct ({aligner})", "Gene": "Overall",
            "Precision_Mean": round(np.average(prec_d_arr.mean(axis=0), weights=support_d), 4),
            "Precision_Std": round(np.average(prec_d_arr.std(axis=0), weights=support_d), 4),
            "Recall_Mean": round(np.average(rec_d_arr.mean(axis=0), weights=support_d), 4),
            "Recall_Std": round(np.average(rec_d_arr.std(axis=0), weights=support_d), 4),
            "Specificity_Mean": round(np.average(spec_d_arr.mean(axis=0), weights=support_d), 4),
            "Specificity_Std": round(np.average(spec_d_arr.std(axis=0), weights=support_d), 4),
        })

    pd.DataFrame(rows_b).to_csv(out_dir / "panel_B_metrics.csv", index=False)
    print(f"Panel B: {out_dir / 'panel_B_metrics.csv'}")


if __name__ == "__main__":
    main()
