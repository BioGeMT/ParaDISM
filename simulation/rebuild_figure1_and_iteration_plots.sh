#!/usr/bin/env bash
set -euo pipefail

# Rebuild all figure1 inputs/outputs and overwrite existing files:
#  - Per-aligner figure1 CSVs (panel_A_confusion_matrix.csv, panel_B_metrics.csv)
#  - Combined all-aligner CSVs (panel_A_all_aligners.csv, panel_B_all_aligners.csv)
#  - Iteration progression CSV/PNG for each aligner
#  - Final 4-panel paper figure SVG/PNG
#
# Usage examples:
#   bash simulation/rebuild_figure1_and_iteration_plots.sh
#   bash simulation/rebuild_figure1_and_iteration_plots.sh --seed-start 1 --seed-end 1000
#   bash simulation/rebuild_figure1_and_iteration_plots.sh --sim-output-base simulation/simulations_outputs_100
#
# Optional env overrides:
#   CONDA_ENV=paradism_env

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

SIM_OUTPUT_BASE="simulation/simulations_outputs_1000"
SEED_START=1
SEED_END=1000
ALIGNERS="bwa-mem2 bowtie2 minimap2"
CONDA_ENV="${CONDA_ENV:-paradism_env}"
FIGURE1_DATA_DIR=""
ITERATION_PLOTS_DIR=""
OUTPUT_SVG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sim-output-base)
            SIM_OUTPUT_BASE="$2"
            shift 2
            ;;
        --seed-start)
            SEED_START="$2"
            shift 2
            ;;
        --seed-end)
            SEED_END="$2"
            shift 2
            ;;
        --aligners)
            ALIGNERS="$2"
            shift 2
            ;;
        --figure1-data-dir)
            FIGURE1_DATA_DIR="$2"
            shift 2
            ;;
        --iteration-plots-dir)
            ITERATION_PLOTS_DIR="$2"
            shift 2
            ;;
        --output-svg)
            OUTPUT_SVG="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        -h|--help)
            sed -n '1,42p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

FIGURE1_DATA_DIR="${FIGURE1_DATA_DIR:-${SIM_OUTPUT_BASE}/figure1_data}"
ITERATION_PLOTS_DIR="${ITERATION_PLOTS_DIR:-${SIM_OUTPUT_BASE}/iteration_plots}"
OUTPUT_SVG="${OUTPUT_SVG:-${FIGURE1_DATA_DIR}/paper_multipanel_figure.svg}"

if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda not found in PATH." >&2
    exit 1
fi

run_py() {
    conda run -n "${CONDA_ENV}" python "$@"
}

echo "Using:"
echo "  SIM_OUTPUT_BASE=${SIM_OUTPUT_BASE}"
echo "  SEED_START=${SEED_START}"
echo "  SEED_END=${SEED_END}"
echo "  ALIGNERS=${ALIGNERS}"
echo "  FIGURE1_DATA_DIR=${FIGURE1_DATA_DIR}"
echo "  ITERATION_PLOTS_DIR=${ITERATION_PLOTS_DIR}"
echo "  OUTPUT_SVG=${OUTPUT_SVG}"
echo "  CONDA_ENV=${CONDA_ENV}"

mkdir -p "${FIGURE1_DATA_DIR}" "${ITERATION_PLOTS_DIR}"
read -r -a ALIGNER_LIST <<< "${ALIGNERS}"

echo
echo "[1/4] Rebuilding per-aligner Figure1 CSVs"
for aligner in "${ALIGNER_LIST[@]}"; do
    out_dir="${FIGURE1_DATA_DIR}/${aligner}"
    mkdir -p "${out_dir}"
    run_py simulation/build_figure1_csv.py \
        --sim-output-base "${SIM_OUTPUT_BASE}" \
        --aligner "${aligner}" \
        --seed-start "${SEED_START}" \
        --seed-end "${SEED_END}" \
        --output-dir "${out_dir}"
done

echo
echo "[2/4] Rebuilding combined all-aligner CSVs"
python - "${FIGURE1_DATA_DIR}" "${ALIGNERS}" <<'PY'
import re
import sys
from pathlib import Path
import pandas as pd

figure1_data_dir = Path(sys.argv[1])
aligners = sys.argv[2].split()

panel_a_frames = []
panel_b_frames = []

for aligner in aligners:
    panel_a = figure1_data_dir / aligner / "panel_A_confusion_matrix.csv"
    panel_b = figure1_data_dir / aligner / "panel_B_metrics.csv"
    if not panel_a.exists():
        raise FileNotFoundError(f"Missing per-aligner panel A CSV: {panel_a}")
    if not panel_b.exists():
        raise FileNotFoundError(f"Missing per-aligner panel B CSV: {panel_b}")

    a = pd.read_csv(panel_a)
    a["Aligner"] = aligner
    a["Approach"] = a["Approach"].map(
        lambda x: "ParaDISM" if str(x).strip() == "ParaDISM" else "Direct"
    )
    panel_a_frames.append(
        a[["Aligner", "Approach", "Origin_Gene", "Mapped_Gene", "Mean_Reads", "Std_Reads"]]
    )

    b = pd.read_csv(panel_b)
    b["Aligner"] = aligner
    b["Approach"] = b["Approach"].map(
        lambda x: "ParaDISM" if str(x).strip() == "ParaDISM" else "Direct"
    )
    panel_b_frames.append(
        b[
            [
                "Aligner",
                "Approach",
                "Gene",
                "Precision_Mean",
                "Precision_Std",
                "Recall_Mean",
                "Recall_Std",
                "Specificity_Mean",
                "Specificity_Std",
            ]
        ]
    )

panel_a_all = pd.concat(panel_a_frames, ignore_index=True)
panel_b_all = pd.concat(panel_b_frames, ignore_index=True)

panel_a_out = figure1_data_dir / "panel_A_all_aligners.csv"
panel_b_out = figure1_data_dir / "panel_B_all_aligners.csv"
panel_a_all.to_csv(panel_a_out, index=False)
panel_b_all.to_csv(panel_b_out, index=False)

print(f"Wrote {panel_a_out}")
print(f"Wrote {panel_b_out}")
PY

echo
echo "[3/4] Rebuilding iteration progression plots/CSVs"
for aligner in "${ALIGNER_LIST[@]}"; do
    run_py simulation/plot_iteration_progression.py \
        --sim-output-base "${SIM_OUTPUT_BASE}" \
        --seed-start "${SEED_START}" \
        --seed-end "${SEED_END}" \
        --output-dir "${ITERATION_PLOTS_DIR}" \
        --aligner "${aligner}"
done

echo
echo "[4/4] Rebuilding 4-panel figure"
run_py simulation/make_paper_multipanel_plot.py \
    --figure1-data-dir "${FIGURE1_DATA_DIR}" \
    --iteration-plots-dir "${ITERATION_PLOTS_DIR}" \
    --output "${OUTPUT_SVG}"

echo
echo "Done. Overwritten outputs:"
echo "  ${FIGURE1_DATA_DIR}/panel_A_all_aligners.csv"
echo "  ${FIGURE1_DATA_DIR}/panel_B_all_aligners.csv"
echo "  ${ITERATION_PLOTS_DIR}/iteration_progression_<aligner>.csv"
echo "  ${ITERATION_PLOTS_DIR}/iteration_progression_<aligner>.png"
echo "  ${OUTPUT_SVG}"
echo "  ${OUTPUT_SVG%.svg}.png"
