#!/bin/bash

# Quick test wrapper for simulate_and_map.sh with 10 seeds and 10 iterations.
# Uses the main script's env overrides so the defaults in simulate_and_map.sh remain unchanged.
# Adjust values below if you want a different seed range or output location.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Override seed range and output location for this test run
export SEED_START=1
export SEED_END=10
export ITERATIONS=10
export SIM_OUTPUT_BASE="${SIM_OUTPUT_BASE:-${SCRIPT_DIR}/sim_output_test_10_seeds}"

# Optional: narrow aligners or threads (uncomment to customize)
# export ALIGNERS="bwa-mem2"
# export THREADS=2

"${SCRIPT_DIR}/simulate_and_map.sh"
