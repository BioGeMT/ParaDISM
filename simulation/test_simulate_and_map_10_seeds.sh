#!/bin/bash

# Quick test wrapper for simulate_and_map.sh with a single seed and focused settings.
# Uses env overrides so the defaults in simulate_and_map.sh remain unchanged.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Override seed range, iterations, reads, and output location for this test run
export SEED_START=1
export SEED_END=1
export ITERATIONS=10
export NUM_READS=5000
export SIM_OUTPUT_BASE="${SIM_OUTPUT_BASE:-${SCRIPT_DIR}/sim_output_test}"

# Optional: narrow aligners or threads (uncomment to customize)
export ALIGNERS="bwa-mem2 bowtie2 minimap2"
# export THREADS=2

"${SCRIPT_DIR}/simulate_and_map.sh"
