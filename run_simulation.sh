#!/bin/bash
# Run simulation analysis: simulate reads and run ParaDISM
# Usage: bash run_simulation.sh [seed_start] [seed_end]

set -euo pipefail

SEED_START="${1:-1}"
SEED_END="${2:-10}"

echo "=========================================="
echo "Simulation Analysis Pipeline"
echo "=========================================="

# Run simulation and mapping
echo ""
echo "Running simulations (seeds $SEED_START to $SEED_END)..."
SEED_START="$SEED_START" SEED_END="$SEED_END" bash simulation/simulate_and_map.sh

echo ""
echo "=========================================="
echo "Simulation Analysis Complete!"
echo "=========================================="
