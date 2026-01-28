#!/bin/bash
# Run GIAB benchmark: download reads, prepare truth, run ParaDISM
# Usage: bash run_giab.sh [threads] [iterations]

set -euo pipefail

THREADS="${1:-8}"
ITERATIONS="${2:-10}"

echo "=========================================="
echo "GIAB HG002 Benchmark Pipeline"
echo "=========================================="

# Download reads if not present
if [ ! -d "giab_benchmark/giab_hg002_reads" ] || [ -z "$(ls -A giab_benchmark/giab_hg002_reads 2>/dev/null)" ]; then
    echo "Downloading GIAB HG002 reads..."
    bash giab_benchmark/download_giab_hg002_reads.sh
else
    echo "GIAB reads already present"
fi

# Prepare ground truth VCF
echo ""
echo "Preparing ground truth VCF..."
bash giab_benchmark/prepare_giab_truth.sh

# Run ParaDISM with G30 and G60 thresholds
echo ""
echo "Running ParaDISM (threads=$THREADS, iterations=$ITERATIONS)..."
bash giab_benchmark/run_giab.sh "$THREADS" "$ITERATIONS"

echo ""
echo "=========================================="
echo "GIAB Benchmark Complete!"
echo "=========================================="
echo ""
echo "Next: Call variants with balanced filters"
echo "  bash giab_benchmark/call_variants_balanced_filters.sh"
