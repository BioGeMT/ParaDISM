#!/bin/bash
# Run HTS analysis: run ParaDISM on HTS samples
# Usage: bash run_hts.sh

set -euo pipefail

echo "=========================================="
echo "HTS Analysis Pipeline"
echo "=========================================="

# Run ParaDISM on samples
echo ""
echo "Running ParaDISM on HTS samples..."
bash hts_analysis/run_hts_samples.sh

echo ""
echo "=========================================="
echo "HTS Analysis Complete!"
echo "=========================================="
