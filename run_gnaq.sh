#!/bin/bash
# Run GNAQ analysis: download reads and run ParaDISM
# Usage: bash run_gnaq.sh

set -euo pipefail

echo "=========================================="
echo "GNAQ Analysis Pipeline"
echo "=========================================="

# Download reads if not present
if [ ! -d "gnaq_analysis/gnaq_reads" ] || [ -z "$(ls -A gnaq_analysis/gnaq_reads 2>/dev/null)" ]; then
    echo "Downloading GNAQ reads..."
    bash gnaq_analysis/download_gnaq_reads.sh
else
    echo "GNAQ reads already present"
fi

# Run ParaDISM on samples
echo ""
echo "Running ParaDISM on GNAQ samples..."
bash gnaq_analysis/run_gnaq_samples.sh

echo ""
echo "=========================================="
echo "GNAQ Analysis Complete!"
echo "=========================================="
