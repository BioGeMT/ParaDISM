#!/bin/bash

set -euo pipefail

# Script to run all 6 GNAQ scripts in parallel

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "Running all 6 GNAQ scripts in parallel"
echo "=========================================="
echo ""

# Array of scripts to run
SCRIPTS=(
    "run_gnaq-gnaqp_bwa.sh"
    "run_gnaq-gnaqp_bowtie2.sh"
    "run_gnaq-gnaqp_minimap2.sh"
    "run_gnaqp-gnapSNPs_bwa.sh"
    "run_gnaqp-gnapSNPs_bowtie2.sh"
    "run_gnaqp-gnapSNPs_minimap2.sh"
)

# Start time
start_time=$(date +%s)
start_date=$(date '+%Y-%m-%d %H:%M:%S')
echo "Start time: $start_date"
echo ""

# Array to store PIDs
declare -a PIDS=()

# Launch all scripts in parallel
for script in "${SCRIPTS[@]}"; do
    if [ ! -f "$script" ]; then
        echo "Warning: Script $script not found, skipping..."
        continue
    fi
    
    echo "Starting: $script"
    bash "$script" &
    PIDS+=($!)
done

echo ""
echo "All scripts started. Waiting for completion..."
echo "PIDs: ${PIDS[@]}"
echo ""

# Wait for all processes to complete
FAILED=0
for i in "${!PIDS[@]}"; do
    pid=${PIDS[$i]}
    script=${SCRIPTS[$i]}
    if wait $pid; then
        echo "✓ Completed: $script"
    else
        echo "✗ Failed: $script"
        FAILED=$((FAILED + 1))
    fi
done

# End time
end_time=$(date +%s)
end_date=$(date '+%Y-%m-%d %H:%M:%S')
duration=$((end_time - start_time))
hours=$((duration / 3600))
minutes=$(((duration % 3600) / 60))
seconds=$((duration % 60))

echo ""
echo "=========================================="
echo "All scripts finished"
echo "=========================================="
echo "Start time: $start_date"
echo "End time: $end_date"
printf "Total duration: %02d:%02d:%02d (%d seconds)\n" $hours $minutes $seconds $duration
echo ""

if [ $FAILED -eq 0 ]; then
    echo "✓ All scripts completed successfully!"
    exit 0
else
    echo "✗ $FAILED script(s) failed"
    exit 1
fi

