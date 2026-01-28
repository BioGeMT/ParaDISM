# Simulation Analysis

Simulated read analysis to evaluate ParaDISM accuracy across multiple aligners.

## Quick Start

```bash
# From repo root (default: seeds 1-10)
bash run_simulation.sh [seed_start] [seed_end]

# Example: run seeds 1-100
bash run_simulation.sh 1 100
```

## Scripts

```bash
# Simulate reads only (using DWGSIM)
bash simulation/run_dwgsim_simulation.sh

# Full pipeline: simulate + map + analyze
bash simulation/simulate_and_map.sh
```

## Configuration

Environment variables:
- `SEED_START`, `SEED_END`: Seed range (default: 1-1000)
- `THREADS`: Threads per run (default: 1)
- `NUM_READS`: Read pairs per seed (default: 100000)
- `ALIGNERS`: Space-separated list (default: "bwa-mem2 bowtie2 minimap2")
- `ITERATIONS`: ParaDISM iterations (default: 10)
- `ERROR_RATE`: Sequencing error rate (default: 0.01)
- `SNP_RATE`: SNP rate (default: 0.005)
- `INDEL_RATE`: Indel rate (default: 0.0005)

## Analysis Scripts

```bash
# Aggregate results across seeds
python aggregate_results.py --sim-output-base <dir> --seed-start 1 --seed-end 100

# Per-seed metrics
python aggregate_per_seed_overall.py --sim-output-base <dir>

# Plot iteration progression
python plot_iteration_progression.py --sim-output-base <dir> --aligner bowtie2

# Plot timing data
python plot_timings.py
```

## Output Structure

```
simulation/
├── dwgsim/               # DWGSIM tool
├── pkd_reads_dwgsim/     # Single-seed simulated reads
└── simulations_outputs/  # Multi-seed results
    ├── seed_*/           # Per-seed outputs
    ├── aggregated_results/
    ├── iteration_plots/
    └── timing_data.csv
```
