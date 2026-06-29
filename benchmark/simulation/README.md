# Simulation Benchmark

This workflow simulates paired-end reads from one reference FASTA with
`dwgsim`, runs ParaDISM and direct aligners on the same reads, and summarizes
read-assignment precision, sensitivity, and specificity.

The reference FASTA should contain one gene/paralog group: the gene of interest
and homologous sequences that reads may be assigned to. Because the reads are
simulated, the true source contig of each read is known.

## Usage

Run one group:

```bash
bash benchmark/simulation/run_simulation.sh \
  --group hba_pair \
  --reference benchmark/references/hba_pair.fa \
  --seeds 2 \
  --reads-per-seed 1000 \
  --aligners bowtie2,bwa-mem2,minimap2 \
  --out results/simulation/hba_pair
```

Or run the default small benchmark:

```bash
bash benchmark/simulation/run_simulation.sh
```

The default uses `pkd1_panel`, one seed, and 1,000 read pairs.

## Options

- `--group`: label for the gene/paralog group
- `--reference`: FASTA to simulate from; defaults to
  `benchmark/references/<group>.fa`
- `--out`: output directory
- `--seeds`: run seeds `1..N`
- `--reads-per-seed`: read pairs generated per seed
- `--aligners`: comma or space separated list
- `--threads`: threads passed to the underlying aligner
- `--workers`: worker processes for ParaDISM read assignment
- `--iterations`: ParaDISM iterations (default: 1)

## Output Layout

```text
results/simulation/<group>/
├── seed_*/                              # simulated reads and per-aligner runs
├── aggregated_results/
│   ├── read_mapping_aggregated_summary.csv
│   └── per_seed_overall_metrics.csv
└── timing_data.csv
```

The aggregate CSVs report per-gene and overall precision, sensitivity, and
specificity across seeds.
