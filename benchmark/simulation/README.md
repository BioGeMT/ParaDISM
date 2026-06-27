# Simulation Benchmark

This workflow simulates paired-end reads from one reference FASTA with
`dwgsim`, runs ParaDISM and direct aligners on the same reads, and summarizes
read-assignment precision/recall/specificity.

Run one group at a time:

```bash
bash benchmark/simulation/run_simulation.sh \
  --group hba_pair \
  --reference benchmark/references/hba_pair.fa \
  --seeds 2 \
  --reads-per-seed 1000 \
  --aligners bowtie2,bwa-mem2,minimap2 \
  --out results/simulation/hba_pair
```

The default command is intentionally small:

```bash
bash benchmark/simulation/run_simulation.sh
```

It runs `pkd1_panel` for one seed with 1,000 read pairs. Full paper-scale runs
used larger seed/read counts and can be reproduced by increasing
`--seeds` and `--reads-per-seed`.

For a faster smoke test:

```bash
bash benchmark/simulation/run_simulation.sh \
  --group hba_pair \
  --reference benchmark/references/hba_pair.fa \
  --seeds 1 \
  --reads-per-seed 20 \
  --aligners bowtie2 \
  --iterations 1 \
  --threads 1 \
  --out results/simulation/smoke_hba_pair
```

## Main Script

`run_simulation.sh`

- generates reads with `dwgsim`
- runs ParaDISM for each requested aligner
- reuses the initial ParaDISM SAM as the direct-aligner comparison
- runs per-seed read-assignment analysis
- writes aggregate metric CSVs and timing data

Important options:

- `--group`: label for the gene/paralog group
- `--reference`: FASTA to simulate from
- `--out`: output directory
- `--seeds`: run seeds `1..N`
- `--reads-per-seed`: read pairs generated per seed
- `--aligners`: comma or space separated list
- `--iterations`: ParaDISM iterations (default: 1)

Outputs are written under the selected `--out` directory:

```text
results/simulation/<group>/
├── seed_*/
├── aggregated_results/
└── timing_data.csv
```
