# ParaDISM Demo

This directory contains a self-contained synthetic example that does not use
manuscript or private data:

- `ref.fa`: PKD1 plus six PKD1 pseudogene reference contigs
- `pkd1_sim_R1.fq` and `pkd1_sim_R2.fq`: 20 committed paired-end read pairs
- `run_demo.sh`: runs ParaDISM with Bowtie2 for one iteration and checks the
  expected outputs
- `generate_demo_reads.sh`: optional script to regenerate the committed reads
  from `ref.fa` with `dwgsim`

From the repository root:

```bash
conda env create -f environment.yml
conda activate paradism
bash demo/run_demo.sh
```

By default, `run_demo.sh` uses the committed reads and writes ParaDISM outputs
to `demo/output/`. Set `OUTPUT_DIR=/path/to/output` before running the script
to choose another ParaDISM output location. The default output directory is
removed and recreated on each run.

Expected output layout:

```text
demo/output/
├── iteration_1/
│   └── mapped_reads.sam
└── final_outputs/
    ├── pkd1_demo_fastq/
    │   └── pkd1_demo_<assigned_contig>.fq
    └── pkd1_demo_bam/
        ├── pkd1_demo_<assigned_contig>.sorted.bam
        └── pkd1_demo_<assigned_contig>.sorted.bam.bai
```

To regenerate the small synthetic read set separately, run:

```bash
bash demo/generate_demo_reads.sh
```

To regenerate the reads before running the demo, run:

```bash
REGENERATE_READS=1 bash demo/run_demo.sh
```
