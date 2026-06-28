# ParaDISM Demo

This directory contains a self-contained synthetic example that does not use
manuscript or private data:

- `ref.fa`: two short homologous reference contigs
- `generate_demo_reads.sh`: generates 20 paired-end read pairs with `dwgsim`
- `run_demo.sh`: generates missing reads, runs ParaDISM with Bowtie2 for one
  iteration, and checks the expected outputs

From the repository root:

```bash
conda env create -f environment.yml
conda activate paradism
bash demo/run_demo.sh
```

By default, `run_demo.sh` writes generated reads to `demo/generated_reads/` and
ParaDISM outputs to `demo/output/`. Set `OUTPUT_DIR=/path/to/output` before
running the script to choose another ParaDISM output location. The default
output directory is removed and recreated on each run.

Expected output layout:

```text
demo/output/
├── iteration_1/
│   └── mapped_reads.sam
└── final_outputs/
    ├── tiny_demo_fastq/
    │   ├── tiny_demo_PARA1.fq
    │   └── tiny_demo_PARA2.fq
    └── tiny_demo_bam/
        ├── tiny_demo_PARA1.sorted.bam
        ├── tiny_demo_PARA1.sorted.bam.bai
        ├── tiny_demo_PARA2.sorted.bam
        └── tiny_demo_PARA2.sorted.bam.bai
```

To regenerate the small synthetic read set separately, run:

```bash
bash demo/generate_demo_reads.sh
```
