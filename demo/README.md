# ParaDISM Demo

This directory contains a self-contained synthetic example that does not use
manuscript or private data:

- `ref.fa`: two short homologous reference contigs
- `reads_R1.fq` and `reads_R2.fq`: four 150 bp paired-end reads with gene-informative bases
- `run_demo.sh`: runs ParaDISM with Bowtie2 for one iteration and checks the expected outputs
- `generate_demo_reads.sh`: optional helper showing how demo reads can be regenerated with `dwgsim`

From the repository root:

```bash
mamba env create -f environment.yml
mamba activate paradism
bash demo/run_demo.sh
```

`conda` can be used instead of `mamba`:

```bash
conda env create -f environment.yml
conda activate paradism
```

By default, the demo writes to `demo/output/`. Set
`OUTPUT_DIR=/path/to/output` before running the script to choose another
location. The default output directory is removed and recreated on each run.

Expected output layout:

```text
demo/output/
├── mapped_reads.sam
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

The demo does not require `dwgsim`; it uses the committed FASTQ files. To
regenerate a small synthetic read set separately, run:

```bash
bash demo/generate_demo_reads.sh
```
