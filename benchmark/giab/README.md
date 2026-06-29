# GIAB Benchmark

This directory contains scripts for public HG002/GIAB benchmarking.

The workflow analyzes PKD1/PKD1P SNP calls from ParaDISM-refined BAMs versus
the base Bowtie2 BAM, using GIAB HG002 truth as the benchmark.

## Inputs

- public HG002 Illumina FASTQ shards, downloaded by
  `download_giab_hg002_reads.sh`
- GIAB HG002 truth VCF and index
- GIAB benchmark BED
- `benchmark/references/pkd1_panel.fa`

## Usage

From the repository root, after installing `environment.yml`:

```bash
bash benchmark/giab/prepare_giab_truth.sh
bash benchmark/giab/download_giab_hg002_reads.sh
bash benchmark/giab/run_giab.sh --threads 8 --workers 2 --iterations 1

RUN_DIR=benchmark/giab/giab_hg002_output_bowtie2_G60_min5_qfilters
OUT_DIR=benchmark/giab/vcf_out_full
bash benchmark/giab/run_variant_calling_to_vcf_out.sh \
  --run-dir "$RUN_DIR" \
  --out-dir "$OUT_DIR" \
  --threads 8
```

For a smaller 10x-style read subset, run `downsample_hg002_fastq.sh` after
downloading the HG002 shards, then pass the downsampled read directory to
`run_giab.sh` with `--reads-dir`.

## Output Layout

```text
benchmark/giab/
├── giab_hg002_reads/                 # downloaded or downsampled FASTQs
├── giab_hg002_vcf/                   # prepared GIAB truth and BED files
├── giab_hg002_output_*/              # ParaDISM run output
│   ├── iteration_1/
│   └── final_outputs/
└── vcf_out*/                         # evaluation CSV/JSON outputs
```

Key output files:

- `variant_calling_metrics.csv` and `.json`
- `confusion_matrices_overall.csv`
- `per_gene_metrics/per_gene_metrics.csv` and `.json`
- `confusion_matrices_by_coverage.csv`
- `coverage_split_confusion_metrics.json`

## Scripts

- `download_giab_hg002_reads.sh`: downloads public HG002 read inputs used by
  the benchmark.
- `downsample_hg002_fastq.sh`: creates a smaller read subset from downloaded
  HG002 shards.
- `prepare_giab_truth.sh`: prepares GIAB truth and benchmark regions in
  ParaDISM gene coordinates.
- `run_giab.sh`: runs ParaDISM on the GIAB read set.
- `run_variant_calling_to_vcf_out.sh`: calls variants from the completed
  ParaDISM run, filters to benchmarkable simple SNPs, and writes metrics.

The workflow uses public data and writes ParaDISM outputs, filtered VCFs, and
evaluation metrics under the selected output directories.
