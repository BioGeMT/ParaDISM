# GIAB Benchmark

This directory contains scripts for public HG002/GIAB benchmarking. GIAB runs
are reproducible but heavyweight because the public HG002 read and truth inputs
are large. Use the top-level `demo/` and `benchmark/simulation/` workflows for
fast checks.

## Inputs

The GIAB workflow needs:

- public HG002 Illumina FASTQ shards, downloaded by
  `download_giab_hg002_reads.sh`
- GIAB HG002 truth VCF and index
- GIAB benchmark BED
- one group reference from `benchmark/references/`

The existing scripts in this directory are the PKD1/GIAB workflow used for the
manuscript benchmark. They are kept here rather than at the repository root so
they are clearly separated from the minimal demo.

## Existing Script Roles

- `download_giab_hg002_reads.sh`: downloads public HG002 read inputs used by
  the benchmark. The read data are large; inspect the script before running on
  a local machine.
- `prepare_giab_truth.sh`: prepares GIAB truth and benchmark regions in the
  gene-coordinate system expected by ParaDISM evaluation.
- `run_giab.sh`: runs ParaDISM on the prepared GIAB read set.
- `run_variant_calling_to_vcf_out.sh`: calls and evaluates variants from
  ParaDISM and direct/base-aligner outputs.
- `call_variants_raw_g60.sh`, `filter_simple_snps_acgt_final.sh`, and helper
  Python scripts: variant calling/filtering utilities used by the benchmark.

## Minimal Command Sequence

From the repository root, after installing `environment.yml`:

```bash
bash benchmark/giab/prepare_giab_truth.sh
bash benchmark/giab/download_giab_hg002_reads.sh
bash benchmark/giab/run_giab.sh --threads 8 --iterations 10

RUN_DIR=benchmark/giab/giab_hg002_output_bowtie2_G60_min5_qfilters
OUT_DIR=benchmark/giab/vcf_out_full
bash benchmark/giab/run_variant_calling_to_vcf_out.sh \
  --run-dir "$RUN_DIR" \
  --out-dir "$OUT_DIR" \
  --threads 8
```

This workflow uses public data, but it is not intended to be a quick smoke
test. The demo and simulation benchmark are the intended lightweight
reproducibility checks.
