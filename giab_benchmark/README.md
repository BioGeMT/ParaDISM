# GIAB HG002 Benchmark

Benchmarking ParaDISM against GIAB HG002 ground truth for the PKD1 gene region.

## Quick Start (Full HG002)

```bash
# From repo root
bash giab_benchmark/prepare_giab_truth.sh
bash giab_benchmark/run_giab.sh [threads] [iterations]

RUN_DIR=giab_benchmark/giab_hg002_output_bowtie2_G60_min5_qfilters
OUT_DIR=giab_benchmark/vcf_out_full
bash giab_benchmark/run_variant_calling_to_vcf_out.sh \
  --run-dir "$RUN_DIR" \
  --out-dir "$OUT_DIR" \
  --threads 8
```

## Scripts

### Download and Preparation

```bash
# Download GIAB HG002 Illumina reads
bash giab_benchmark/download_giab_hg002_reads.sh

# Prepare ground truth VCF/BED in gene coordinates
bash giab_benchmark/prepare_giab_truth.sh
```

### Run ParaDISM

```bash
# Run with G60 threshold
bash giab_benchmark/run_giab.sh [threads] [iterations]
```

Default: 8 threads, 10 iterations, min-alternate-count 5, bowtie2

### Post-processing (variant calling + filtering + plots)

```bash
RUN_DIR=giab_benchmark/giab_hg002_output_bowtie2_G60_min5_qfilters
OUT_DIR=giab_benchmark/vcf_out_full
bash giab_benchmark/run_variant_calling_to_vcf_out.sh \
  --run-dir "$RUN_DIR" \
  --out-dir "$OUT_DIR" \
  --threads 8
```

### Downsampled runs (optional)

```bash
# Example: create 10x subset (2x250 => 500 bp/pair)
bash giab_benchmark/downsample_hg002_fastq.sh \
  --coverage 10 \
  --bp-per-pair 500 \
  --out-dir giab_benchmark/giab_hg002_reads/subset_genome10x_bp500_balanced_v1

# Special experiment: 10x without iterative quality filters
bash giab_benchmark/run_giab_10x_noqual.sh --threads 8 --iterations 10
```

## Filtering Summary

- `call_variants_raw_g60.sh` generates raw calls and simple biallelic SNP A/C/G/T calls.
- Calls are normalized/atomized (`bcftools norm -m -any` + `atomize_equal_length_substitutions.py`).
- Final benchmark VCFs are produced by `filter_simple_snps_acgt_final.sh` using:
  - `INFO/AO>=3 && QUAL>=10 && INFO/SAF>=1 && INFO/SAR>=1`
  - benchmark BED mask (`HG002_PKD1_genes_benchmarkable_gene_coords.bed`)

## Output Structure

```
giab_benchmark/
├── giab_hg002_reads/                              # Downloaded reads
├── giab_hg002_vcf/                                # Ground truth VCF
├── giab_hg002_output_bowtie2_G60_min5_qfilters/   # G60 results
    └── variant_calling/                           # Raw/simple/final VCFs
└── vcf_out_full/                                  # Metrics + confusion matrices + plots
```
