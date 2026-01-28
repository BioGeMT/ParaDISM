# GIAB HG002 Benchmark

Benchmarking ParaDISM against GIAB HG002 ground truth for the PKD1 gene region.

## Quick Start

```bash
# From repo root
bash run_giab.sh [threads] [iterations]
```

## Scripts

### Download and Preparation

```bash
# Download GIAB HG002 Illumina reads
bash giab_benchmark/download_giab_hg002_reads.sh

# Prepare ground truth VCF (downloads benchmark, extracts PKD1 SNPs)
bash giab_benchmark/prepare_giab_truth.sh
```

### Run ParaDISM

```bash
# Run with G30 and G60 thresholds
bash giab_benchmark/run_giab.sh [threads] [iterations]
```

Default: 8 threads, 10 iterations, min-alternate-count 5, bowtie2

### Variant Calling and Comparison

```bash
# Call variants with balanced filters
bash giab_benchmark/call_variants_balanced_filters.sh

# Compare to ground truth
python compare_variants.py \
    --truth giab_benchmark/giab_hg002_vcf/HG002_PKD1_genes_SNPs_illumina_only.vcf.gz \
    --called <output_dir>/final_outputs/variant_calling_balanced/variants.vcf.gz

# Plot precision/recall comparison
python giab_benchmark/plot_balanced_comparison.py
```

## Filters

Balanced variant filters: QUAL>=500, DP>=30, AF>=0.90, RO<=5, MQMR<=10

## Output Structure

```
giab_benchmark/
├── giab_hg002_reads/                              # Downloaded reads
├── giab_hg002_vcf/                                # Ground truth VCF
├── giab_hg002_output_bowtie2_G30_min5_qfilters/   # G30 results
└── giab_hg002_output_bowtie2_G60_min5_qfilters/   # G60 results
```
