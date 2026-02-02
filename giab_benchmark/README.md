# GIAB HG002 Benchmark

Benchmarking ParaDISM against GIAB HG002 ground truth for the PKD1 gene region.

## Quick Start

```bash
# From giab_benchmark/
bash run_giab.sh [threads] [iterations]
```

## Scripts

### Download and Preparation

```bash
# Download GIAB HG002 Illumina reads
bash download_giab_hg002_reads.sh

# Merge lane chunks into HG002_R1.fq.gz / HG002_R2.fq.gz
bash merge_giab_hg002_reads.sh

# Prepare ground truth VCF (downloads benchmark, extracts PKD1 SNPs)
bash prepare_giab_truth.sh
```

### Run ParaDISM

```bash
# Run with G60 threshold
bash run_giab.sh [threads] [iterations]
```

Default: 8 threads, 10 iterations, min-alternate-count 5, bowtie2

### Variant Calling and Comparison

```bash
# Call raw variants (G60 only)
bash call_variants_raw_g60.sh

# Filter variants (SNPs only, GT=1/1)
bash filter_variants_g60.sh

# Compare to ground truth
python compare_variants.py \
    --truth giab_hg002_vcf/HG002_PKD1_genes_SNPs_exact.vcf.gz \
    --called <output_dir>/final_outputs/variant_calling_balanced/variants.vcf.gz

# Plot precision/recall comparison (G60 only)
python plot_balanced_comparison.py
```

## Filters

Final filtering used: SNPs only, GT=1/1

## Output Structure

```
giab_benchmark/
├── giab_hg002_reads/                              # Downloaded reads
├── giab_hg002_vcf/                                # Ground truth VCF
└── giab_hg002_output_bowtie2_G60_min5_qfilters/   # G60 results
```
