# GNAQ Analysis

Analysis of GNAQ/GNAQ pseudogene samples using ParaDISM.

## Quick Start

```bash
# From gnaq_analysis/
bash download_gnaq_reads.sh
bash run_gnaq_samples_bwa_180.sh
```

## Scripts

```bash
# Download GNAQ reads from EBI
bash download_gnaq_reads.sh

# Run ParaDISM on all samples with bwa-mem2 threshold=180
bash run_gnaq_samples_bwa_180.sh

# Base counts at selected positions (outputs to base_counts_180/)
bash make_base_counts_csvs.sh bwa-mem2_180_output
```

## Samples

- SRR5602384
- SRR5602389
- SRR5602393
- SRR5602414
- SRR5602419

## Configuration

- Aligner: bwa-mem2 (threshold=180)
- Threads: 12
- Iterations: 10
- Reference: gnaq-gnaqp_ref.fa

## Output Structure

```
gnaq_analysis/
├── GNAQ_reads/                  # Downloaded FASTQ files
├── bwa-mem2_180_output/         # ParaDISM results
└── base_counts_180/             # Base counts at selected positions
```
