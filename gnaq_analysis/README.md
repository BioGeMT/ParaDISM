# GNAQ Analysis

Analysis of GNAQ/GNAQ pseudogene samples using ParaDISM.

## Quick Start

```bash
# From repo root
bash run_gnaq.sh
```

## Scripts

```bash
# Download GNAQ reads from EBI
bash gnaq_analysis/download_gnaq_reads.sh

# Run ParaDISM on all samples
bash gnaq_analysis/run_gnaq_samples.sh
```

## Samples

- SRR5602384
- SRR5602389
- SRR5602393
- SRR5602414
- SRR5602419

## Configuration

- Aligner: minimap2 (short profile)
- Threads: 2
- Iterations: 10
- Reference: gnaq_analysis/gnaq-gnaqp_ref.fa

## Output Structure

```
gnaq_analysis/
├── gnaq_reads/           # Downloaded FASTQ files
└── GNAQ_minimap2_output/ # ParaDISM results per sample
```
