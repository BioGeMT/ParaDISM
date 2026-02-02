# GNAQ Analysis

Analysis of GNAQ/GNAQ pseudogene samples using ParaDISM.

## Quick Start

```bash
# From gnaq_analysis/
bash download_gnaq_reads.sh
bash run_gnaq_samples_bwa_160.sh
```

## Scripts

```bash
# Download GNAQ reads from EBI
bash download_gnaq_reads.sh

# Run ParaDISM on all samples with bwa-mem2 threshold=160 (100bp equivalent)
bash run_gnaq_samples_bwa_160.sh

# (Optional) Base counts at selected positions for any BAM dir
bash make_base_counts_csvs.sh bwa-mem2_160_output/SRR5602384/final_outputs
```

## Samples

- SRR5602384
- SRR5602389
- SRR5602393
- SRR5602414
- SRR5602419

## Configuration

- Aligner: bwa-mem2 (threshold=160)
- Threads: 12
- Iterations: 10
- Reference: gnaq-gnaqp_ref.fa

## Output Structure

```
gnaq_analysis/
├── GNAQ_reads/                  # Downloaded FASTQ files
└── bwa-mem2_160_output/          # ParaDISM results per sample
```
