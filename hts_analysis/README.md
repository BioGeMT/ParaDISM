# HTS Clinical Samples Analysis

Analysis of HTS clinical samples for PKD1 using ParaDISM.

## Quick Start

```bash
# From repo root
bash run_hts.sh
```

## Scripts

```bash
# Run ParaDISM on all HTS samples
bash run_hts_samples.sh
```

## Samples

20 samples from shared storage (`/mnt/STORAGE-BioGeMT-01/pkd_data`):
- HTS009
- HTS131_L1, HTS131_L2
- HTS132_L1, HTS132_L2
- HTS158-171
- IndexCHKPEPI00000968_HTS002

## Configuration

- Aligner: bowtie2
- Threads: 8
- Iterations: 10
- Reference: ref.fa

## Output Structure

```
hts_analysis/
├── HTS_bowtie2_output/  # New runs
├── HTS_CNTRL/           # Control outputs (multiple aligners)
└── HTS_output/          # Previous run outputs
```
