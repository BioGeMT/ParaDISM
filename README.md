# Homology Mapper Pipeline

Read-mapping and refinement workflow for highly homologous genomic regions. The pipeline aligns paired-end reads with MAFFT/Bowtie2/BWA-MEM2/Minimap2, maps alignments back to multiple-sequence alignments (MSAs), filters for uniquely supported reference sequences, and emits gene-specific FASTQ/BAM outputs.

## Prerequisites

### Environment Setup
To create a conda environment with all required dependencies, run:

```bash
conda env create -f mapper_env.yml
conda activate mapper_env
```

Key packages installed by the environment:
- Python 3.11
- MAFFT, Bowtie2, BWA-MEM2, Minimap2, Samtools
- Biopython, PySAM, Rich

## Directory Structure

```
.
├── mapper.py
├── mapper_env.yml
├── src/
│   ├── pipeline/
│   │   ├── executor.py
│   │   ├── mapper_algo.py
│   │   ├── mapper_algo_snp_only.py
│   │   ├── output.py
│   │   ├── read_2_gene.py
│   │   └── ref_2_msa.py
│   ├── ui/
│   │   ├── interactive.py
│   │   └── ui_components.py
│   └── utils/
│       ├── file_scanner.py
│       ├── logger.py
│       ├── progress.py
│       └── validators.py
└── output/                        # Created after pipeline run
```

## Usage

### Interactive Mode
Run without arguments to launch the guided CLI:

```bash
python mapper.py
```

You’ll be prompted to:
- Select FASTQ pairs (R1/R2)
- Choose reference FASTA
- Optionally reuse an existing SAM (skipping alignment)
- Pick an aligner (`bowtie2`, `bwa-mem2`, `minimap2`)
- Configure threads and minimap2 profile (if applicable)

### Argument-Driven Mode
```bash
python mapper.py --read1 <forward_reads.fq> \
                 --read2 <reverse_reads.fq> \
                 --reference <reference.fa> \
                 [--aligner ALIGNER] \
                 [--threads N] \
                 [--minimap2-profile PROFILE] \
                 [--sam existing.sam] \
                 [--output-dir OUTPUT]
```

Required:
- `--read1`, `--read2`: Paired FASTQ files
- `--reference`: Reference FASTA

Optional:
- `--aligner`: `bowtie2` (default), `bwa-mem2`, or `minimap2`
- `--threads`: Number of alignment threads (default: 4)
- `--minimap2-profile`: `short`, `pacbio`, or `nanopore` (required with minimap2)
- `--sam`: Existing alignment (skips alignment stage)
- `--output-dir`: Destination directory (default: `./output`)

### SAM Requirements
When skipping alignment via `--sam`, ensure the SAM:
- Contains a valid header
- Has at least one mapped read
- Includes MD tags (`samtools calmd` can add them)

## Pipeline Outline

1. **MSA generation** – MAFFT builds reference MSAs.
2. **Indexer** – Builds Bowtie2/BWA-MEM2/Minimap2 reference indexes.
3. **Read alignment** – Selected aligner maps paired reads.
4. **Reference→MSA mapping** – Captures reference base positions across the MSA.
5. **Read→reference mapping** – Parses SAM alignments into per-base TSV (PySAM).
6. **Unique mapping** – SNP-focused mapper identifies uniquely supported references.
7. **Output creation** – Writes gene-specific FASTQ/BAM outputs.

## Output Layout

```
output/
├── ref_seq_msa.aln                # MAFFT alignment
├── ref_seq_msa.tsv                # Reference→MSA coordinate map
├── mapped_reads.sam               # Aligner output (if generated)
├── mapped_reads.tsv               # Read→reference TSV
├── unique_mappings.tsv            # Final unique mapping table
├── fastq/                         # Gene-specific FASTQs
└── bam/                           # Gene-specific BAMs (+ .bai)
```

## Development Tips

- Run the pipeline from the repo root (`python mapper.py ...`).
- To add new stages, place modules under `src/pipeline/` and import from `mapper.py`.
- Keep the Conda spec (`mapper_env.yml`) in sync when adding binary dependencies.
