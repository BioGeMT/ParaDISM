# Homology Mapper Pipeline

Read-mapping and refinement workflow for highly homologous genomic regions. The pipeline aligns paired-end reads, maps alignments back to multiple-sequence alignments, refines the mapping for uniquely supported reads and reference sequences, and produces gene-specific FASTQ/BAM outputs.

## Prerequisites

### Environment Setup
To create a conda environment with all required dependencies, run:

```bash
conda env create -f mapper_env.yml
conda activate mapper_env
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

> **Note:** Interactive mode scans the current working directory for input files. Place your inputs in the project root (or run the command from the directory that contains them) so they appear in the selection lists.

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
