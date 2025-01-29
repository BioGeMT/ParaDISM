# PKD1 Read Mapping Pipeline

This pipeline performs unique read mapping, and generates gene-specific outputs in FASTQ and/or BAM format. It processes paired-end sequencing data against reference sequences.

## Prerequisites

### Environment Setup
To create a conda environment with all the required dependencies run the following:

```bash
conda env create -f mapper_env.yml
conda activate mapper_env
```

Required dependencies:
- Python 3.10
- pysam
- biopython
- mafft
- bowtie2
- pandas
- minimap2
- samtools

## Directory Structure

```
.
├── scripts/
│   ├── bam_output.py
│   ├── find_mappings.py
│   ├── fq_output.py
│   ├── read_2_gene.py
│   ├── reads_2_msa.py
│   ├── ref_2_msa.py
│   └── refine.py
├── mapper.sh
├── mapper_env.yml
└── README.md
```

## Usage

### Basic Command
```bash
./mapper.sh -r1 <forward_reads.fq> -r2 <reverse_reads.fq> -ref <reference.fasta> [-out <fastq,bam>]
```

### Arguments
- `-r1`: Forward reads FASTQ file (required)
- `-r2`: Reverse reads FASTQ file (required)
- `-ref`: Reference sequences in FASTA format (required)
- `-out`: Output format(s), comma-separated (optional)
  - `fastq`: Generate gene-specific FASTQ files
  - `bam`: Generate gene-specific BAM files
  - Default: both formats if not specified

## Pipeline Steps

1. **MSA Generation**: Uses MAFFT to create multiple sequence alignment of reference sequences
2. **Read Mapping**: Aligns reads to references using Bowtie2
3. **Reference to MSA Mapping**: Maps reference sequences to MSA positions
4. **Read to Gene Mapping**: Processes SAM alignments to map reads to genes
5. **MSA Integration**: Maps reads to MSA coordinates
6. **Refinement**: Refines read mappings
7. **Unique Mapping**: Identifies uniquely mapped reads
8. **Output Generation**: Creates gene-specific files in requested format

## Output Structure

```
output/
├── ref_seq_msa.aln                 # MAFFT MSA output
├── PKD_index*                      # Bowtie2 index files
├── mapped_reads.sam                # Initial Bowtie2 alignment
├── ref_seq_msa.tsv                 # Reference to MSA mapping
├── mapped_reads.tsv                # Read to reference mapping
├── reads/                          # Read to MSA mapping files
├── results/                        # Refinement results for each read
├── unique_mappings.tsv             # Unique mapping results
├── unique_mappings_fastq_files/    # Gene-specific FASTQ files (if requested)
└── unique_mappings_bam_files/      # Gene-specific BAM files (if requested)
```
