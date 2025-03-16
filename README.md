# PKD1 Read Mapping Pipeline

This pipeline performs unique read mapping, and generates gene-specific outputs in FASTQ and BAM formats. It processes paired-end sequencing data against reference sequences.

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
- samtools
## Directory Structure

```
.
├── scripts/
│   ├── output.py
│   ├── read_2_gene.py 
│   ├── ref_2_msa.py
│   └── mapper_algo.py
│      
├── mapper.sh
├── simulated_r1.fq
├── simulated_r2.fq
└── mapper_env.yml

```

## Usage

First run the read simulator notebook to generate simulated paired reads. This produces two fq files.

### Basic Command
```bash
./mapper.sh -r1 <forward_reads.fq> -r2 <reverse_reads.fq> -ref <reference.fasta>
```

### Arguments
- `-r1`: Forward reads FASTQ file (required)
- `-r2`: Reverse reads FASTQ file (required)
- `-ref`: Reference sequences in FASTA format (required)

  


## Pipeline Steps

1. **MSA Generation**: Uses MAFFT to create multiple sequence alignment of reference sequences
2. **Read Mapping**: Aligns reads to references using Bowtie2
3. **Reference to MSA Mapping**: Maps reference sequences to MSA positions
4. **Read to Gene Mapping**: Processes SAM alignments to map reads to genes
5. **MSA Integration**: Maps reads to MSA coordinates
6. **Refinement**: Refines read mappings
7. **Unique Mapping**: Identifies uniquely mapped reads
8. **Output Generation**: Creates gene-specific files in fastq and bam formats

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
├── fastq/                          # Gene-specific FASTQ files 
└── bam/                            # Gene-specific BAM files 
```
