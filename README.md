# PKD1 Read Mapping Pipeline

This pipeline performs unique read mapping and generates gene-specific outputs in FASTQ and BAM formats. It processes paired-end sequencing data against reference sequences using multiple alignment algorithms.

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
- bowtie2 (default aligner)
- bwa-mem2 (optional)
- minimap2 (optional)
- pandas
- samtools

## Directory Structure

```
.
├── scripts/
│   ├── output.py
│   ├── read_2_gene.py
│   ├── ref_2_msa.py
│   └── mapper_algo_snp_only.py
│
├── mapper.sh
├── mapper_interactive.sh
├── simulated_r1.fq
├── simulated_r2.fq
└── mapper_env.yml

```

## Usage

First run the read simulator notebook to generate simulated paired reads. This produces two fq files.

### Interactive Mode (Recommended)
Run without arguments to launch interactive mode with guided prompts:

```bash
./mapper.sh
```

The interactive mode will guide you through:
- Selecting read files (R1 and R2)
- Choosing a reference sequence
- Option to use an existing SAM file (skips alignment)
- Picking an aligner (Bowtie2, BWA-MEM2, or minimap2)
- Setting thread count
- Configuring minimap2 profile (if applicable)

### Command-Line Mode
```bash
./mapper.sh --read1 <forward_reads.fq> --read2 <reverse_reads.fq> --reference <reference.fasta> [OPTIONS]
```

### Arguments
- `--read1`: Forward reads FASTQ file (required)
- `--read2`: Reverse reads FASTQ file (required)
- `--reference`: Reference sequences in FASTA format (required)
- `--sam`: Pre-existing SAM/BAM alignment file (optional, skips alignment step)
- `--aligner`: Alignment tool to use: `bowtie2` (default), `bwa-mem2`, or `minimap2`
- `--threads`: Number of threads for alignment (default: 4)
- `--minimap2-profile`: Profile for minimap2: `short` (default), `pacbio`, or `nanopore`

### Examples

**Using default Bowtie2 aligner:**
```bash
./mapper.sh --read1 simulated_r1.fq --read2 simulated_r2.fq --reference ref.fa
```

**Using BWA-MEM2 with 8 threads:**
```bash
./mapper.sh --read1 simulated_r1.fq --read2 simulated_r2.fq --reference ref.fa --aligner bwa-mem2 --threads 8
```

**Using minimap2 for short reads:**
```bash
./mapper.sh --read1 simulated_r1.fq --read2 simulated_r2.fq --reference ref.fa --aligner minimap2 --minimap2-profile short
```

**Using minimap2 for PacBio data:**
```bash
./mapper.sh --read1 pacbio_r1.fq --read2 pacbio_r2.fq --reference ref.fa --aligner minimap2 --minimap2-profile pacbio --threads 16
```

**Using your own alignment (skip alignment step):**
```bash
./mapper.sh --read1 simulated_r1.fq --read2 simulated_r2.fq --reference ref.fa --sam my_custom_alignment.sam
```
Note: When `--sam` is provided, the alignment step is skipped and `--aligner`, `--threads`, and `--minimap2-profile` options are ignored.

**SAM File Requirements:**
Your SAM file must include:
- Valid SAM header
- Aligned reads (not all unmapped)
- MD tags (use `samtools calmd` if missing)

  


## Pipeline Steps

1. **MSA Generation**: Uses MAFFT to create multiple sequence alignment of reference sequences
2. **Reference Indexing**: Builds index for the selected aligner (Bowtie2, BWA-MEM2, or minimap2)
3. **Read Alignment**: Aligns paired-end reads to reference sequences using the selected aligner
4. **Reference to MSA Mapping**: Maps reference sequences to MSA positions using BioPython
5. **Read to Gene Mapping**: Processes SAM alignments to extract read-to-gene mappings using pysam
6. **Unique Mapping Algorithm**: Identifies uniquely mapped reads using dynamic programming (SNP-optimized)
7. **Output Generation**: Creates gene-specific FASTQ and BAM files from unique mappings

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
