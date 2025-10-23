# Homologous-Region Mapping Pipeline

Read-mapping and refinement workflow for highly homologous genomic regions. The pipeline aligns reads (paired-end or single-end), maps alignments back to a multiple‑sequence alignment (MSA), refines to uniquely supported mappings, and produces gene‑specific FASTQ/BAM outputs.

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
- Select sequencing mode (Paired‑End or Single‑End)
- Pick FASTQ file(s) accordingly (R1/R2 or a single FASTQ)
- Choose reference FASTA
- Optionally use an existing SAM (skips alignment)
- Pick an aligner (`bowtie2`, `bwa-mem2`, `minimap2`)
- Configure threads and minimap2 profile (if using minimap2)

 **Note:** Interactive mode scans the current working directory for input files. Place your inputs in the project root (or run the command from the directory that contains them) so they appear in the selection lists.

### Argument-Driven Mode
```bash
python mapper.py --read1 <forward_reads.fq> \
                 [--read2 <reverse_reads.fq>] \
                 --reference <reference.fa> \
                 [--aligner ALIGNER] \
                 [--threads N] \
                 [--minimap2-profile PROFILE] \
                 [--sam existing.sam] \
                 [--output-dir OUTPUT]
```

Required:
- `--read1`: FASTQ (R1 for paired‑end, or single‑end reads)
- `--reference`: Reference FASTA

Optional:
- `--read2`: R2 FASTQ for paired‑end mode
- `--aligner`: `bwa-mem2` (default), `bowtie2`, or `minimap2`
- `--threads`: Number of alignment threads (default: 4)
- `--minimap2-profile` (required with minimap2): one of
  - `short` (short‑read Illumina)
  - `pacbio-hifi`
  - `pacbio-clr`
  - `ont-q20`
  - `ont-standard`
- `--sam`: Existing alignment (skips alignment stage)
- `--output-dir`: Destination directory (default: `./output`)

Paired‑end vs single‑end
- For paired‑end runs, provide both `--read1` and `--read2`.
- If `--read2` is omitted, the pipeline runs in single‑end mode.

Examples:
```bash
# Paired-end, default aligner (BWA-MEM2)
python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa

# Single-end with minimap2 (PacBio HiFi)
python mapper.py --read1 hifi.fq --reference ref.fa \
  --aligner minimap2 --minimap2-profile pacbio-hifi

# Single-end with minimap2 (ONT Q20+)
python mapper.py --read1 ont_q20.fq --reference ref.fa \
  --aligner minimap2 --minimap2-profile ont-q20
```

### SAM Requirements
When skipping alignment via `--sam`, ensure the SAM:
- Contains a valid header
- Has at least one mapped read
- Includes MD tags (`samtools calmd` can add them). The pipeline requests MD tags during alignment (e.g., `minimap2 --MD`).

## Output Layout

```
output/
├── ref_seq_msa.aln                # MAFFT alignment
├── ref_seq_msa.tsv                # Reference→MSA coordinate map
├── mapped_reads.sam               # Aligner output (if generated)
├── mapped_reads.tsv               # Read→reference TSV
├── unique_mappings.tsv            # Final unique mapping table
├── fastq/                         # Gene-specific FASTQs
├── bam/                           # Gene-specific BAMs (+ .bai)
└── pipeline_YYYYMMDD_HHMMSS.log  # Run log (sections + streamed output)

Optional variant-calling outputs (if you run `var_calling.sh`):

```
output/
└── variant_calling/
    ├── freebayes/*.vcf
    ├── bcftools/*.vcf
    └── gatk/*.vcf
```

