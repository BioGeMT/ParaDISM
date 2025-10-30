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

Or specify a custom input directory:

```bash
python mapper.py --input-dir /path/to/data
```

You'll be prompted to:
- Select sequencing mode (Paired‑End or Single‑End)
- Pick FASTQ file(s) accordingly (R1/R2 or a single FASTQ)
- Choose reference FASTA
- Optionally use an existing SAM (skips alignment)
- Pick an aligner (`bowtie2`, `bwa-mem2`, `minimap2`)
- Configure threads and minimap2 profile (if using minimap2)

 **Note:** By default, interactive mode scans the current working directory for input files. Use `--input-dir` to specify a different directory to scan. Place your inputs in the specified directory (or current directory if not specified) so they appear in the selection lists.

### Argument-Driven Mode
```bash
python mapper.py --read1 <forward_reads.fq> \
                 [--read2 <reverse_reads.fq>] \
                 --reference <reference.fa> \
                 [--aligner ALIGNER] \
                 [--threads N] \
                 [--minimap2-profile PROFILE] \
                 [--sam existing.sam] \
                 [--output-dir OUTPUT] \
                 [--prefix PREFIX]
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
- `--prefix`: Prefix for output files (default: derived from output directory name)

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

# Using custom output directory (prefix auto-derived)
python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa \
  --output-dir sample_001
  # Output files will be prefixed with "sample_001_"

# Explicitly specify a different prefix
python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa \
  --output-dir results --prefix custom_prefix
  # Output files will be prefixed with "custom_prefix_"
```

### SAM Requirements
When skipping alignment via `--sam`, ensure the SAM:
- Contains a valid header
- Has at least one mapped read
- Includes MD tags (`samtools calmd` can add them). The pipeline requests MD tags during alignment (e.g., `minimap2 --MD`).

## Output Layout

By default, output files are prefixed with the output directory name. For example, if `--output-dir SAMPLE_NAME`, all final outputs will be prefixed with `SAMPLE_NAME_`.

```
output_dir/                                   # Output directory
├── prefix_pipeline_YYYYMMDD_HHMMSS.log      # Run log (sections + streamed output)
├── prefix_unique_mappings.tsv               # Final unique mapping table
├── prefix_fastq/                            # Gene-specific FASTQs directory
│   ├── prefix_gene1.fq                      # Gene 1 reads
│   ├── prefix_gene2.fq                      # Gene 2 reads
│   ├── prefix_gene3.fq                      # Gene 3 reads
│   └── ...                                  # Additional genes
└── prefix_bam/                              # Gene-specific BAMs directory
    ├── prefix_gene1.sorted.bam              # Gene 1 alignments
    ├── prefix_gene1.sorted.bam.bai          # BAM index
    ├── prefix_gene2.sorted.bam              # Gene 2 alignments
    ├── prefix_gene2.sorted.bam.bai
    └── ...                                  # Additional genes
```

**Note:** Intermediate files (MSA, SAM, indices) are created during processing but cleaned up automatically at the end. Only the final outputs listed above are retained.

### Customizing the Prefix

You can control the output file prefix using the `--prefix` option:

```bash
# Auto-derive prefix from output directory name (default)
python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa --output-dir sample_001

# Explicitly specify custom prefix
python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa \
  --output-dir results --prefix custom_prefix

# Use empty prefix (no prefix on output files)
python mapper.py --read1 r1.fq --read2 r2.fq --reference ref.fa \
  --output-dir output --prefix ""
```
