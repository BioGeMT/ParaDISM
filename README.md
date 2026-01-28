# ParaDISM: Paralog Disambiguating Mapper

A read-mapping and refinement workflow for highly homologous genomic regions. The pipeline aligns reads (paired-end or single-end), maps alignments back to a multiple-sequence alignment (MSA), refines to uniquely supported mappings, and produces gene-specific FASTQ/BAM outputs.

## Prerequisites

### Environment Setup
To create a conda environment with all required dependencies, run:

```bash
conda env create -f paradism_env.yml
conda activate paradism_env
```

## Usage

### Interactive Mode
Run without arguments to launch the guided CLI:

```bash
python paradism.py
```

Or specify a custom input directory:

```bash
python paradism.py --input-dir /path/to/data
```

You'll be prompted to:
- Select sequencing mode (Paired-End or Single-End)
- Pick FASTQ file(s) accordingly (R1/R2 or a single FASTQ)
- Choose reference FASTA
- Optionally use an existing SAM (skips alignment)
- Pick an aligner (`bowtie2`, `bwa-mem2`, `minimap2`)
- Configure threads and minimap2 profile (if using minimap2)
- Set alignment score threshold
- Configure iterations for iterative refinement
- Set variant calling options (when iterations > 1):
  - Minimum alternate allele count
  - Quality filters (QUAL, DP, AF thresholds)

**Note:** By default, interactive mode scans the current working directory for input files. Use `--input-dir` to specify a different directory to scan.

### Argument-Driven Mode
```bash
python paradism.py --read1 <forward_reads.fq> \
                 [--read2 <reverse_reads.fq>] \
                 --reference <reference.fa> \
                 [--aligner ALIGNER] \
                 [--threads N] \
                 [--minimap2-profile PROFILE] \
                 [--sam existing.sam] \
                 [--output-dir OUTPUT] \
                 [--prefix PREFIX] \
                 [--iterations N] \
                 [--threshold THRESHOLD] \
                 [--min-alternate-count N] \
                 [--add-quality-filters] \
                 [--qual-threshold N] \
                 [--dp-threshold N] \
                 [--af-threshold F]
```

**Required:**
- `--read1`: R1 FASTQ file (or single-end reads)
- `--reference`: Reference FASTA

**Optional - General:**
- `--read2`: R2 FASTQ file (for paired-end mode)
- `--aligner`: `bowtie2` (default), `bwa-mem2`, or `minimap2`
- `--threads`: Number of alignment threads (default: 4)
- `--minimap2-profile` (required with minimap2): one of
  - `short` (short-read Illumina)
  - `pacbio-hifi`
  - `pacbio-clr`
  - `ont-q20`
  - `ont-standard`
- `--sam`: Existing alignment (skips alignment stage)
- `--output-dir`: Destination directory (default: `./output`)
- `--prefix`: Prefix for output files (default: derived from output directory name)
- `--threshold`: Alignment score threshold (aligner-specific). For example: `G,40,40` (bowtie2) or `240` (bwa-mem2/minimap2).

**Optional - Iterative Refinement:**
- `--iterations`: Total ParaDISM runs (default: 1). Use >1 to enable iterative refinement.
- `--min-alternate-count`: Minimum alternate allele count for FreeBayes variant calling (default: 5)
- `--add-quality-filters`: Enable quality filtering during variant calling
- `--qual-threshold`: Minimum QUAL score for quality filtering (default: 20)
- `--dp-threshold`: Minimum depth (DP) for quality filtering (default: 10)
- `--af-threshold`: Minimum allele frequency (AF) for quality filtering (default: 0.05)

**Paired-end vs single-end:**
- For paired-end runs, provide both `--read1` and `--read2`.
- If `--read2` is omitted, the pipeline runs in single-end mode.

### Examples

```bash
# Paired-end short reads, default aligner (Bowtie2)
python paradism.py --read1 r1.fq --read2 r2.fq --reference ref.fa

# Single-end long reads with minimap2 (PacBio HiFi)
python paradism.py --read1 hifi.fq --reference ref.fa \
  --aligner minimap2 --minimap2-profile pacbio-hifi

# Using custom output directory (prefix auto-derived)
python paradism.py --read1 r1.fq --read2 r2.fq --reference ref.fa \
  --output-dir sample_001

# Two pipeline runs with iterative refinement
python paradism.py --read1 r1.fq --read2 r2.fq --reference ref.fa \
  --iterations 2

# Iterative refinement with custom variant calling parameters
python paradism.py --read1 r1.fq --read2 r2.fq --reference ref.fa \
  --iterations 5 \
  --min-alternate-count 5 \
  --add-quality-filters \
  --qual-threshold 20 \
  --dp-threshold 10 \
  --af-threshold 0.05
```

### SAM Requirements
When skipping alignment via `--sam`, ensure the SAM:
- Contains a valid header
- Has at least one mapped read
- Includes MD tags (`samtools calmd` can add them)

## Output Layout

By default, output files are prefixed with the output directory name. For example, if `--output-dir SAMPLE_NAME`, all final outputs will be prefixed with `SAMPLE_NAME_`.

```
output_dir/
├── prefix_pipeline_YYYYMMDD_HHMMSS.log      # Run log
├── iteration_1/                             # Iteration outputs (when iterations > 1)
├── iteration_2/                             # Additional iterations
└── final_outputs/
    ├── prefix_fastq/                        # Gene-specific FASTQs
    │   ├── prefix_gene1.fq
    │   ├── prefix_gene2.fq
    │   └── ...
    └── prefix_bam/                          # Gene-specific BAMs
        ├── prefix_gene1.sorted.bam
        ├── prefix_gene1.sorted.bam.bai
        └── ...
```

**Note:** Intermediate files (MSA, SAM, indices) are created during processing but cleaned up automatically. Only final outputs are retained.

## Iterative Refinement

Set `--iterations` to run ParaDISM multiple times. The first run produces mappings and per-gene BAMs. Each refinement iteration:
1. Calls variants from previously mapped (non-NONE) reads using FreeBayes
2. Optionally applies quality filters (QUAL, DP, AF thresholds)
3. Updates the reference with called variants
4. Re-aligns only reads that were labeled `NONE` against the updated reference
5. Merges results so prior successful mappings remain unchanged

The loop stops early if:
- No reads remain to rescue
- No variants to apply
- No reads were reassigned in the latest iteration

## Analysis and Experiments

See individual README files in each directory:
- `giab_benchmark/` - GIAB HG002 benchmarking
- `gnaq_analysis/` - GNAQ samples analysis
- `hts_analysis/` - HTS clinical samples
- `simulation/` - Simulated reads analysis
