# ParaDISM: Paralog Disambiguating Mapper

ParaDISM has been published in BMC Bioinformatics: [link to published paper](https://doi.org/10.1186/s12859-026-06701-y)

ParaDISM is a read-mapping pipeline with optional refinement for highly
homologous genomic regions. It aligns short reads to a multi-sequence
reference, uses gene-informative positions in the multiple-sequence alignment
to assign reads, and writes gene-specific FASTQ/BAM outputs.

## Installation

```bash
conda env create -f environment.yml
conda activate paradism
```

## Quick Demo

Run the included synthetic demo. It does not require downloads or private data:

```bash
bash demo/run_demo.sh
```

The demo uses committed synthetic reads from a PKD1/pseudogene reference:

- `demo/ref.fa`
- `demo/pkd1_sim_R1.fq`
- `demo/pkd1_sim_R2.fq`

and writes checked outputs to `demo/output/`. The deterministic result contains
9 correctly assigned pairs, 11 unassigned pairs, and no incorrect assignments;
the demo prints and validates these counts and writes
`demo_assignment_summary.tsv`. The demo permits up to three iterations and
converges when no variants are found as iteration 2 begins.

To regenerate the tiny demo reads from `demo/ref.fa`, run:

```bash
bash demo/generate_demo_reads.sh
```

## Usage

### Interactive Mode

Run without arguments to launch the guided CLI:

```bash
python paradism.py
```

Interactive mode scans the current directory for input files and guides you
through file selection, aligner choice, output directory, run parameters, and
optional settings.
If the current directory has no FASTQ files, ParaDISM searches nested
directories for reads. It auto-selects a single reads directory, or prompts you
to choose when multiple candidates are found. Use
`--input-dir` to start from a different directory, and provide `--reference`
when the reference FASTA is stored outside that directory.
When run interactively, reference selection includes FASTA files in the selected
reads directory and FASTA files under the original starting directory.

### Non-Interactive Mode

The reference input file should contain the sequence(s) you want ParaDISM to distinguish in FASTA format.
The FASTA record names are used as the gene/local-contig names in the outputs.

```bash
python paradism.py \
  --read1 reads_R1.fq \
  --read2 reads_R2.fq \
  --reference ref.fa \
  --aligner bowtie2 \
  --threads 8 \
  --workers 8 \
  --iterations 1 \
  --output-dir output
```

Flags:

- `--read1 READ1`: R1 FASTQ file, or single-end FASTQ file.
- `--read2 READ2`: R2 FASTQ file for paired-end mode.
- `--reference REFERENCE`: FASTA reference containing the gene/paralog
  sequence(s) to distinguish.
- `--sam SAM`: Existing SAM file to use instead of running a new alignment.
- `--aligner ALIGNER`: Short-read base aligner, one of `bowtie2`, `bwa-mem2`,
  or `minimap2`.
- `--threads THREADS`: Threads passed to the base aligner.
- `--workers WORKERS`: Worker processes for ParaDISM read assignment after
  alignment.
- `--compress-intermediate-sam`: After assignment, validate and replace each
  `mapped_reads.sam` with `mapped_reads.bam`. This reduces retained disk usage
  but not the temporary peak while the aligner is writing SAM.
- `--output-dir OUTPUT_DIR`: Output directory.
- `--prefix PREFIX`: Prefix for output files.
- `--iterations N`: Number of ParaDISM runs; `1` is a single run, larger values
  enable iterative refinement.
- `--n_anchors N`: Minimum number of distinct gene-unique C1 anchor positions
  required for read assignment.
- `--input-dir INPUT_DIR`: Reads directory scanned by interactive mode.
- `--threshold THRESHOLD`: Minimum alignment score threshold or Bowtie2 score
  function.
- `--min-alternate-count N`: Minimum alternate allele count for FreeBayes during
  iterative variant calling.
- `--add-quality-filters`: Apply quality filters during iterative variant
  calling.
- `--qual-threshold N`: Minimum QUAL score when quality filters are enabled.
- `--dp-threshold N`: Minimum depth when quality filters are enabled.
- `--af-threshold F`: Minimum allele frequency when quality filters are enabled.

## Output Layout

```text
output/
├── .paradism_complete               # written only after successful completion
├── <prefix>_pipeline_<time>.log     # stage log and command output
├── iteration_1/
│   ├── mapped_reads.sam             # initial direct-alignment records
│   ├── <prefix>_fastq/              # assignments made in this iteration
│   └── <prefix>_bam/                # per-contig BAMs for this iteration
├── iteration_<n>/                   # present when refinement is requested
│   ├── mapped_reads.sam             # realigned previously unresolved reads
│   └── variant_calling/             # variants and updated refinement reference
└── final_outputs/
    ├── ref_seq_msa.aln              # MSA used for the final assignments
    ├── <prefix>_fastq/              # assigned reads, one FASTQ per contig
    ├── <prefix>_bam/                # sorted/indexed BAMs on the input reference
    └── <prefix>_none/               # unresolved FASTQs and inspection BAM
```

The FASTQ and BAM filenames contain the corresponding FASTA record name.
`<prefix>_none/` retains reads that did not uniquely support one reference
sequence; these reads are not silently forced onto a paralog. Iteration
directories are diagnostic intermediates, whereas `final_outputs/` is the
stable location intended for downstream analysis. Do not treat a run as
complete unless `.paradism_complete` exists.

The initial and refinement alignments are retained so assignments can be
audited and the base-aligner result can be evaluated. By default they remain
as SAM. With `--compress-intermediate-sam`, each consumed SAM is validated and
atomically replaced by a BAM with the same stem; the original SAM is preserved
if conversion fails. See `benchmark/giab/README.md` before running the full
HG002 benchmark.

## Liftover

ParaDISM reports VCF/BED positions on the local contigs in `ref.fa`. Use
`liftover` only when you need to convert those local coordinates back to
chromosomal coordinates:

```bash
python paradism.py liftover \
  --bed benchmark/references/pkd1_exons.bed \
  --ref benchmark/references/pkd1_panel.fa \
  --positions benchmark/references/pkd1_panel_positions.txt \
  --output lifted_exons.bed
```

The `--positions` file maps each FASTA contig name to the exact chromosomal
interval it was cut from, including any flanking sequence, and its strand
(`-1` when the contig is the reverse complement of the interval). Chromosome
names are written to the output as given:

```text
PKD1 chr16:2088108-2136498:-1
PKD1P1 chr16:16309741-16334790:1
```

`liftover` stops with an error if a contig length in `--ref` differs from its
interval. With `--genome genome.fa`, it also checks each contig sequence
against the genome and writes chromosome lengths to the VCF header. On
minus-strand contigs, VCF alleles are reverse complemented and indels are
left-aligned again with `bcftools norm`.

`liftover` does not convert between genome assemblies. For that, use a chain
file with a tool such as transanno, CrossMap or Picard `LiftoverVcf`.

## Reproducible Benchmarks

- `benchmark/references/`: committed FASTA references for the included
  gene/paralog groups.
- `benchmark/simulation/`: synthetic read simulation with `dwgsim`, followed
  by ParaDISM/Bowtie2/BWA-MEM2/minimap2 read-assignment evaluation. It writes
  per-seed metrics, aggregate read-mapping CSVs, and timing data.
- `benchmark/giab/`: public HG002/GIAB benchmark scripts that evaluate SNP
  calls against GIAB truth.

Private data workflows are not part of the public benchmark.
