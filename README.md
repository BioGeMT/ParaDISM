# ParaDISM: Paralog Disambiguating Mapper

ParaDISM is a read-mapping workflow with optional refinement for highly
homologous genomic regions. It aligns short or long reads to a multi-sequence
reference, uses gene-informative positions in the multiple-sequence alignment
to assign reads, and writes gene-specific FASTQ/BAM outputs.

## Installation

The supported installation is the pinned conda environment. ParaDISM uses
external bioinformatics binaries in addition to Python packages, so this
repository does not use `uv` as the primary reproducibility mechanism.

```bash
conda env create -f environment.yml
conda activate paradism
```

The environment pins ParaDISM's non-Python tools (`bowtie2`,
`bwa-mem2`, `minimap2`, `samtools`, `bcftools`, `freebayes`, `mafft`, and
`dwgsim`) as well as the required Python packages.

## Quick Demo

Run the included synthetic demo first. It does not require downloads or private
data:

```bash
bash demo/run_demo.sh
```

The demo uses:

- `demo/ref.fa`
- 20 generated paired-end read pairs under `demo/generated_reads/`

and writes checked outputs to `demo/output/`.

## Basic Usage

```bash
python paradism.py \
  --read1 reads_R1.fq \
  --read2 reads_R2.fq \
  --reference ref.fa \
  --aligner bowtie2 \
  --threads 4 \
  --workers 2 \
  --iterations 1 \
  --output-dir output
```

Supported aligners are `bowtie2`, `bwa-mem2`, and `minimap2`. For minimap2,
also provide `--minimap2-profile`, for example `--minimap2-profile short`.

Use `--iterations 1` for a single ParaDISM run. Larger values enable iterative
refinement of reads initially assigned to `NONE`.

Use `--n_anchors N` to require at least `N` distinct gene-unique C1 positions
for read assignment. The default is `--n_anchors 1`, matching the original
behavior.

Use `--workers N` to process read assignment in parallel after alignment. This
is separate from `--threads`, which controls the underlying aligner.

## Output Layout

```text
output/
├── iteration_1/
│   └── mapped_reads.sam              # initial direct-alignment SAM
└── final_outputs/
    ├── <prefix>_fastq/               # gene-specific FASTQs
    ├── <prefix>_bam/                 # gene-specific sorted BAMs
    └── <prefix>_none/                # optional unresolved-read outputs
```

## Liftover

ParaDISM outputs use gene-local contig coordinates. To convert VCF or BED files
back to chromosomal coordinates, use the liftover subcommand:

```bash
python paradism.py liftover \
  --bed benchmark/references/pkd1_exons.bed \
  --positions positions.txt \
  --output lifted_exons.bed
```

The `--positions` file must define the chromosomal interval and strand for
each gene/reference contig being lifted. Each line should end with
`CHR:START-END:STRAND`, for example:

```text
PKD1 16:2088708-2135898:1
PKD1P1 16:16310341-16334190:1
```

## Reproducible Benchmarks

Reviewer-facing benchmark workflows are under `benchmark/`:

- `benchmark/references/`: committed FASTA references for the included
  paralog groups.
- `benchmark/simulation/`: synthetic read simulation with `dwgsim`, followed
  by ParaDISM/direct-aligner read-assignment evaluation.
- `benchmark/giab/`: public HG002/GIAB benchmark scripts. These require large
  public GIAB read inputs and are intentionally separate from the quick demo.

Private or unpublished-data workflows are not part of the public benchmark
path.
