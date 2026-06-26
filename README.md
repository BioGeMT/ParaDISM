# ParaDISM: Paralog Disambiguating Mapper

ParaDISM is a read-mapping and refinement workflow for highly homologous
genomic regions. It aligns short or long reads to a multi-sequence reference,
uses gene-informative positions in the multiple-sequence alignment to assign
reads, and writes gene-specific FASTQ/BAM outputs.

## Installation

The supported installation is the pinned conda/mamba environment. ParaDISM
uses external bioinformatics binaries in addition to Python packages, so this
repository does not use `uv` as the primary reproducibility mechanism.

```bash
mamba env create -f environment.yml
mamba activate paradism
```

`conda` can be used instead of `mamba`:

```bash
conda env create -f environment.yml
conda activate paradism
```

The environment pins ParaDISM's non-Python tools (`bowtie2`,
`bwa-mem2`, `minimap2`, `samtools`, `bcftools`, `freebayes`, `mafft`, and
`dwgsim`) as well as the required Python packages.

The repository also includes a GitHub Actions smoke test on Ubuntu 24.04 that
creates this environment and runs the demo plus the small simulation check.

## Quick Demo

Run the included synthetic demo first. It does not require downloads or private
data:

```bash
bash demo/run_demo.sh
```

The demo uses:

- `demo/ref.fa`
- `demo/reads_R1.fq`
- `demo/reads_R2.fq`

and writes checked outputs to `demo/output/`.

For a slightly larger synthetic smoke test that also exercises `dwgsim` and
the benchmark aggregation code, run:

```bash
bash benchmark/simulation/run_simulation.sh \
  --group hba_pair \
  --reference benchmark/references/hba_pair.fa \
  --seeds 1 \
  --reads-per-seed 20 \
  --aligners bowtie2 \
  --iterations 1 \
  --threads 1 \
  --out results/simulation/smoke_hba_pair
```

## Basic Usage

```bash
python paradism.py \
  --read1 reads_R1.fq \
  --read2 reads_R2.fq \
  --reference ref.fa \
  --aligner bowtie2 \
  --threads 4 \
  --iterations 2 \
  --output-dir output
```

Supported aligners are `bowtie2`, `bwa-mem2`, and `minimap2`. For minimap2,
also provide `--minimap2-profile`, for example `--minimap2-profile short`.

Use `--iterations 1` for a single ParaDISM run. Larger values enable iterative
refinement of reads initially assigned to `NONE`.

## Output Layout

```text
output/
├── mapped_reads.sam                  # present for one-iteration runs
├── iteration_1/mapped_reads.sam       # present for iterative runs
└── final_outputs/
    ├── <prefix>_fastq/               # gene-specific FASTQs
    ├── <prefix>_bam/                 # gene-specific sorted BAMs
    └── <prefix>_none/                # optional unresolved-read outputs
```

## Reproducible Benchmarks

Reviewer-facing benchmark workflows are under `benchmark/`:

- `benchmark/references/`: committed FASTA references for the included
  paralog groups, plus a helper for building additional Ensembl-based FASTAs.
- `benchmark/simulation/`: synthetic read simulation with `dwgsim`, followed
  by ParaDISM/direct-aligner read-assignment evaluation.
- `benchmark/giab/`: public HG002/GIAB benchmark scripts. These require large
  public GIAB read inputs and are intentionally separate from the quick demo.

Private or unpublished-data workflows are not part of the public benchmark
path.
