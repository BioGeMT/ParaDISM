# GIAB Benchmark

This directory contains scripts for public HG002/GIAB benchmarking.

The workflow analyzes PKD1/PKD1P SNP calls from ParaDISM-refined BAMs versus
the base Bowtie2 BAM, using GIAB HG002 truth as the benchmark.

## Inputs

- public HG002 Illumina FASTQ shards, downloaded by
  `download_giab_hg002_reads.sh`
- GIAB HG002 truth VCF and index
- GIAB benchmark BED
- `benchmark/references/pkd1_panel.fa`

## Resource requirements

The full public HG002 download is approximately 322.5 GB of compressed FASTQ
data. Plan for at least 1 TiB of free disk space for the download, the initial
SAM alignment, ParaDISM FASTQ/BAM outputs, and variant-calling intermediates.
The download script retains the 68 source shards and creates merged R1/R2
files, so two compressed copies of the reads coexist unless the source shards
are archived or removed after the merged files have been verified.

The GIAB launcher enables `--compress-intermediate-sam`. ParaDISM therefore
validates and replaces each consumed `mapped_reads.sam` with
`mapped_reads.bam`, while GIAB post-processing accepts either format. The raw
SAM still exists while the aligner is producing it, so compression reduces
retained disk usage rather than the peak space required during alignment.

This is a long-running whole-genome benchmark, not a quick demo. Alignment and
read assignment can each take hours; reserve a multi-hour or longer batch
window and use a persistent session or scheduler. Exact runtime depends on CPU,
available memory, storage throughput, worker count, and system load. The
launcher prints the measured FASTQ footprint, available output-disk space, a
warning when less than 1 TiB is available, elapsed time for each pipeline
stage, and the growing SAM size during alignment.

For a quick installation and correctness check, use `bash demo/run_demo.sh`
instead. The committed 20-pair demo normally completes in under one minute on
the tested macOS and Linux systems and validates its expected assignments.

### Observed runtimes

Archived full-depth HG002 measurements on a Linux server with two Intel Xeon
Gold 6342 CPUs provide practical reference points; they are not guaranteed
runtimes for other systems.

- A targeted `PKD1`/`PKD1P` family workflow using pre-extracted reads,
  Bowtie2, 16 threads, a single assignment process, and early convergence took
  14 min 54 s including FASTQ preparation, ParaDISM, variant calling,
  filtering, and final callset generation.
- Sequential targeted processing of all 13 configured families took
  2 h 21 min 25 s. Extracting all family-specific read sets from the indexed
  HG002 BAM was a separate approximately 1 min 7 s preparation step.
- An archived run that supplied the complete HG002 FASTQs directly to the
  `PKD1`/`PKD1P` reference took 16 h 33 min 25 s and converged after five
  iterations. Its iteration-1 SAM was 497,903,998,666 bytes (463.7 GiB).
- A fresh one-iteration reproduction with the complete HG002 FASTQs, 8
  alignment threads, and 2 assignment workers took 26 h 9 min 51 s and had a
  peak resident set size of 78.5 GiB. Subsequent variant calling and filtering
  took 1 h 23 min 9 s and had a peak resident set size of 6.7 GiB. These are
  single wall-clock measurements from a shared server.

The targeted-read workflow is therefore the practical reference for routine
family analysis. Direct raw-WGS processing requires a substantially longer
batch window and enough temporary space for the uncompressed SAM.

## Usage

From the repository root, after installing `environment.yml`:

```bash
bash benchmark/giab/prepare_giab_truth.sh
bash benchmark/giab/download_giab_hg002_reads.sh
bash benchmark/giab/run_giab.sh --threads 8 --workers 2 --iterations 10

RUN_DIR=benchmark/giab/giab_hg002_output_bowtie2_G60_min5_qfilters
OUT_DIR=benchmark/giab/vcf_out_full
bash benchmark/giab/run_variant_calling_to_vcf_out.sh \
  --run-dir "$RUN_DIR" \
  --out-dir "$OUT_DIR" \
  --threads 8
```

For a smaller 10x-style read subset, run `downsample_hg002_fastq.sh` after
downloading the HG002 shards, then pass the downsampled read directory to
`run_giab.sh` with `--reads-dir`.

`--iterations 10` reproduces the refinement setting used for the manuscript;
the pipeline stops early if no additional reads can be rescued. Use
`--iterations 1` only when a single, unrefined ParaDISM assignment pass is the
intended comparison.

### Interrupted runs

A run is complete only when `.paradism_complete` and the expected final BAMs
are present. If an interrupted run leaves a nonempty output directory,
`run_giab.sh` exits instead of treating that directory as complete or mixing a
new run with partial files. Preserve or move the partial directory for
diagnosis, then choose a new `--output-dir`. GIAB post-processing also exits
when the completion marker or expected BAMs are absent. It does not skip the
missing ParaDISM result and continue with only the base aligner.

## Output Layout

```text
benchmark/giab/
├── giab_hg002_reads/                 # downloaded or downsampled FASTQs
├── giab_hg002_vcf/                   # prepared GIAB truth and BED files
├── giab_hg002_output_*/              # ParaDISM run output
│   ├── .paradism_complete
│   ├── iteration_1/mapped_reads.bam  # compressed direct Bowtie2 alignment
│   ├── iteration_<n>/                # refinement intermediates
│   ├── final_outputs/                # final per-gene FASTQ/BAM outputs
│   └── variant_calling/              # ParaDISM and baseline VCFs
└── vcf_out*/                         # final evaluation CSV/JSON metrics
```

The post-processing script compares exact `CHROM`, `POS`, `REF`, and `ALT`
matches after restricting both call sets to simple biallelic A/C/G/T SNPs and
the GIAB benchmark regions. It maps GRCh38 truth positions to the exact panel
contig offsets and reverse-complements alleles for reverse-strand contigs. TP,
FP, and FN denote true-positive,
false-positive, and false-negative SNP calls. Precision is `TP/(TP+FP)`, recall
or sensitivity is `TP/(TP+FN)`, and F1 is their harmonic mean. True-negative
counts and specificity are not reported because the large number of invariant
sites makes specificity uninformative for this comparison.

Key evaluation files:

- `variant_calling_metrics.csv`: one row per method (`ParaDISM` and
  `BaseAligner`) with TP, FP, FN, precision, recall, and F1 for all seven
  PKD1/PKD1P contigs pooled together.
- `variant_calling_metrics.json`: the same pooled metrics in a nested,
  machine-readable structure, together with the truth-variant count.
- `confusion_matrices_overall.csv`: the pooled TP/FP/FN counts without derived
  rates, suitable for plotting a call-count matrix.
- `per_gene_metrics/per_gene_metrics.csv`: truth count, TP, FP, FN, precision,
  recall, and F1 for each method and each reference contig.
- `per_gene_metrics/per_gene_metrics.json`: the same per-contig rows in
  machine-readable form.
- `confusion_matrices_by_coverage.csv`: TP/FP/FN counts for each method below
  and above the selected read-depth threshold.
- `coverage_split_confusion_metrics.json`: the coverage threshold and how it
  was selected, plus pooled and coverage-stratified TP/FP/FN counts. By
  default, the threshold is the rounded pooled median truth-site depth across
  both methods.

## Reference results

The manuscript reports the following results for the archived, completed
workflows. These values are reference targets for interpreting a regenerated
run, not a substitute for checking the completion marker and the generated
metric files.

| Dataset | Method | TP | FP | FN | Precision | Recall | F1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Fresh full-depth HG002, one iteration | ParaDISM | 18 | 1 | 30 | 0.947 | 0.375 | 0.537 |
| Fresh full-depth HG002, one iteration | Bowtie2 baseline | 40 | 2 | 8 | 0.952 | 0.833 | 0.889 |
| Full-depth HG002 | ParaDISM | 17 | 5 | 31 | 0.773 | 0.354 | 0.486 |
| Full-depth HG002 | Bowtie2 baseline | 39 | 13 | 9 | 0.750 | 0.812 | 0.780 |
| Downsampled ~10x HG002 | ParaDISM | 11 | 0 | 37 | 1.000 | 0.229 | 0.373 |
| Downsampled ~10x HG002 | Bowtie2 baseline | 29 | 8 | 19 | 0.784 | 0.604 | 0.682 |

The default coverage threshold selected for the archived full-depth run was 48
reads; the threshold for the archived ~10x run was 6 reads. These results show
the intended precision-sensitivity trade-off: under the tested settings,
ParaDISM produced fewer false-positive calls but also recovered fewer true
variants than the Bowtie2 baseline.

The fresh one-iteration rows are a reproduction check of the reviewer-facing
command. The remaining rows are the archived manuscript workflows and use the
iteration settings described in the manuscript.

## Scripts

- `download_giab_hg002_reads.sh`: downloads public HG002 read inputs used by
  the benchmark.
- `downsample_hg002_fastq.sh`: creates a smaller read subset from downloaded
  HG002 shards.
- `prepare_giab_truth.sh`: prepares GIAB truth and benchmark regions in
  ParaDISM gene coordinates.
- `run_giab.sh`: runs ParaDISM on the GIAB read set.
- `run_variant_calling_to_vcf_out.sh`: calls variants from the completed
  ParaDISM run, filters to benchmarkable simple SNPs, and writes metrics.

The workflow uses public data and writes ParaDISM outputs, filtered VCFs, and
evaluation metrics under the selected output directories.
