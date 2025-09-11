#!/bin/bash

# Run VCFCompare per iteration for a method and plot SNV precision/recall
# Usage: ./plot_iterative_vcfcompare.sh <freebayes|bcftools> <iterations> <output_dir> <truth_vcf> [vcfcompare_dir]

set -euo pipefail

METHOD="${1:?method (freebayes|bcftools) required}"
ITERATIONS="${2:?number of iterations required}"
OUT_DIR="${3:-./output}"
TRUTH_VCF="${4:?ground truth VCF required}"
VCFCOMPARE_DIR="${5:-./VCFCompare}"

echo "=== VCFCompare Iterative: $METHOD ==="
echo "Iterations: $ITERATIONS"
echo "Output dir: $OUT_DIR"
echo "Truth VCF: $TRUTH_VCF"
echo "VCFCompare: $VCFCOMPARE_DIR"

PLOT_ROOT="plots_iterative_vcfcompare_${METHOD}"
mkdir -p "$PLOT_ROOT"

merge_vcfs() {
  # Merge multiple VCFs into one (keep first header, append records)
  # Args: output_vcf, file1, file2, ...
  local out_vcf="$1"; shift
  local first=
  : > "$out_vcf"
  for f in "$@"; do
    [[ -s "$f" ]] || continue
    if [[ -z "$first" ]]; then
      cat "$f" > "$out_vcf"
      first=1
    else
      awk 'BEGIN{FS=OFS="\t"} !/^#/' "$f" >> "$out_vcf"
    fi
  done
}

for ((i=1; i<=ITERATIONS; i++)); do
  ITER_DIR="$OUT_DIR/iter_$i"
  VC_DIR="$ITER_DIR/variant_calling"
  PLOT_DIR="$PLOT_ROOT/iter_$i"
  mkdir -p "$PLOT_DIR"

  if [[ ! -d "$VC_DIR" ]]; then
    echo "Warning: $VC_DIR not found, skipping iteration $i"
    continue
  fi

  echo "Processing iteration $i..."

  if [[ "$METHOD" == "freebayes" ]]; then
    pattern_mapper="${VC_DIR}/mapper*_freebayes_norm.vcf"
    bowtie_vcf="${VC_DIR}/bowtie2_freebayes_norm.vcf"
  else
    pattern_mapper="${VC_DIR}/mapper*_bcftools_norm.vcf"
    bowtie_vcf="${VC_DIR}/bowtie2_bcftools_norm.vcf"
  fi

  mapper_merged="$PLOT_DIR/mapper_merged.vcf"
  merge_vcfs "$mapper_merged" $pattern_mapper || true

  if [[ ! -s "$mapper_merged" ]]; then
    echo "Warning: no mapper VCFs to merge for iteration $i; skipping"
    continue
  fi
  if [[ ! -s "$bowtie_vcf" ]]; then
    echo "Warning: bowtie VCF missing for iteration $i; skipping"
    continue
  fi

  mapper_csv="$PLOT_DIR/vcfcompare_mapper.csv"
  bowtie_csv="$PLOT_DIR/vcfcompare_bowtie2.csv"

  python3 "$VCFCOMPARE_DIR/src/python/VCFCompare.py" "$TRUTH_VCF" "$mapper_merged" -o "$PLOT_DIR/vcfcompare_mapper" >/dev/null
  # output already at $mapper_csv

  python3 "$VCFCOMPARE_DIR/src/python/VCFCompare.py" "$TRUTH_VCF" "$bowtie_vcf" -o "$PLOT_DIR/vcfcompare_bowtie2" >/dev/null
  # output already at $bowtie_csv

  python3 analysis/plot_vcfcompare_iter.py \
    --iteration $i \
    --method "$METHOD" \
    --mapper-csv "$mapper_csv" \
    --bowtie-csv "$bowtie_csv" \
    --output-dir "$PLOT_DIR"

  # Compute per-gene confusion and plot horizontal bars
  python3 analysis/compute_variant_confusion_per_gene.py \
    --truth "$TRUTH_VCF" \
    --mapper "$mapper_merged" \
    --bowtie "$bowtie_vcf" \
    --out-tsv "$PLOT_DIR/variant_calling_confusion_matrix.tsv"

  python3 analysis/plot_variant_per_gene.py \
    --tsv "$PLOT_DIR/variant_calling_confusion_matrix.tsv" \
    --output "$PLOT_DIR/variant_per_gene_metrics.png" \
    --title "${METHOD^} Iteration $i - Per-Gene SNP Metrics"
done

python3 analysis/plot_vcfcompare_summary.py \
  --method "$METHOD" \
  --iterations "$ITERATIONS" \
  --root-dir "$PLOT_ROOT" \
  --output "$PLOT_ROOT/vcfcompare_summary.png"

echo "Done. Results in: $PLOT_ROOT"
