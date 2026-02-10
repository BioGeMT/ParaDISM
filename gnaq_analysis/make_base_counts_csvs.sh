#!/usr/bin/env bash
# Generate base-count CSVs at selected positions from BAMs, collect into base_counts_<threshold>/ dir.
# Usage: bash gnaq_analysis/make_base_counts_csvs.sh <output_base_dir>
#
# Example:
#   bash gnaq_analysis/make_base_counts_csvs.sh gnaq_analysis/bwa-mem2_180_output

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <output_base_dir>"
    echo "  e.g. $0 gnaq_analysis/bwa-mem2_180_output"
    exit 1
fi

OUTPUT_BASE="$1"

if [[ ! -d "$OUTPUT_BASE" ]]; then
    echo "Error: directory $OUTPUT_BASE not found"
    exit 1
fi

# Extract threshold from dir name (e.g. bwa-mem2_180_output -> 180)
DIR_NAME="$(basename "$OUTPUT_BASE")"
THRESHOLD="$(echo "$DIR_NAME" | grep -oP '\d+(?=_output)')"
if [[ -z "$THRESHOLD" ]]; then
    exit 1
fi

BASE_COUNTS_DIR="$SCRIPT_DIR/base_counts_${THRESHOLD}"
mkdir -p "$BASE_COUNTS_DIR"

python - <<PY
import subprocess, shutil
from pathlib import Path
from collections import Counter

output_base = Path(r"${OUTPUT_BASE}")
base_counts_dir = Path(r"${BASE_COUNTS_DIR}")

positions = [
    ('GNAQ', 206083),
    ('GNAQ', 206100),
    ('GNAQP1', 765),
    ('GNAQP1', 786),
    ('GNAQP1', 826),
]

def parse_bases(bases, ref_base=None):
    counts = Counter()
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == '^':
            i += 2
            continue
        if c == '\$':
            i += 1
            continue
        if c in '+-':
            i += 1
            num = ''
            while i < len(bases) and bases[i].isdigit():
                num += bases[i]
                i += 1
            if num:
                i += int(num)
            continue
        if c in '.,':
            if ref_base:
                counts[ref_base.upper()] += 1
            else:
                counts['REF'] += 1
            i += 1
            continue
        if c in 'ACGTNacgtn':
            counts[c.upper()] += 1
            i += 1
            continue
        i += 1
    return counts

def mpileup_line(bam, contig, pos):
    cmd = ['samtools', 'mpileup', '-aa', '-r', f'{contig}:{pos}-{pos}', str(bam)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(res.stderr)
    line = res.stdout.strip()
    if not line:
        return None
    parts = line.split('\t')
    if len(parts) < 5:
        return None
    chrom, pos_s, ref_base, depth, bases = parts[:5]
    return chrom, int(pos_s), ref_base, int(depth), bases

# Find all sample dirs (directories containing final_outputs/)
sample_dirs = sorted(d for d in output_base.iterdir()
                     if d.is_dir() and (d / "final_outputs").is_dir())

if not sample_dirs:
    raise SystemExit(f"No sample directories with final_outputs/ found in {output_base}")

all_rows = []

for sample_dir in sample_dirs:
    final_outputs = sample_dir / "final_outputs"
    for bam in sorted(final_outputs.rglob("*.bam")):
        if bam.name.endswith(".bai"):
            continue
        if "_GNAQP1" in bam.name:
            target_contig = "GNAQP1"
        elif "_GNAQ" in bam.name:
            target_contig = "GNAQ"
        else:
            continue

        rows = []
        for contig, pos in positions:
            if contig != target_contig:
                continue
            line = mpileup_line(bam, contig, pos)
            if not line:
                rows.append((contig, pos, 'N', 0, 0, 0, 0, 0, 0))
                continue
            chrom, p, ref_base, depth, bases = line
            counts = parse_bases(bases, ref_base=ref_base)
            rows.append((chrom, p, ref_base.upper(), depth,
                         counts.get('A', 0), counts.get('C', 0),
                         counts.get('G', 0), counts.get('T', 0),
                         counts.get('N', 0)))

        # Write per-BAM CSV in the BAM directory
        per_bam_csv = bam.with_suffix('.base_counts.csv')
        with open(per_bam_csv, 'w') as f:
            f.write('contig,pos,ref,depth,A,C,G,T,N\n')
            for r in rows:
                f.write(','.join(map(str, r)) + '\n')
        print(f"  {per_bam_csv}")

        # Copy to base_counts/ dir
        dest = base_counts_dir / f"{bam.stem}_base_counts.csv"
        shutil.copy2(per_bam_csv, dest)

        # Collect for combined CSV
        for r in rows:
            all_rows.append((dest.name,) + r)

# Write combined CSV
combined = base_counts_dir / "base_counts_all.csv"
with open(combined, 'w') as f:
    f.write('file,contig,pos,ref,depth,A,C,G,T,N\n')
    for r in all_rows:
        f.write(','.join(map(str, r)) + '\n')
print(f"\nCombined: {combined}")
PY
