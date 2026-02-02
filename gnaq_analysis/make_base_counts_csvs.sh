#!/usr/bin/env bash
# Generate base-count CSVs at selected positions from BAMs.
# Usage: bash gnaq_analysis/make_base_counts_csvs.sh <bam_dir>
#
# Example:
#   bash gnaq_analysis/make_base_counts_csvs.sh gnaq_analysis/bwa-mem2_160_output/SRR5602384/final_outputs

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <bam_dir>"
    exit 1
fi

BAM_DIR="$1"

python - <<PY
import subprocess
from pathlib import Path
from collections import Counter

bam_dir = Path(r"${BAM_DIR}")
if not bam_dir.exists():
    raise SystemExit(f"Missing BAM dir: {bam_dir}")

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
        if c == '$':
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

for bam in sorted(bam_dir.rglob("*.bam")):
    if bam.name.endswith(".bai"):
        continue
    # Infer contig from filename
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
                     counts.get('A', 0), counts.get('C', 0), counts.get('G', 0), counts.get('T', 0), counts.get('N', 0)))

    out = bam.with_suffix('.base_counts.csv')
    with open(out, 'w') as f:
        f.write('contig,pos,ref,depth,A,C,G,T,N\n')
        for r in rows:
            f.write(','.join(map(str, r)) + '\n')
    print(out)
PY
