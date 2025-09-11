#!/usr/bin/env python3
import argparse
import os

def parse_vcf_snvs(vcf_path):
    per_gene = {}
    with open(vcf_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            chrom, pos, _, ref, alt = parts[:5]
            # Only simple SNVs (single ALT, 1bp ref/alt)
            if ',' in alt:
                continue
            if len(ref) == 1 and len(alt) == 1:
                pos_i = int(pos)
                per_gene.setdefault(chrom, set()).add((pos_i, ref, alt))
    return per_gene

def main():
    ap = argparse.ArgumentParser(description='Compute per-gene variant confusion (SNV) for mapper and bowtie2')
    ap.add_argument('--truth', required=True)
    ap.add_argument('--mapper', required=True)
    ap.add_argument('--bowtie', required=True)
    ap.add_argument('--out-tsv', required=True)
    args = ap.parse_args()

    truth = parse_vcf_snvs(args.truth)
    mapper = parse_vcf_snvs(args.mapper)
    bowtie = parse_vcf_snvs(args.bowtie)

    genes = sorted(set(truth.keys()) | set(mapper.keys()) | set(bowtie.keys()))

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    with open(args.out_tsv, 'w') as out:
        out.write('\t'.join(['Gene','Mapper_TP','Mapper_FP','Mapper_FN','Bowtie2_TP','Bowtie2_FP','Bowtie2_FN'])+'\n')
        for g in genes:
            t = truth.get(g, set())
            m = mapper.get(g, set())
            b = bowtie.get(g, set())
            mtp = len(t & m); mfp = len(m - t); mfn = len(t - m)
            btp = len(t & b); bfp = len(b - t); bfn = len(t - b)
            out.write(f"{g}\t{mtp}\t{mfp}\t{mfn}\t{btp}\t{bfp}\t{bfn}\n")

if __name__ == '__main__':
    main()

