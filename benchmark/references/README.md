# Benchmark References

This directory contains committed FASTA references for the reviewer-facing
demo and benchmark workflows. Each FASTA represents one multi-sequence
gene/paralog group and can be passed directly to `paradism.py` or
`benchmark/simulation/run_simulation.sh`.

Each FASTA should contain the gene of interest and the homologous paralog or
pseudogene sequences that reads may be assigned to. The FASTA record names are
used as ParaDISM local contig/gene names.

Additional reference annotation files:

- `pkd1_exons.bed`: PKD1/PKD1P gene-local exon intervals for the committed
  `pkd1_panel.fa` reference.

Included groups:

- `cfh_cfhr_cluster.fa`
- `cyp21_pair.fa`
- `cyp2d_cluster.fa`
- `fcgr_cluster.fa`
- `gba_pair.fa`
- `gnaq_pair.fa`
- `hba_pair.fa`
- `ncf1_triple.fa`
- `pkd1_panel.fa`
- `pms2_pair.fa`
- `sbds_pair.fa`
- `smn_pair.fa`
- `strc_pair.fa`

Example:

```bash
bash benchmark/simulation/run_simulation.sh \
  --group hba_pair \
  --reference benchmark/references/hba_pair.fa \
  --out results/simulation/hba_pair
```
