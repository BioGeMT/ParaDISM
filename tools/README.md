# Tools

Small command-line utilities that support ParaDISM outputs.

- `liftover.py`: converts ParaDISM VCF or BED files from gene-local reference
  coordinates back to chromosomal coordinates using a user-provided positions
  file.

The same liftover functionality is exposed through:

```bash
python paradism.py liftover --bed input.bed --positions positions.txt --output lifted.bed
```

The positions file is whitespace-delimited. The first field is the ParaDISM
contig name and the last field is `CHR:START-END:STRAND`, where coordinates
are 1-based inclusive and `STRAND` is `1` or `-1`.
