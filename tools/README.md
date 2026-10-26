# Tools

Implementation modules used by ParaDISM commands.

- `liftover.py`: implementation for `python paradism.py liftover`.

Run liftover through the main ParaDISM entry point:

```bash
python paradism.py liftover --bed input.bed --ref ref.fa --positions positions.txt --output lifted.bed
```

The positions file is whitespace-delimited. The first field is the ParaDISM
contig name and the last field is `CHR:START-END:STRAND`, where coordinates
are 1-based inclusive and `STRAND` is `1` or `-1`. The interval must be the
exact span of the contig in `--ref`, including any flanking sequence.
