# Benchmarks

The public benchmark workflows are organized by input type:

- `references/`: committed FASTA references for each included gene/paralog
  group.
- `simulation/`: fully synthetic `dwgsim` read simulations. This is the
  fastest reproducible benchmark path after the demo.
- `giab/`: HG002/GIAB benchmark scripts using public but large GIAB read and
  truth inputs.

Run `demo/run_demo.sh` before running benchmarks. The demo is the shortest
end-to-end check that the environment and output paths are working.

A small simulation smoke test is documented in the top-level `README.md`.
