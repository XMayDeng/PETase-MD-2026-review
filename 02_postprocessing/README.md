# 02 — Delivered-data analysis

Run `python -B 02_postprocessing/run.py --output-dir ../analysis_results`.
This performs the declared numerical recalculations and table reconstruction
without modifying the supplied reference data. Each step has a log and a PASS
or FAIL result in the new directory. No GROMACS or raw XTC is required for this
entry. Python 3.12 and the recorded analysis libraries are required.

`inputs/` contains annotations, curated table cells and added per-trajectory
measurements. `results/` contains the original delivered frame observations,
trajectory summaries, aggregate references and publication tables. `code/`
contains calculation, validation and table-generation routines. Original source
hashes and the exact selected numerical-function provenance are in `metadata/`.

The initial measurement/extraction stage is not silently substituted by reading
final summary tables. See [REPRODUCIBILITY.md](../REPRODUCIBILITY.md) for the
starting level of each supported calculation and the steps requiring original
full trajectories. Supplied curated structural mappings and structural QA are
reference inputs, not newly inferred annotations.
