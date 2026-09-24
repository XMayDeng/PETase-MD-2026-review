# Reproduction scope and verification

## Supported starting points

| Task | Included starting material | Entry point |
|---|---|---|
| Replay production MD | 126 original protein–PET production TPRs | [01_simulation/run.py](01_simulation/README.md) |
| Recalculate statistics | Trajectory-level measurements and specified frame observables | [02_postprocessing/run.py](02_postprocessing/README.md) |
| Generate six tables and Data S1/S2 | Measurements, curated cells, and table templates | [02_postprocessing/run.py](02_postprocessing/README.md) |
| Redraw 15 numbered figures | Numerical data, selected coordinates, and drawing code | [03_figures/run.py](03_figures/README.md) |
| Redraw the TOC graphic | Python drawing sources and a baseline image check | [TOC command](03_figures/README.md#toc-graphic) |

The numerical pipeline recalculates retention summaries, pooled contact
profiles, contact-class composition and uncertainty, the 95 preselected contact
candidates, PET radius of gyration, W-loop RMSF, PHL7 rotamer/co-engagement
statistics and mixture-model sensitivity, hydrogen-bond and angle-sensitivity
statistics, local-interaction statistics, and Table S4 aggregates.

Table 1 and Table S1 use curated cells. Initial candidate screening, all-atom
sensitivity, structural correspondences, coordinate fitting, and structural QC
remain supplied inputs where identified in the
[data guide](02_postprocessing/DATA_GUIDE.md#table-and-supplementary-data-coverage).
Whole trajectories, not individual frames, are the statistical units.

## What requires the original trajectories?

The original full XTC files are needed to extract new atom-level measurements,
change time windows or selections, repeat initial screening, or recluster the
original ensembles. Frame-observable CSVs contain derived measurements, not
Cartesian coordinates. Figure snapshots illustrate selected configurations and
do not replace an ensemble.

This repository includes all 126 production TPRs, 378 effective NVT/NPT/MD MDPs,
recorded preparation MDPs, source/reference structures, frozen PET conformers,
readable topologies, force-field includes, and system/preparation metadata.
It omits full original trajectories, the full-panel cluster archive, earlier-stage
TPRs, full solvated-coordinate archives, original final checkpoints, and the
complete historical docking outputs. Protocol source is included for inspection;
the supported simulation command starts at an existing production TPR, not at
preparation, docking, or equilibration.

TPR replay starts a new execution of the recorded production stage. It does not
continue an original final checkpoint or recover the original trajectory exactly.
The delivered chemical states and simulation parameters are preserved. Execution
and numerical agreement checks do not establish sampling convergence or validate
the scientific suitability of a model.

## Verification commands

Run from the repository root using the [tested environment](environment/README.md):

```bash
# File integrity and declared input dependencies
python -B verify.py --checksums-only

# Preparation metadata, statistical calculations, tables, and Data S1/S2
python -B verify.py

# Also parse all production TPRs and regenerate the numbered figures
python -B verify.py --gmx /path/to/gmx --render --chimerax /path/to/ChimeraX --output-dir ../verified_figures
```

Use a new external output directory. Numerical checks use the declared
tolerances, identities, and missing-value rules. Table text and supplementary
CSV exports are compared exactly with the references. Figure checks require
exact PNG agreement in the recorded Python/font/ChimeraX environment; this is
not a promise of pixel identity with arbitrary software or graphics drivers.
Use the separate TOC command in 03 to validate that graphic as well. A bounded
MD execution check is documented in 01; parsing a TPR alone is not MD execution.

## Machine-readable records

- `metadata/file_manifest.csv`: delivered file paths, sizes, and checksums.
- `metadata/trajectory_index.csv`: trajectory identities and retained populations.
- `metadata/analyses.json`, `tables.json`, and `figures.json`: inputs, references,
  and entry points used by the checking and generation scripts.
- `metadata/sources.csv`: original source identities and documented adaptations,
  including extracted statistical functions. Source names identify historical
  inputs, not paths that must exist on the reader's machine.
- Simulation-specific indexes remain alongside their inputs in 01.

The verifier ignores Git internals, Python caches, `.venv/`, `work/`, and
`validation_runs/`; other unlisted files fail the inventory check. Routine outputs
should remain outside the checkout. Manuscript typesetting is outside this
data/code repository's scope.
