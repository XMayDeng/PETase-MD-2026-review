# Reproduction scope

This repository is a self-contained review subset of the manuscript's
simulation and analysis materials. It combines original production inputs,
the numerical data used by the supported analyses, and figure sources.
Data for the commands described here are local to this
repository; no Hugging Face account or upstream working directory is needed.

## Supported starting points

| Task | Starting material | Entry |
|---|---|---|
| Check files and dependencies | File manifest, figure/analysis registries | `python -B verify.py --checksums-only` |
| Inspect or replay production MD | 126 original production TPRs | `01_simulation/run.py` |
| Recalculate reported statistics | Delivered trajectory-level measurements and specified frame observations | `02_postprocessing/run.py` |
| Rebuild tables and supplementary CSVs | Measurements, curated cells, and table templates | `02_postprocessing/run.py` |
| Redraw the 15 numbered figures | Delivered data, selected PDB coordinates, and rendering code | `03_figures/run.py` |
| Redraw the separate TOC | Self-contained Python drawing sources and baseline check | `03_figures/toc/make_figure.py` |

See the corresponding folder README for full commands and dependencies.
[Validation results](metadata/TEST_RESULTS.md) distinguish file checks,
numerical recalculation, rendering and bounded MD execution.
Manuscript and SI typesetting are outside this repository's scope.

## Numerical coverage

The analysis entry recalculates retention summaries, pooled top-eight and
top-20 contact profiles, length-stratified contact composition and bootstrap
uncertainty, statistics for the 95 preselected contact candidates, PET radius
of gyration comparisons, W-loop RMSF summaries, PHL7 rotamer/co-engagement
summaries and mixture-model sensitivity, hydrogen-bond and angle-sensitivity
statistics, local-interaction statistics, and the Table S4 aggregate.

Calculations start from delivered measurements, not from reconstructed
coordinates. Initial candidate screening, all-atom sensitivity, manually
curated residue correspondences, experimental-reference fitting and structural
QC remain supplied evidence inputs where identified in
[the table/analysis coverage map](metadata/TABLES_AND_STRUCTURAL_SENSITIVITY.md).
Table 1 and Table S1 preserve curated cells. Each stage's inputs, expected
outputs, population and comparison rules are recorded in `metadata/analyses.json`
and `metadata/tables.json`.

Whole trajectories are the independent units. Primary contact/local-geometry
measurements use the registered 20–100 ns window, generally at 10 ps spacing;
coordinate-sensitivity measurements use their recorded 100 ps spacing.
Cutoffs, atom selections, retained denominators, zero-event handling, rounding
and tolerances are preserved from the source calculations. The full-precision
data, source hashes and original numerical-function provenance remain available
in the data files and `metadata/statistical_core_provenance.csv`.

## Included and omitted simulation material

Included: all 126 protein–PET production TPRs; all 378 effective NVT/NPT/MD MDPs;
recorded preparation MDPs; source/reference structures; frozen PET conformers;
readable protein/PET topologies and force-field includes; system/seed/preparation
metadata; and the PDB snapshots actually needed to redraw the structural figures.

Not included: full original XTC trajectories; the separate full-panel cluster
archive; EM, NVT, NPT and PET-only stage TPRs; full solvated-coordinate archives;
original final checkpoints; or the complete historical docking output archive.
These exclusions reduce download size without substituting cluster medoids
for the numerical observations used in the paper.

A production TPR can start a new execution of its recorded stage and can be
used to inspect/extract its input structure. It cannot recover the original
100 ns trajectory. Recomputing original atom-level distances, changing time
windows or cutoffs, repeating the original screening, or reclustering original
trajectories requires the matching original trajectories and selections.
The displayed structure snapshots are illustrative inputs, not a substitute
for an ensemble or evidence of convergence.

The original simulation histories, source models and chemical preparation
states are preserved rather than altered during packaging. Fifteen case YAMLs
are historical deployed copies; the three TfCut1 cases retain their explicit
fallback provenance. The presence or readability of an input does not validate
its scientific appropriateness.

## File integrity and provenance

`metadata/file_manifest.csv` records the exact included paths, sizes and
SHA-256 checksums. It is the upload/export allowlist together with the manifest
itself. `metadata/review_source_provenance.csv` records each review file's
source identity or documented adaptation. No source experiment or HF-package
file is edited by this extraction.

The verifier excludes Git internals, Python caches, a local virtual environment,
`work/` and `validation_runs/`. Other unlisted files cause a file-set failure.

## Review access and remaining author decisions

A private GitHub link is usable only by authorized readers. Arrange reviewer
access before relying on it in a submission; visibility is not changed by the
reproduction commands. Fixed-version citation/archiving and any DOI should
identify the actual delivered revision, not an unversioned promise.

Manuscript submission and author approval are handled separately. Reviewer
access, fixed-version archiving and third-party redistribution permissions
remain the authors' responsibility.
Reproduction tests do not certify sampling convergence, scientific validity
or journal acceptance. This package provides the explicitly listed review
capabilities, not an assertion that every possible reanalysis is supported.
