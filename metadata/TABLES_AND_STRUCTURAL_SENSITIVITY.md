# Table assembly and the STR-SENS-01 statistical branch

The active manuscript endpoints are unchanged. All eight registered table/data
artifacts have executable checks, with different scientific starting points:

| Artifact | Implemented starting point | What remains upstream |
|---|---|---|
| Table 1 | Approved six-row enzyme-panel TeX cells and layout | Source-organism, thermal-class and citation claims are curated, not recalculated |
| Table S1 | Approved six-row residue-landmark TeX cells and layout | Landmark/numbering assignments are curated, not a new sequence/structure analysis |
| Table S2 | Existing 29 selected CON-03 candidate rows | Initial candidate selection and all-atom sensitivity calculation |
| Table S3 | Six full-precision STR-EQ-01 protein summaries and named-residue measurements | Coordinate fitting, sequence mapping and RMSD extraction |
| Table S4 | 252 target/trajectory/window summaries and 65 interaction/trajectory rows, plus four static reference measurements | Coordinate/RMSD extraction and framewise state/event assignment |
| Table S5 | Five LOC-01 rows, also checked from trajectory metric values | Extraction of the original geometry measurements |
| Data S1 | Executed residue-annotation table | Generation of the annotation map |
| Data S2 | Complete 95-candidate results, 29-row selection and Table S2 membership | Initial selection and full upstream measurement generation |

Table assembly preserves all approved wording, captions, notes, labels and row
order. Table 1 was originally extracted from V8. On 2026-09-14 the complete delivered table environment was verified unchanged in the active `main_V9.tex`.
Table S3 retains its current numbering; historical source filenames referring
to structural concordance as S2 are provenance, not additional tables.
The current TeX layout is a template with a single row placeholder, not a
second complete table used as a fallback. Curated-cell CSVs preserve TeX spacing
needed for byte identity and are explicitly document data.

Generated `.tex` files are table environments, not standalone LaTeX documents.
Typesetting them requires the manuscript's preamble, packages, custom units and
bibliography. The table entry checks their exact text and data without LaTeX.
Complete main/SI typesetting is handled separately from this data/code
repository; no standalone table PDF bundle is generated.

## Commands

From the package root, check the five newly connected table environments:

```bash
python -B 02_postprocessing/code/tables/build_tables.py --check
```

Generate them in a new directory, or select identifiers with `--tables`:

```bash
python -B 02_postprocessing/code/tables/build_tables.py --tables Table_01 Table_S01 Table_S02 Table_S03 Table_S04 --output-dir /tmp/petase_tables_rebuilt
```

Table S4 first recomputes its statistics from the bundled trajectory summaries.
It does not merely format the supplied aggregate CSV. The other three existing
Data S1/Data S2/Table S5 entry points remain under `02_postprocessing/code/tables/`.
The root `verify.py` runs all eight artifact checks and the statistical checks.
Existing differing outputs are never overwritten.

To recompute STR-SENS-01 independently:

```bash
python -B 02_postprocessing/code/statistics/recompute_structural_sensitivity.py --check
python -B 02_postprocessing/code/statistics/recompute_structural_sensitivity.py --output-dir /tmp/petase_structural_statistics
```

The output-directory form writes the actual recalculated numerical CSVs, not
copies of the supplied reference files. Floating-point serialization may differ
in final digits; validation uses absolute tolerance 1e-12, zero relative
tolerance, identical identities/categories and identical undefined-value masks.
The approved Table S4 TeX remains byte-exact after rounding.

## STR-SENS-01 fields and calculations

- `trajectory_target_summary.csv`: 63 unique trajectories across TfCut1,
  FoCut5a and HiC. Four targets give 84 trajectory–target pairs (the 21 HiC
  trajectories contribute both T164 and T166), each with start, early and primary
  windows, for 252 rows. The windows contain 1, 101 and 801 samples at 0 ns,
  0–10 ns and 20–100 ns, respectively, using the original 100-ps sampling.
- `*_rmsd_angstrom` fields are distances in Å. `mean_`, `median_`, `q05_` and
  `q95_` are the original within-window summaries, not independent frame
  replicates. This branch consumes the primary means, reference-closeness mean,
  and overlap/closer fractions; other source fields are preserved unchanged.
- `primary_reference_closeness_delta` is experimental-reference RMSD minus
  production-reference RMSD. The analysis-minus-start statistic subtracts each
  trajectory's 0-ns experimental-reference RMSD from its primary-window mean.
- `fraction_primary_experimental_overlap` uses the source diagnostic of 1.5 Å
  for side-chain conformation or 1.0 Å for FoCut5a L184 local backbone. The
  target aggregation counts trajectories at overlap fractions ≥0.01, ≥0.05 and
  ≥0.10; the Table S4 count uses ≥0.05. Subsets are all, retained and retained
  L10/L20. Each target has 21 trajectories in the all subset.
- Target intervals use 20,000 whole-trajectory mean resamples. The original
  seed is 20260726, incremented twice per target/subset in the original order.
  The bootstrap function is AST-identical to the historical implementation.
- `trajectory_interaction_state_coupling.csv`: 65 trajectory/target/indicator
  rows. Counts are sampled frames, event fractions use [0,1], and the signed
  state contrast is experimental-closer minus production-closer event fraction.
  All zero-event trajectories are retained. An absent state has an undefined
  conditional event fraction; a trajectory without both states has an undefined
  contrast. H185 has no eligible within-trajectory contrasts, so its group
  contrast and interval remain undefined, displayed as “Not estimable”.
- Interaction intervals resample finite trajectory contrasts only, with initial
  seed offset 500, incremented only for groups with estimable contrasts. Group
  order is preserved. Pooled frame counts are reconstructed from validated
  integer counts and exact count fractions, not by treating frames as replicates.
- Every trajectory identity, retention flag and stability label is checked
  against the existing 126-row retention manifest. Missing windows, duplicate
  identities, malformed fractions and missing required measurements fail.

The check covers all 12 target/subset rows, five interaction-summary rows and
four full-precision Table S4 CSV rows. The four static model/reference RMSDs are
checked inputs, not newly measured distances. Table S3 additionally reproduces
all 78 cells of its historical numerical projection before using the approved
current layout; its protein-summary coordinate-path fields are provenance
labels and are never opened by the table generator.

Exact inputs, hashes, output references and entry points are in
[tables.json](tables.json), [analyses.json](analyses.json) and
[provenance.csv](provenance.csv). These checks do not change scientific evidence
status or imply that full trajectory processing and MD setup are complete.
