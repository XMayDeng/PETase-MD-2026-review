# Data guide

This reference describes the data and executed statistical conventions used by
the reproduction commands. For installation and execution, start with the
[module README](README.md).

- [Shared conventions](#shared-conventions)
- [Trajectory index](#trajectory-index-126-rows)
- [Data S1](#data-s1-1532-rows-15-columns) and [Data S2](#data-s2-95-rows-67-columns)
- [Table S5](#table-s5-five-rows-26-columns)
- [Figure data and coordinates](#figure-data-and-structure-coordinates)
- [Statistical definitions](#statistical-definitions)
- [Table coverage](#table-and-supplementary-data-coverage)
- [Table S4 calculations](#table-s4-structural-sensitivity-calculations)

## Shared conventions

- Identifiers and categorical strings retain their source spelling. Residue
  numbers are simulation numbering. PDB reference identifiers do not convert
  simulation numbers to PDB-native numbers (notably FoCut5a versus 5AJH).
- L10/L20 denote PET chain lengths. Differences labeled L20-minus-L10 retain
  that sign. Fractions use the 0–1 scale unless explicitly stated otherwise.
- Whole trajectories are the independent analysis units. Repeated frames are
  not independent replicates. The registered primary window is 20–100 ns.
- Empty CSV cells remain empty; they are not converted into zeros. Source
  categorical flags and numeric precision are preserved as strings in S2/S5.

## Trajectory index: 126 rows

File: `metadata/trajectory_index.csv` (paths in this guide are repository-relative).

| Fields | Meaning |
|---|---|
| trajectory_id, run_id | Original trajectory and run identifiers, retained for provenance |
| enzyme, protein_id, pet_length | Panel identity and PET chain length |
| direction, replica_id | Initial-presentation and replicate labels |
| planned_production_ns | Production duration calculated from recorded MDP settings, in ns |
| final_logged_step, final_logged_time_ns | Last step/time recorded in production log; not a fresh trajectory-frame check |
| log_finished_mdrun, time_check_scope | Log completion flag and explicit scope of that check |
| retention_fraction | Hydrogen-inclusive Protein–PETL minimum distance ≤1.0 nm fraction over 20–100 ns |
| retention_class | retained ≥0.80; partial ≥0.50 and <0.80; detached <0.50 |
| analysis_window_ns | Registered analysis window, not the full production duration |
| production_gate_verdict, stability_verdict, pbc_verdict, clash_verdict | Preserved source QA categories; production FAIL must not be equated with a crashed run or used to delete detached trajectories |
| nvt_gen_seed | Recorded NVT velocity-generation seed |
| gromacs_version_registry, gromacs_version_log | Versions from the registry and the executed production log |
| case_manifest_provenance | Historical audit category: deployed run copy, or missing historical case configuration. The latter label is preserved, but is never used to select a replacement configuration in this package |

The historical `production_tpr_path` and `production_start_coordinate_path`
columns were omitted from this trajectory-summary index. Remaining field
values and rows are unchanged. Actual delivered production TPR paths are
indexed in `01_simulation/configurations/stages.csv`; GRO starting coordinates
can be extracted from those TPRs when required.

## Data S1: 1,532 rows, 15 columns

| Fields | Meaning |
|---|---|
| protein_id, enzyme, pdb_reference | Panel identifier, enzyme name, structural reference accession |
| simulation_residue_number, residue_name, residue_label | Residue identity in simulation numbering |
| structural_region, assigned_contact_class | Executed region and contact-class annotations, not a mechanistic classification |
| annotation_confidence, assignment_basis | Preserved source confidence category and annotation rationale |
| nearest_configured_landmark, nearest_configured_landmark_number | Source landmark identifier and residue number |
| nearest_configured_landmark_role, nearest_landmark_class | Source landmark role/class labels |
| nearest_landmark_ca_distance_A | Source Cα landmark distance in Å, rounded to six decimal places by the original projection |

The direct input is `full_panel_contact_landmark_annotation.csv`. The exporter
renames selected fields, adds enzyme/reference metadata and applies the original
ordering. It does not rerun a structural annotation analysis.

## Data S2: 95 rows, 67 columns

Direct inputs are `trajectory_level_contact_shift_candidates.csv`,
`manuscript_candidate_residues.csv` and the current Table S2 TeX membership check.
The full candidate population is retained: 29 tabulated entries and 66 others.
The 29 comprise eight supported L20 increases and 21 tabulated stable cores.

| Fields or family | Verified meaning / scope |
|---|---|
| protein_id, enzyme, pdb_reference, simulation_residue_number, residue_name, residue_label | Enzyme/reference/residue identity as above |
| structural_region, assigned_contact_class, annotation_confidence | Preserved source annotation labels |
| n_l10, n_l20 | Retained-trajectory counts for the compared groups |
| mean_contact_*, median_contact_*, sd_contact_* | Source contact-fraction summaries; the canonical residue-contact observable is a heavy-atom distance ≤0.45 nm |
| mean_delta_l20_minus_l10 | Source difference of L20 and L10 group means |
| mean_contact_*_ci95_*, delta_ci95_* | Source 95% interval endpoints, recalculated by the statistical pipeline |
| all_atom_*, heavy_*, in_all_atom_*, in_heavy_* | Preserved screen/rank and cross-selection fields; these do not relabel heavy-atom contact as a hydrogen bond |
| included_in_table_s3 | Historical column name: membership in the **current Table S2**, not current Table S3 |
| table_s3_classification | Historical column name: current Table S2 classification; empty for nonmembers |
| underpowered, stable_contact_core, shift_class, reporting_status, stability_evidence, stability_pass_* | Preserved source flags/labels; exact decision rules are given below |
| bootstrap_*, cliffs_delta_*, loo_*, n_cases_*, support_*, n_common_directions, common_directions, direction_*, n_direction_* | Source-traced statistical/stratification diagnostics; see the definitions below for estimators, weighting, thresholds and zero handling |

The CSV exporter preserves source field values and checks the exact 29-row
Table S2 values. The separate statistical pipeline recalculates the supported
bootstrap statistics; initial rankings and structural-QC flags remain inputs. The legacy `table_s3_rows` key in verifier output also
means current Table S2 membership. No legacy column was renamed.

## Table S5: five rows, 26 columns

The full-precision CSV and its TeX rendering both live in
`02_postprocessing/results/tables/`. Both are byte-identical to the approved sources.
The direct source is `interaction_shift_statistics.csv`.

- Identity fields: `analysis_target`, `protein_id`, `protein_name`, `metric_id`,
  `metric_label`, `scope`, `criterion_group` retain registered metric identifiers.
- `scope=per_pet_ring`: the fraction of PET aromatic rings satisfying the
  criterion is averaged over frames within each retained trajectory, then across
  trajectories. Zero-valued trajectories are included.
- `n_l10`, `n_l20`: TfCut1 counts 6/6 and HiC counts 5/8 for these selected rows.
- `mean_l10`, `mean_l20` and `delta_l20_minus_l10`: group means and their signed
  difference; fractions, not percentages. Mean/interval fields remain full precision.
- `mean_l10_ci95_low/high`, `mean_l20_ci95_low/high`, `delta_ci95_low/high`:
  source percentile-bootstrap 95% endpoints. The supplied Table S5 note documents
  20,000 whole-trajectory resamples within length groups for its displayed interval.
- `ci_excludes_zero`: source interval classification, checked against endpoints
  by the original generator. All five difference intervals include zero.
- `bootstrap_probability_delta_positive`, `bootstrap_direction_probability`,
  `loo_sign_agreement`, `direction_adjusted_delta`, `direction_sign_agreement`,
  `n_common_directions`, `evidence_class`: retained diagnostics/labels;
  exact statistical definitions and interval labels are documented below.

The displayed metrics are F210 aromatic-ring proximity and pi-like geometry,
and L138 nonpolar contact at 0.40 nm (primary), 0.45/0.50 nm (sensitivity).
Geometry definitions, post-selection and interval scope are retained in the
supplied Table S5 note. This projection is not evidence of reaction chemistry.

## Figure data and structure coordinates

`metadata/figures.json` maps each figure to its exact inputs and reference image.
Immediate source tables retain their additional columns; this guide describes
the fields used by the supported analyses, not every unused source field.

Distances with `_nm` use nanometers, and `_A`/`_a` distance fields use angstroms.
Published display conversions stay in the relevant drawing code. Angular and
periodic-rotamer calculations retain their source units and boundaries. The
frame-observable files contain derived measurements, not Cartesian trajectories.
Do not use their frame count as an independent-replicate count.

PDB files contain Cartesian coordinates in angstroms. Reference structures,
production-model protein structures, and selected protein/PET frame snapshots
remain separate input roles. Historical names such as `Figure_05_HiC.pdb` are
source identities, not the current figure number. `metadata/figures.json` gives
the exact current consumers, with one copy of each coordinate file.

The retention manifest omits only 11 original filesystem-location columns and
the audit timestamp. The reference-angle table omits only `source_path`.
The two representative-frame tables omit recorded local trajectory/input paths.
All these projections preserve every retained field string and every source row.
Source identities and transformations are in `metadata/sources.csv`; current
file hashes are in `metadata/file_manifest.csv`.

The existing name `candidate_case_hbond_values_zero_filled.csv` records the
upstream analysis convention of keeping zero-event retained trajectories. It
does not authorize filling newly missing records with zeros at run time.

## Statistical definitions

These are the executed conventions for contact-shift analysis (CON-03) and
local-interaction analysis (LOC-01). Their post-selection statistics are
descriptive/exploratory and unadjusted for multiplicity. A bootstrap probability
is neither a p-value nor a posterior probability. Whole-trajectory resampling
does not by itself establish strict independence between trajectories.

`02_postprocessing/code/validation/verify_statistics.py` checks the 95 contact
candidates and five Table S5 metrics from supplied measurements. Original
source hashes and extracted function names are recorded in `metadata/sources.csv`.

## Shared bootstrap and sign rules

Let x be the L10 trajectory values and y the L20 values. Delta is mean(y) minus
mean(x). Each group is separately resampled with replacement at its original
sample size, 20,000 times. The difference distribution subtracts the corresponding
resampled group means. Mean and difference intervals are NumPy quantiles 0.025
and 0.975 (the default linear quantile interpolation in these implementations).

- `bootstrap_probability_delta_positive` is the fraction of resampled
  differences strictly greater than zero.
- `bootstrap_direction_probability` equals that fraction when the observed
  delta is nonnegative, and one minus that fraction when delta is negative.
  Consequently, for a negative observed delta the implementation includes
  zero-valued bootstrap differences in this complement. Preserve this exact
  convention; do not relabel it as strictly P(delta < 0).
- Resampling uses `numpy.random.default_rng`. CON-03 seed is the first eight
  SHA-256 bytes of `20260712|protein_id|residue`, interpreted little-endian.
  LOC-01 uses `20260716|metric_id` in the same construction. The L10 draw precedes
  the L20 draw on the same seeded generator.
- `loo_sign_agreement` removes one trajectory at a time from either group
  (only from groups with more than one member), with the other group unchanged.
  Among nonzero leave-one-out differences, it reports the fraction with the
  sign of the original delta. It is missing when the reference delta is zero
  or no nonzero leave-one-out difference exists.

These definitions do not treat frames as resampling units and do not imply
paired L10/L20 trajectories.

## Data S2 statistics and screen annotations

| Field(s) | Exact meaning |
|---|---|
| n_l10, n_l20 | Number of trajectory contact fractions in each group |
| underpowered | True when min(n_l10,n_l20) < 3; a screen flag, not a formal power calculation |
| mean_contact_l10/l20, median_contact_l10/l20 | Arithmetic mean/median of the group trajectory fractions |
| sd_contact_l10/l20 | Sample standard deviation with ddof=1; missing for group size ≤1 |
| mean_contact_*_ci95_low/high, delta_ci95_low/high | Bootstrap interval endpoints using the shared rules above |
| mean_delta_l20_minus_l10 | Unrounded difference of group means |
| cliffs_delta_l20_vs_l10 | (number of y>x pairs minus number of y<x pairs)/(n_l20*n_l10); ties contribute zero, remain in the denominator |
| loo_min_delta, loo_max_delta | Minimum/maximum over all leave-one-out differences, including zero; missing if none |
| n_cases_ge_0_2_l10/l20, n_cases_ge_0_5_l10/l20 | Counts of trajectory fractions ≥0.2 or ≥0.5, respectively |
| support_ge_0_5_fraction_l10/l20 | Fraction of all trajectories in that group with contact fraction ≥0.5 |
| all_atom_mean_delta_20_100ns | L20-minus-L10 difference computed from the all-atom sensitivity contact definition over the same window |
| all_atom_heavy_delta_sign_same | Equality of NumPy signs of the all-atom and heavy-atom deltas, including zero-sign comparisons |

The primary values are heavy-atom contact fractions at 0.45 nm. The all-atom
columns retain a distinct sensitivity/screen role and do not define hydrogen bonds.

Candidate membership is the union of the all-atom and heavy-atom L10/L20 top
sets marked `is_top_n` by their supplied ranking tables, plus the three largest
and three smallest filtered L20-minus-L10 changes per enzyme in each definition.
The source's `*_pet_l10_rank`/`*_pet_l20_rank` fields are imported
`rank_within_pet` values only for selected top-set entries. Missing rank means
not present in that selected table, not zero contact. Corresponding `in_*_top8`
flags indicate a nonmissing imported rank; `in_*_top3_gain/loss` flags mark
membership in the respective three sorted extremes. “Gain/loss” here describes
the screen label, not proof of a supported positive/negative effect.

## Direction-stratum fields

Only initial-presentation directions present in both length groups are used.
Within each common direction, delta_direction is its L20 mean minus L10 mean.

| Field(s) | Exact meaning |
|---|---|
| n_common_directions | Number of shared directions |
| common_directions | Sorted shared labels joined with semicolons (S2) |
| direction_delta_values | Those labels and signed deltas printed to four decimals (S2); these display values do not replace full precision used in other fields |
| direction_adjusted_delta | Equal-weight arithmetic mean of the common-direction deltas, not weighted by trajectory counts |
| direction_sign_agreement | Fraction of nonzero common-direction deltas matching the overall delta sign; missing for zero reference or no nonzero deltas |
| direction_delta_range | Largest minus smallest common-direction delta, including zeros (S2) |
| n_direction_strata_effect_consistent | Number of common directions with absolute delta ≥0.10 and the same NumPy sign as the overall delta (S2) |
| direction_effect_consistent_fraction | Previous count divided by all common directions, including zero-delta strata (S2) |

When no common direction exists, mean/agreement/range/fraction fields are missing
and the common-direction and effect-consistent counts are zero. Thus the
denominators of direction_sign_agreement and direction_effect_consistent_fraction
are intentionally not identical when a stratum delta equals zero.

## Data S2 classification rules

A supported shift requires all of the following:

1. At least three trajectories per group.
2. Absolute overall delta ≥0.10 and a bootstrap difference interval strictly
   above or below zero.
3. Bootstrap direction probability ≥0.95 and leave-one-out sign agreement ≥0.80.
4. At least two common directions and direction sign agreement ≥2/3.
5. Absolute direction-adjusted delta ≥0.10.
6. At least max(2, ceil(2*n_common_directions/3)) common directions independently
   satisfying the same-sign absolute 0.10 effect rule.
7. Matching signs for the all-atom and heavy-atom differences.

`shift_class` applies this ordered decision:

- Underpowered groups: `underpowered_descriptive`.
- Supported positive/negative differences: `L20_enriched_supported_candidate`
  or `L10_enriched_supported_candidate`.
- Otherwise, absolute delta ≥0.10 and direction probability ≥0.80:
  `directional_but_uncertain_L20` or `_L10` by sign.
- Remaining rows: `no_resolved_shift`.

`stable_contact_core` requires median contact ≥0.5 in each length group,
at least half of trajectories reaching contact ≥0.5 in each group, and absolute
delta ≤0.20. This is an operational contact-core label, not temporal structural
stability or a catalytic classification. The flag alone does not impose n≥3.

`reporting_status` is the exported name of source `manuscript_status`. A supported
shift receives `candidate_with_caveats` first; otherwise a stable core with n≥3
per group receives `descriptive_core`; all remaining rows receive
`sensitivity_or_hypothesis_only`. A stable-core flag can coexist with another
reporting status. The existing summary and current Table S2 contain the first
two reporting categories. Legacy `included_in_table_s3` and
`table_s3_classification` preserve that Table S2 membership and display label.

`stability_pass_l10/l20` are counts, not Boolean values: the number of manifest
trajectories with stability_pass_bool true in each enzyme–length group.
`stability_evidence` is `stability-supported` only when every manifest trajectory
in both groups passes, otherwise `retained-contact_with_stability_caveat`.
This field aggregates the existing manifest flags; it does not recompute their
underlying structural-QC tests.

## Table S5-specific labels

The shared bootstrap, leave-one-out and common-direction definitions above apply
to trajectory-level ring-normalized geometry fractions instead of residue-contact
fractions. No S2 supported-candidate threshold rule is applied to these rows.

`ci_excludes_zero` is true exactly when delta_ci95_low > 0 or delta_ci95_high < 0.
`evidence_class` is `bootstrap_interval_above_zero`,
`bootstrap_interval_below_zero` or `bootstrap_interval_includes_zero`, using that
ordered interval test. An endpoint equal to zero belongs to the includes-zero
category. These are interval-position labels, not independent evidence grades.

## Table and supplementary data coverage

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

Table environments preserve reference wording, labels, and row order. They
require the manuscript preamble and are not standalone LaTeX documents.
Curated-cell CSVs retain the TeX spacing used to reproduce the reference text.

Table S4 is recalculated from trajectory summaries before formatting. Its
statistical check uses absolute tolerance 1e-12, zero relative tolerance,
identical identities/categories, and identical missing-value masks. The
formatted table must match exactly. To run this calculation separately:

```bash
python -B 02_postprocessing/code/statistics/recompute_structural_sensitivity.py --check
```

## Table S4: structural-sensitivity calculations

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
  The included bootstrap function preserves the original implementation.
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
