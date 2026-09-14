# Data dictionary and preservation notes

Status: working-copy dictionary. Field meanings are documented below; the
CON-03/LOC-01 statistical diagnostics have now been traced to their source code.
See [STATISTICAL_DEFINITIONS.md](STATISTICAL_DEFINITIONS.md) for exact rules,
zero handling, thresholds and verification scope. This does not certify the
entire upstream workflow or establish publication rights.

The additional Table 1/S1–S4 inputs, 63-trajectory STR-SENS-01 aggregation,
undefined-state rules and statistical units are documented in
[TABLES_AND_STRUCTURAL_SENSITIVITY.md](TABLES_AND_STRUCTURAL_SENSITIVITY.md).

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

## trajectory_index.csv: 126 rows

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
| mean_contact_*_ci95_*, delta_ci95_* | Source 95% interval endpoints, not recomputed here |
| all_atom_*, heavy_*, in_all_atom_*, in_heavy_* | Preserved screen/rank and cross-selection fields; these do not relabel heavy-atom contact as a hydrogen bond |
| included_in_table_s3 | Historical column name: membership in the **current Table S2**, not current Table S3 |
| table_s3_classification | Historical column name: current Table S2 classification; empty for nonmembers |
| underpowered, stable_contact_core, shift_class, reporting_status, stability_evidence, stability_pass_* | Preserved source flags/labels; exact decision rules are given in the linked statistical definitions |
| bootstrap_*, cliffs_delta_*, loo_*, n_cases_*, support_*, n_common_directions, common_directions, direction_*, n_direction_* | Source-traced statistical/stratification diagnostics; see the linked definitions for estimators, weighting, thresholds and zero handling |

The exporter preserves all retained source field values and checks the exact
29-row Table S2 values. It does not regenerate the bootstrap, rankings or
stability diagnostics. The legacy `table_s3_rows` key in verifier output also
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
  exact statistical definitions and interval labels are documented in the linked file.

The displayed metrics are F210 aromatic-ring proximity and pi-like geometry,
and L138 nonpolar contact at 0.40 nm (primary), 0.45/0.50 nm (sensitivity).
Geometry definitions, post-selection and interval scope are retained in the
supplied Table S5 note. This projection is not evidence of reaction chemistry.

## MDP and source tables

MDP names and values are GROMACS resolved settings, not newly designed templates.
Only full comment lines with private-root path markers were removed in copies.
All non-comment lines were compared against the original and are identical.
No MD production command is supplied because the complete simulation-input
dependency set is not included. The small structure subset is for figure rendering,
not for restarting the 126 production systems. Immediate source tables are
preserved, including unused columns;
their additional fields are not newly interpreted by this dictionary.

## Figure data and structure coordinates

`dataset_inventory.csv` enumerates every packaged CSV/CSV.gz, its row count,
exact field names, hash and figure consumers. The file inventory is a schema
inventory, not a claim that every additional source field has been audited.

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
Their exact source/copy hashes and transformations are in `provenance.csv`.

The existing name `candidate_case_hbond_values_zero_filled.csv` records the
upstream analysis convention of keeping zero-event retained trajectories. It
does not authorize filling newly missing records with zeros at run time.
