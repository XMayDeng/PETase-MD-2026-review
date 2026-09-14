# Source-traced statistical definitions for Data S2 and Table S5

Verified 2026-09-09. These definitions describe the executed analysis conventions,
not newly recommended estimators. The original data, labels and source scripts
were not changed. Statistical outputs are post-selection, descriptive/exploratory
and unadjusted for multiplicity. A bootstrap probability is not a p-value or a
posterior probability. Strict independence is not asserted merely because whole
trajectories, rather than frames, are the resampling units.

The separately ported STR-SENS-01 trajectory-summary statistics, exact bootstrap
ordering and undefined-state rules are documented in
[TABLES_AND_STRUCTURAL_SENSITIVITY.md](TABLES_AND_STRUCTURAL_SENSITIVITY.md).

## Source identities and verification scope

| Analysis | Source code snapshot | SHA-256 |
|---|---|---|
| CON-03 | analyze_heavy_atom_contact_shift_uncertainty.py | c76cc56111b83fd3298f5cbdaeda499d13e23ba56543aeb5e8307e1763a9d583 |
| LOC-01 | run_nonpolar_candidate_interactions.executed.py | b34e12326979f6312b0a5e013043479e06aff80acdf04cdf8dbe6bc361b86c53 |

Source-code inspection covered the relevant statistics, classification and
screen-annotation functions. All 95 CON-03 classification triplets matched the
existing source table. For the five S5 metrics, 63 existing trajectory-metric
records reproduced all numerical statistics to absolute tolerance 1e-12
(largest observed difference about 1.1e-16); categorical fields matched exactly.
This checks derived-table-to-statistic consistency, not raw-trajectory extraction.
The complete original analysis programs are not included. This package provides
AST-preserving extracts of their inspected pure statistical functions and two
existing trajectory-summary tables. The entry point
`python -B 02_postprocessing/code/validation/verify_statistics.py` verifies CON-03 heavy-contact statistics
for 95 candidates and these five S5 metrics. Raw-trajectory extraction, initial
candidate selection/ranking, all-atom sensitivity calculation and structural-QC
generation remain outside this bundled verification.

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
