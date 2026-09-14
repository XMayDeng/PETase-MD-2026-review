"""STR-SENS-01 bootstrap and target aggregation from trajectory summaries.
Does not extract coordinates or infer absent measurements.
"""
import numpy as np
import pandas as pd
TARGET_ORDER = (('Pro00057', 185), ('Pro00075', 184), ('Pro00121', 164), ('Pro00121', 166))

TARGET_LABELS = {('Pro00057', 185): 'TfCut1 H185', ('Pro00075', 184): 'FoCut5a L184', ('Pro00121', 164): 'HiC T164', ('Pro00121', 166): 'HiC T166'}

PRIMARY_METRIC = {('Pro00057', 185): 'sidechain_conformation', ('Pro00075', 184): 'target_backbone_after_pocket_fit', ('Pro00121', 164): 'sidechain_conformation', ('Pro00121', 166): 'sidechain_conformation'}

BOOTSTRAP_SAMPLES = 20000

BOOTSTRAP_SEED = 20260726

def bootstrap_mean_ci(values: np.ndarray, seed: int) -> tuple[float, float, float]:
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    estimates = np.empty(BOOTSTRAP_SAMPLES, dtype=float)
    for index in range(BOOTSTRAP_SAMPLES):
        estimates[index] = rng.choice(values, size=len(values), replace=True).mean()
    low, high = np.quantile(estimates, [0.025, 0.975])
    return (float(values.mean()), float(low), float(high))

def target_summaries(trajectory: pd.DataFrame, spec: dict) -> pd.DataFrame:
    primary = trajectory[trajectory['window'] == 'primary'].copy()
    start = trajectory[trajectory['window'] == 'start'][['trajectory_id', 'target_resid', 'mean_primary_experimental_rmsd_angstrom']].rename(columns={'mean_primary_experimental_rmsd_angstrom': 'start_experimental_rmsd_angstrom'})
    primary = primary.merge(start, on=['trajectory_id', 'target_resid'], how='left')
    primary['analysis_minus_start_experimental_rmsd_angstrom'] = primary['mean_primary_experimental_rmsd_angstrom'] - primary['start_experimental_rmsd_angstrom']
    rows: list[dict[str, object]] = []
    sensitivity = [float(value) for value in spec['diagnostic_thresholds']['trajectory_overlap_fraction_sensitivity']]
    primary_overlap = float(spec['diagnostic_thresholds']['primary_trajectory_overlap_fraction'])
    seed_index = 0
    for protein_id, target in TARGET_ORDER:
        target_frame = primary[(primary['protein_id'] == protein_id) & (primary['target_resid'] == target)].copy()
        for subset_name, subset in (('all', target_frame), ('retained', target_frame[target_frame['retained'].astype(bool)]), ('retained_l10_l20', target_frame[target_frame['retained'].astype(bool) & target_frame['pet_kind'].isin(['PET_L10', 'PET_L20'])])):
            if subset.empty:
                raise ValueError('Missing required target/subset')
            mean_delta, delta_low, delta_high = bootstrap_mean_ci(subset['mean_primary_reference_closeness_delta_angstrom'].to_numpy(), BOOTSTRAP_SEED + seed_index)
            seed_index += 1
            mean_change, change_low, change_high = bootstrap_mean_ci(subset['analysis_minus_start_experimental_rmsd_angstrom'].to_numpy(), BOOTSTRAP_SEED + seed_index)
            seed_index += 1
            row: dict[str, object] = {'protein_id': protein_id, 'protein': subset.iloc[0]['protein'], 'target_resid': target, 'target_label': TARGET_LABELS[protein_id, target], 'primary_metric': PRIMARY_METRIC[protein_id, target], 'subset': subset_name, 'n_trajectories': len(subset), 'mean_primary_production_rmsd_angstrom': float(subset['mean_primary_production_rmsd_angstrom'].mean()), 'mean_primary_experimental_rmsd_angstrom': float(subset['mean_primary_experimental_rmsd_angstrom'].mean()), 'median_primary_experimental_overlap_fraction': float(subset['fraction_primary_experimental_overlap'].median()), 'median_primary_experimental_closer_fraction': float(subset['fraction_primary_experimental_closer'].median()), 'mean_reference_closeness_delta_angstrom': mean_delta, 'reference_closeness_delta_bootstrap_low_angstrom': delta_low, 'reference_closeness_delta_bootstrap_high_angstrom': delta_high, 'mean_analysis_minus_start_experimental_rmsd_angstrom': mean_change, 'analysis_minus_start_bootstrap_low_angstrom': change_low, 'analysis_minus_start_bootstrap_high_angstrom': change_high}
            for threshold in sensitivity:
                label = str(threshold).replace('.', 'p')
                qualifies = subset[subset['fraction_primary_experimental_overlap'] >= threshold]
                row[f'n_trajectories_exp_overlap_ge_{label}'] = len(qualifies)
                row[f'n_strata_exp_overlap_ge_{label}'] = (qualifies['pet_kind'].astype(str) + '/' + qualifies['direction'].astype(str)).nunique()
            rows.append(row)
    return pd.DataFrame(rows)
