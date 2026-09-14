"""Pure statistical extract; source SHA-256 b34e12326979f6312b0a5e013043479e06aff80acdf04cdf8dbe6bc361b86c53.
No trajectory, plotting, filesystem or historical main entry point.
"""
import hashlib
import math
import numpy as np
import pandas as pd

N_BOOTSTRAP = 20000
BOOTSTRAP_SEED = 20260716
EXPECTED_GROUP_COUNTS = {('Pro00057', 'PET_L10'): 6, ('Pro00057', 'PET_L20'): 6, ('Pro00121', 'PET_L10'): 5, ('Pro00121', 'PET_L20'): 8}

def stable_seed(*tokens: object) -> int:
    payload = '|'.join((str(token) for token in (BOOTSTRAP_SEED, *tokens))).encode('ascii')
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], 'little')

def bootstrap_means(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    indices = rng.integers(0, len(values), size=(N_BOOTSTRAP, len(values)))
    return values[indices].mean(axis=1)

def ci95(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, [0.025, 0.975])
    return (float(low), float(high))

def leave_one_out_deltas(l10: np.ndarray, l20: np.ndarray) -> list[float]:
    values: list[float] = []
    for index in range(len(l10)):
        if len(l10) > 1:
            values.append(float(l20.mean() - np.delete(l10, index).mean()))
    for index in range(len(l20)):
        if len(l20) > 1:
            values.append(float(np.delete(l20, index).mean() - l10.mean()))
    return values

def sign_agreement(values: list[float], reference: float) -> float:
    nonzero = [value for value in values if value != 0]
    if reference == 0 or not nonzero:
        return float('nan')
    return float(np.mean([np.sign(value) == np.sign(reference) for value in nonzero]))

def direction_statistics(group: pd.DataFrame, overall_delta: float) -> tuple[float, float, int]:
    common = sorted(set(group.loc[group['pet_kind'].eq('PET_L10'), 'direction']) & set(group.loc[group['pet_kind'].eq('PET_L20'), 'direction']))
    deltas: list[float] = []
    for direction in common:
        l10 = group.loc[group['pet_kind'].eq('PET_L10') & group['direction'].eq(direction), 'value']
        l20 = group.loc[group['pet_kind'].eq('PET_L20') & group['direction'].eq(direction), 'value']
        deltas.append(float(l20.mean() - l10.mean()))
    return (float(np.mean(deltas)) if deltas else float('nan'), sign_agreement(deltas, overall_delta), len(deltas))

def build_shift_statistics(metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = ['analysis_target', 'protein_id', 'protein_name', 'metric_id', 'metric_label', 'scope', 'criterion_group']
    for keys, group in metrics.groupby(group_columns, sort=False):
        base = dict(zip(group_columns, keys))
        l10 = group.loc[group['pet_kind'].eq('PET_L10'), 'value'].to_numpy(float)
        l20 = group.loc[group['pet_kind'].eq('PET_L20'), 'value'].to_numpy(float)
        expected = (EXPECTED_GROUP_COUNTS[str(base['protein_id']), 'PET_L10'], EXPECTED_GROUP_COUNTS[str(base['protein_id']), 'PET_L20'])
        if (len(l10), len(l20)) != expected:
            raise ValueError(f'Unexpected trajectory counts for {base['metric_id']}')
        rng = np.random.default_rng(stable_seed(base['metric_id']))
        bootstrap_l10 = bootstrap_means(l10, rng)
        bootstrap_l20 = bootstrap_means(l20, rng)
        bootstrap_delta = bootstrap_l20 - bootstrap_l10
        l10_low, l10_high = ci95(bootstrap_l10)
        l20_low, l20_high = ci95(bootstrap_l20)
        delta_low, delta_high = ci95(bootstrap_delta)
        delta = float(l20.mean() - l10.mean())
        probability_positive = float(np.mean(bootstrap_delta > 0))
        direction_probability = probability_positive if delta >= 0 else 1.0 - probability_positive
        loo_agreement = sign_agreement(leave_one_out_deltas(l10, l20), delta)
        direction_adjusted, direction_agreement, n_directions = direction_statistics(group, delta)
        if delta_low > 0:
            evidence = 'bootstrap_interval_above_zero'
        elif delta_high < 0:
            evidence = 'bootstrap_interval_below_zero'
        else:
            evidence = 'bootstrap_interval_includes_zero'
        rows.append({**base, 'n_l10': len(l10), 'n_l20': len(l20), 'mean_l10': float(l10.mean()), 'mean_l10_ci95_low': l10_low, 'mean_l10_ci95_high': l10_high, 'mean_l20': float(l20.mean()), 'mean_l20_ci95_low': l20_low, 'mean_l20_ci95_high': l20_high, 'delta_l20_minus_l10': delta, 'delta_ci95_low': delta_low, 'delta_ci95_high': delta_high, 'bootstrap_probability_delta_positive': probability_positive, 'bootstrap_direction_probability': direction_probability, 'loo_sign_agreement': loo_agreement, 'direction_adjusted_delta': direction_adjusted, 'direction_sign_agreement': direction_agreement, 'n_common_directions': n_directions, 'ci_excludes_zero': bool(delta_low > 0 or delta_high < 0), 'evidence_class': evidence})
    return pd.DataFrame(rows).sort_values(['analysis_target', 'metric_id'])
