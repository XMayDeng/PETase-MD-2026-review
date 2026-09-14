"""Pure statistical extract; source SHA-256 c76cc56111b83fd3298f5cbdaeda499d13e23ba56543aeb5e8307e1763a9d583.
No trajectory, plotting, filesystem or historical main entry point.
"""
import hashlib
import math
import numpy as np
import pandas as pd

PET_FROM = 'PET_L10'
PET_TO = 'PET_L20'
N_BOOTSTRAP = 20000
BOOTSTRAP_SEED = 20260712
MIN_GROUP_N = 3
MIN_EFFECT = 0.1
MIN_DIRECTIONAL_PROBABILITY = 0.95
MIN_LOO_SIGN_AGREEMENT = 0.8
MIN_DIRECTION_SIGN_AGREEMENT = 2.0 / 3.0

def stable_seed(protein_id: str, residue: int) -> int:
    token = f'{BOOTSTRAP_SEED}|{protein_id}|{residue}'.encode('ascii')
    return int.from_bytes(hashlib.sha256(token).digest()[:8], 'little')

def bootstrap_means(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    indices = rng.integers(0, len(values), size=(N_BOOTSTRAP, len(values)))
    return values[indices].mean(axis=1)

def ci95(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, [0.025, 0.975])
    return (float(low), float(high))

def cliffs_delta(to_values: np.ndarray, from_values: np.ndarray) -> float:
    differences = to_values[:, None] - from_values[None, :]
    return float((np.sum(differences > 0) - np.sum(differences < 0)) / differences.size)

def leave_one_out_deltas(from_values: np.ndarray, to_values: np.ndarray) -> list[float]:
    deltas: list[float] = []
    if len(from_values) > 1:
        for index in range(len(from_values)):
            deltas.append(float(to_values.mean() - np.delete(from_values, index).mean()))
    if len(to_values) > 1:
        for index in range(len(to_values)):
            deltas.append(float(np.delete(to_values, index).mean() - from_values.mean()))
    return deltas

def sign_agreement(values: list[float], reference: float) -> float:
    if not values or reference == 0:
        return float('nan')
    nonzero = [value for value in values if value != 0]
    if not nonzero:
        return float('nan')
    target = math.copysign(1.0, reference)
    return float(np.mean([math.copysign(1.0, value) == target for value in nonzero]))

def direction_metrics(rows: pd.DataFrame, overall_delta: float) -> tuple[int, str, str, float, float, float, int, float]:
    from_dirs = set(rows[rows['pet_kind'] == PET_FROM]['direction'])
    to_dirs = set(rows[rows['pet_kind'] == PET_TO]['direction'])
    common = sorted(from_dirs & to_dirs)
    deltas: list[float] = []
    labels: list[str] = []
    for direction in common:
        from_mean = rows[(rows['pet_kind'] == PET_FROM) & (rows['direction'] == direction)]['contact_fraction'].mean()
        to_mean = rows[(rows['pet_kind'] == PET_TO) & (rows['direction'] == direction)]['contact_fraction'].mean()
        delta = float(to_mean - from_mean)
        deltas.append(delta)
        labels.append(f'{direction}:{delta:+.4f}')
    adjusted = float(np.mean(deltas)) if deltas else float('nan')
    agreement = sign_agreement(deltas, overall_delta)
    delta_range = float(max(deltas) - min(deltas)) if deltas else float('nan')
    target_sign = np.sign(overall_delta)
    effect_consistent = sum((abs(delta) >= MIN_EFFECT and np.sign(delta) == target_sign for delta in deltas))
    effect_consistent_fraction = effect_consistent / len(deltas) if deltas else float('nan')
    return (len(common), ';'.join(common), ';'.join(labels), adjusted, agreement, delta_range, effect_consistent, effect_consistent_fraction)

def classify_shift(row: dict[str, object]) -> tuple[str, str, bool]:
    delta = float(row['mean_delta_l20_minus_l10'])
    underpowered = bool(row['underpowered'])
    ci_excludes_zero = float(row['delta_ci95_low']) > 0 or float(row['delta_ci95_high']) < 0
    bootstrap_probability = float(row['bootstrap_direction_probability'])
    loo_agreement = float(row['loo_sign_agreement'])
    direction_agreement = float(row['direction_sign_agreement'])
    n_common_directions = int(row['n_common_directions'])
    n_direction_effect_consistent = int(row['n_direction_strata_effect_consistent'])
    required_direction_strata = max(2, math.ceil(MIN_DIRECTION_SIGN_AGREEMENT * n_common_directions))
    direction_adjusted_delta = float(row['direction_adjusted_delta'])
    contact_definition_sign_stable = bool(row['all_atom_heavy_delta_sign_same'])
    supported = bool(not underpowered and abs(delta) >= MIN_EFFECT and ci_excludes_zero and (bootstrap_probability >= MIN_DIRECTIONAL_PROBABILITY) and (loo_agreement >= MIN_LOO_SIGN_AGREEMENT) and (n_common_directions >= 2) and (direction_agreement >= MIN_DIRECTION_SIGN_AGREEMENT) and (abs(direction_adjusted_delta) >= MIN_EFFECT) and (n_direction_effect_consistent >= required_direction_strata) and contact_definition_sign_stable)
    if underpowered:
        shift_class = 'underpowered_descriptive'
    elif supported and delta > 0:
        shift_class = 'L20_enriched_supported_candidate'
    elif supported and delta < 0:
        shift_class = 'L10_enriched_supported_candidate'
    elif abs(delta) >= MIN_EFFECT and bootstrap_probability >= 0.8:
        shift_class = 'directional_but_uncertain_L20' if delta > 0 else 'directional_but_uncertain_L10'
    else:
        shift_class = 'no_resolved_shift'
    stable_core = bool(float(row['median_contact_l10']) >= 0.5 and float(row['median_contact_l20']) >= 0.5 and (float(row['support_ge_0_5_fraction_l10']) >= 0.5) and (float(row['support_ge_0_5_fraction_l20']) >= 0.5) and (abs(delta) <= 0.2))
    if supported:
        manuscript_status = 'candidate_with_caveats'
    elif stable_core and (not underpowered):
        manuscript_status = 'descriptive_core'
    else:
        manuscript_status = 'sensitivity_or_hypothesis_only'
    return (shift_class, manuscript_status, stable_core)
