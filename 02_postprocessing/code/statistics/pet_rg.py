"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import itertools

import math

import numpy as np

import pandas as pd

BOOTSTRAP_SEED = 20_260_718

BOOTSTRAP_SAMPLES = 100_000

def exact_permutation_p(lcc: np.ndarray, tfcut1: np.ndarray) -> tuple[float, int]:
    values = np.concatenate([lcc, tfcut1])
    observed = float(lcc.mean() - tfcut1.mean())
    differences: list[float] = []
    for selected in itertools.combinations(range(len(values)), len(lcc)):
        selected_set = set(selected)
        group_a = values[list(selected)]
        group_b = values[[i for i in range(len(values)) if i not in selected_set]]
        differences.append(float(group_a.mean() - group_b.mean()))
    extreme = sum(abs(value) >= abs(observed) - 1e-12 for value in differences)
    return extreme / len(differences), len(differences)

def cliffs_delta(lcc: np.ndarray, tfcut1: np.ndarray) -> float:
    signs = [np.sign(x - y) for x in lcc for y in tfcut1]
    return float(np.mean(signs))

def hedges_g(lcc: np.ndarray, tfcut1: np.ndarray) -> float:
    n1, n2 = len(lcc), len(tfcut1)
    pooled_variance = (
        (n1 - 1) * np.var(lcc, ddof=1) + (n2 - 1) * np.var(tfcut1, ddof=1)
    ) / (n1 + n2 - 2)
    if pooled_variance <= 0:
        return math.nan
    cohen_d = (lcc.mean() - tfcut1.mean()) / math.sqrt(pooled_variance)
    correction = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    return float(correction * cohen_d)

def trajectory_bootstrap_ci(
    lcc: np.ndarray,
    tfcut1: np.ndarray,
    samples: int = BOOTSTRAP_SAMPLES,
    seed: int = BOOTSTRAP_SEED,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    lcc_samples = rng.choice(lcc, size=(samples, len(lcc)), replace=True).mean(axis=1)
    tf_samples = rng.choice(tfcut1, size=(samples, len(tfcut1)), replace=True).mean(axis=1)
    return tuple(float(value) for value in np.quantile(lcc_samples - tf_samples, [0.025, 0.975]))

def summarize_groups(trajectory: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    group_rows: list[dict[str, object]] = []
    for enzyme, group in trajectory.groupby("enzyme", sort=False):
        values = group["mean_rg_nm"].to_numpy(float)
        group_rows.append(
            {
                "enzyme": enzyme,
                "n_trajectories": len(values),
                "mean_of_trajectory_means_rg_nm": float(values.mean()),
                "sd_of_trajectory_means_rg_nm": float(values.std(ddof=1)),
                "median_trajectory_mean_rg_nm": float(np.median(values)),
                "min_trajectory_mean_rg_nm": float(values.min()),
                "max_trajectory_mean_rg_nm": float(values.max()),
            }
        )
    group_summary = pd.DataFrame(group_rows)
    lcc = trajectory.loc[trajectory["enzyme"] == "LCC", "mean_rg_nm"].to_numpy(float)
    tfcut1 = trajectory.loc[trajectory["enzyme"] == "TfCut1", "mean_rg_nm"].to_numpy(float)
    difference = float(lcc.mean() - tfcut1.mean())
    ci_low, ci_high = trajectory_bootstrap_ci(lcc, tfcut1)
    permutation_p, permutation_count = exact_permutation_p(lcc, tfcut1)
    comparison = pd.DataFrame(
        [
            {
                "contrast": "LCC_minus_TfCut1",
                "n_lcc_trajectories": len(lcc),
                "n_tfcut1_trajectories": len(tfcut1),
                "lcc_mean_of_trajectory_means_rg_nm": float(lcc.mean()),
                "tfcut1_mean_of_trajectory_means_rg_nm": float(tfcut1.mean()),
                "mean_difference_rg_nm": difference,
                "relative_difference_percent_of_tfcut1": float(100.0 * difference / tfcut1.mean()),
                "trajectory_bootstrap_ci95_low_nm": ci_low,
                "trajectory_bootstrap_ci95_high_nm": ci_high,
                "bootstrap_samples": BOOTSTRAP_SAMPLES,
                "bootstrap_seed": BOOTSTRAP_SEED,
                "exact_two_sided_permutation_p": permutation_p,
                "exact_permutation_assignments": permutation_count,
                "cliffs_delta": cliffs_delta(lcc, tfcut1),
                "probability_of_superiority": (cliffs_delta(lcc, tfcut1) + 1.0) / 2.0,
                "hedges_g": hedges_g(lcc, tfcut1),
                "complete_replica_mean_separation_lcc_above_tfcut1": bool(lcc.min() > tfcut1.max()),
            }
        ]
    )
    return group_summary, comparison
