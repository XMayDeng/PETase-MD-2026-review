"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import math

import numpy as np

import pandas as pd

PROTEIN_ORDER = [
    "Pro00083",
    "Pro00062",
    "Pro00057",
    "Pro00137",
    "Pro00075",
    "Pro00121",
]

BOOTSTRAP_SEED = 20260721

N_BOOTSTRAP = 20_000

PRIMARY_LENGTHS = ("PET_L10", "PET_L20")

def percentile_interval(values: np.ndarray) -> tuple[float, float]:
    return tuple(float(x) for x in np.quantile(values, [0.025, 0.975]))

def length_balanced_bootstrap(
    l10: np.ndarray,
    l20: np.ndarray,
    seed: int,
    n_bootstrap: int = N_BOOTSTRAP,
) -> tuple[float, float, float, np.ndarray]:
    l10 = np.asarray(l10, dtype=float)
    l20 = np.asarray(l20, dtype=float)
    if l10.size == 0 or l20.size == 0:
        return math.nan, math.nan, math.nan, np.asarray([], dtype=float)
    point = float((l10.mean() + l20.mean()) / 2.0)
    rng = np.random.default_rng(seed)
    draw10 = l10[rng.integers(0, l10.size, size=(n_bootstrap, l10.size))].mean(axis=1)
    draw20 = l20[rng.integers(0, l20.size, size=(n_bootstrap, l20.size))].mean(axis=1)
    draws = (draw10 + draw20) / 2.0
    low, high = percentile_interval(draws)
    return point, low, high, draws

def summarize_profile(loop: pd.DataFrame) -> pd.DataFrame:
    rows = []
    primary = loop.loc[loop["pet_kind"].isin(PRIMARY_LENGTHS)]
    for (pid, relative), group in primary.groupby(["protein_id", "relative_position"], sort=False):
        l10 = group.loc[group["pet_kind"] == "PET_L10", "rmsf_nm"].to_numpy()
        l20 = group.loc[group["pet_kind"] == "PET_L20", "rmsf_nm"].to_numpy()
        point, low, high, _ = length_balanced_bootstrap(l10, l20, BOOTSTRAP_SEED + int(relative) + 100 * PROTEIN_ORDER.index(pid))
        first = group.iloc[0]
        rows.append(
            {
                "protein_id": pid,
                "protein": first["protein"],
                "mapping_class": first["mapping_class"],
                "relative_position": int(relative),
                "residue": int(first["residue"]),
                "residue_label": first["residue_label"],
                "is_marker": bool(first["is_marker"]),
                "n_l10": int(l10.size),
                "n_l20": int(l20.size),
                "mean_l10_rmsf_nm": float(l10.mean()),
                "mean_l20_rmsf_nm": float(l20.mean()),
                "length_balanced_mean_rmsf_nm": point,
                "bootstrap_ci_low_nm": low,
                "bootstrap_ci_high_nm": high,
            }
        )
    return pd.DataFrame(rows)

def summarize_metric_by_protein(table: pd.DataFrame, value_column: str, metric: str) -> pd.DataFrame:
    rows = []
    primary = table.loc[table["pet_kind"].isin(PRIMARY_LENGTHS)]
    for pid in PROTEIN_ORDER:
        group = primary.loc[primary["protein_id"] == pid]
        l10 = group.loc[group["pet_kind"] == "PET_L10", value_column].to_numpy()
        l20 = group.loc[group["pet_kind"] == "PET_L20", value_column].to_numpy()
        point, low, high, _ = length_balanced_bootstrap(l10, l20, BOOTSTRAP_SEED + 1000 + PROTEIN_ORDER.index(pid))
        rows.append(
            {
                "protein_id": pid,
                "protein": group.iloc[0]["protein"],
                "metric": metric,
                "n_l10": int(l10.size),
                "n_l20": int(l20.size),
                "mean_l10_nm": float(l10.mean()),
                "mean_l20_nm": float(l20.mean()),
                "length_balanced_mean_nm": point,
                "bootstrap_ci_low_nm": low,
                "bootstrap_ci_high_nm": high,
                "stability_verdict": group.iloc[0]["stability_verdict"],
            }
        )
    return pd.DataFrame(rows)

def summarize_pair_sites(pair: pd.DataFrame) -> pd.DataFrame:
    rows = []
    primary = pair.loc[pair["pet_kind"].isin(PRIMARY_LENGTHS)]
    for (pid, role), group in primary.groupby(["protein_id", "site_role"], sort=False):
        l10 = group.loc[group["pet_kind"] == "PET_L10", "rmsf_nm"].to_numpy()
        l20 = group.loc[group["pet_kind"] == "PET_L20", "rmsf_nm"].to_numpy()
        point, low, high, _ = length_balanced_bootstrap(l10, l20, BOOTSTRAP_SEED + 2000 + 10 * PROTEIN_ORDER.index(pid) + int(role[-1]))
        first = group.iloc[0]
        rows.append(
            {
                "protein_id": pid,
                "protein": first["protein"],
                "site_role": role,
                "residue": int(first["residue"]),
                "residue_label": first["residue_label"],
                "n_l10": int(l10.size),
                "n_l20": int(l20.size),
                "mean_l10_nm": float(l10.mean()),
                "mean_l20_nm": float(l20.mean()),
                "length_balanced_mean_nm": point,
                "bootstrap_ci_low_nm": low,
                "bootstrap_ci_high_nm": high,
                "stability_verdict": first["stability_verdict"],
            }
        )
    return pd.DataFrame(rows)
