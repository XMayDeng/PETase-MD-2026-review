"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import hashlib

import numpy as np

import pandas as pd

N_BOOTSTRAP = 20_000

BOOTSTRAP_SEED = 20260715

EXPECTED_GROUP_COUNTS = {
    ("Pro00057", "PET_L10"): 6,
    ("Pro00057", "PET_L20"): 6,
    ("Pro00121", "PET_L10"): 5,
    ("Pro00121", "PET_L20"): 8,
}

def stable_seed(*tokens: object) -> int:
    text = "|".join(str(token) for token in (BOOTSTRAP_SEED, *tokens)).encode("ascii")
    return int.from_bytes(hashlib.sha256(text).digest()[:8], "little")

def bootstrap_means(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    indices = rng.integers(0, len(values), size=(N_BOOTSTRAP, len(values)))
    return values[indices].mean(axis=1)

def ci95(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, [0.025, 0.975])
    return float(low), float(high)

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
        return float("nan")
    return float(np.mean([np.sign(value) == np.sign(reference) for value in nonzero]))

def direction_statistics(group: pd.DataFrame, metric: str, overall_delta: float) -> tuple[float, float, int]:
    common = sorted(
        set(group.loc[group["pet_kind"].eq("PET_L10"), "direction"])
        & set(group.loc[group["pet_kind"].eq("PET_L20"), "direction"])
    )
    deltas: list[float] = []
    for direction in common:
        l10 = group.loc[
            group["pet_kind"].eq("PET_L10") & group["direction"].eq(direction), metric
        ]
        l20 = group.loc[
            group["pet_kind"].eq("PET_L20") & group["direction"].eq(direction), metric
        ]
        deltas.append(float(l20.mean() - l10.mean()))
    return (
        float(np.mean(deltas)) if deltas else float("nan"),
        sign_agreement(deltas, overall_delta),
        len(deltas),
    )

def evidence_classification(
    delta: float,
    low: float,
    high: float,
    direction_probability: float,
    loo_agreement: float,
    direction_agreement: float,
) -> str:
    ci_excludes_zero = low > 0 or high < 0
    robust = bool(
        ci_excludes_zero
        and direction_probability >= 0.95
        and loo_agreement >= 0.80
        and direction_agreement >= 2.0 / 3.0
    )
    if robust and delta > 0:
        return "resolved_L20_increase"
    if robust and delta < 0:
        return "resolved_L20_decrease"
    if direction_probability >= 0.80:
        return "directional_but_unresolved"
    return "no_resolved_hbond_shift"

def build_shift_statistics(case_values: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    keys = [
        "criterion",
        "angle_definition",
        "protein_id",
        "protein_name",
        "residue_number",
        "residue_label",
    ]
    for group_keys, group in case_values.groupby(keys, sort=False):
        base = dict(zip(keys, group_keys))
        for metric, metric_label in (
            ("any_hbond_occupancy", "any_residue_atom"),
            ("sidechain_hbond_occupancy", "sidechain_only"),
        ):
            l10 = group.loc[group["pet_kind"].eq("PET_L10"), metric].to_numpy(float)
            l20 = group.loc[group["pet_kind"].eq("PET_L20"), metric].to_numpy(float)
            if (len(l10), len(l20)) != (
                EXPECTED_GROUP_COUNTS[(str(base["protein_id"]), "PET_L10")],
                EXPECTED_GROUP_COUNTS[(str(base["protein_id"]), "PET_L20")],
            ):
                raise ValueError(f"Unexpected L10/L20 counts for {base}")
            rng = np.random.default_rng(
                stable_seed(
                    base["criterion"],
                    base["protein_id"],
                    base["residue_number"],
                    metric_label,
                )
            )
            bootstrap_delta = bootstrap_means(l20, rng) - bootstrap_means(l10, rng)
            delta = float(l20.mean() - l10.mean())
            low, high = ci95(bootstrap_delta)
            probability_positive = float(np.mean(bootstrap_delta > 0))
            direction_probability = (
                probability_positive if delta >= 0 else 1.0 - probability_positive
            )
            loo_agreement = sign_agreement(leave_one_out_deltas(l10, l20), delta)
            direction_adjusted, direction_agreement, n_directions = direction_statistics(
                group, metric, delta
            )
            classification = evidence_classification(
                delta,
                low,
                high,
                direction_probability,
                loo_agreement,
                direction_agreement,
            )
            rows.append(
                {
                    **base,
                    "metric": metric_label,
                    "n_l10": len(l10),
                    "n_l20": len(l20),
                    "mean_occupancy_l10": float(l10.mean()),
                    "mean_occupancy_l20": float(l20.mean()),
                    "median_occupancy_l10": float(np.median(l10)),
                    "median_occupancy_l20": float(np.median(l20)),
                    "delta_l20_minus_l10": delta,
                    "delta_ci95_low": low,
                    "delta_ci95_high": high,
                    "bootstrap_probability_delta_positive": probability_positive,
                    "bootstrap_direction_probability": direction_probability,
                    "loo_sign_agreement": loo_agreement,
                    "direction_adjusted_delta": direction_adjusted,
                    "direction_sign_agreement": direction_agreement,
                    "n_common_directions": n_directions,
                    "ci_excludes_zero": bool(low > 0 or high < 0),
                    "evidence_class": classification,
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["metric", "protein_id", "residue_number", "criterion"]
    )
