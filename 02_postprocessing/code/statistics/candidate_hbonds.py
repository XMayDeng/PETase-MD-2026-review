"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import hashlib

import numpy as np

import pandas as pd

N_BOOTSTRAP = 20_000

BOOTSTRAP_SEED = 20260713

def stable_seed(protein_id: str, residue_number: int, suffix: str) -> int:
    token = f"{BOOTSTRAP_SEED}|{protein_id}|{residue_number}|{suffix}".encode("ascii")
    return int.from_bytes(hashlib.sha256(token).digest()[:8], "little")

def bootstrap_means(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    indices = rng.integers(0, len(values), size=(N_BOOTSTRAP, len(values)))
    return values[indices].mean(axis=1)

def ci95(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, [0.025, 0.975])
    return float(low), float(high)

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
    if reference == 0 or not values:
        return float("nan")
    target = np.sign(reference)
    return float(np.mean([np.sign(value) == target for value in values if value != 0]))

def summarize_candidate_pet(case_values: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "any_hbond_occupancy",
        "sidechain_hbond_occupancy",
        "backbone_hbond_occupancy",
        "protein_donor_hbond_occupancy",
        "protein_acceptor_hbond_occupancy",
        "mean_triplet_count_per_frame",
        "heavy_atom_contact_fraction",
    ]
    rows: list[dict[str, object]] = []
    for keys, group in case_values.groupby(
        ["protein_id", "protein_name", "residue_number", "residue_label", "pet_kind"],
        sort=False,
    ):
        row: dict[str, object] = dict(
            zip(["protein_id", "protein_name", "residue_number", "residue_label", "pet_kind"], keys)
        )
        row["n_trajectories"] = len(group)
        for metric in metrics:
            values = group[metric].to_numpy(dtype=float)
            rng = np.random.default_rng(stable_seed(str(keys[0]), int(keys[2]), f"{keys[4]}|{metric}"))
            boot = bootstrap_means(values, rng)
            low, high = ci95(boot)
            row[f"mean_{metric}"] = float(values.mean())
            row[f"median_{metric}"] = float(np.median(values))
            row[f"mean_{metric}_ci95_low"] = low
            row[f"mean_{metric}_ci95_high"] = high
        occupancy = group["any_hbond_occupancy"].to_numpy(dtype=float)
        sidechain_occupancy = group["sidechain_hbond_occupancy"].to_numpy(dtype=float)
        for threshold in (0.01, 0.10, 0.50):
            token = str(threshold).replace(".", "_")
            row[f"n_trajectories_ge_{token}"] = int(np.sum(occupancy >= threshold))
            row[f"n_sidechain_trajectories_ge_{token}"] = int(
                np.sum(sidechain_occupancy >= threshold)
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["protein_id", "residue_number", "pet_kind"])

def direction_metrics(
    group: pd.DataFrame,
    overall_delta: float,
    metric: str,
) -> tuple[list[dict[str, object]], float, float]:
    common = sorted(
        set(group.loc[group["pet_kind"].eq("PET_L10"), "direction"])
        & set(group.loc[group["pet_kind"].eq("PET_L20"), "direction"])
    )
    rows: list[dict[str, object]] = []
    deltas: list[float] = []
    for direction in common:
        l10 = group[
            group["pet_kind"].eq("PET_L10") & group["direction"].eq(direction)
        ][metric]
        l20 = group[
            group["pet_kind"].eq("PET_L20") & group["direction"].eq(direction)
        ][metric]
        delta = float(l20.mean() - l10.mean())
        deltas.append(delta)
        rows.append(
            {
                "direction": direction,
                "n_l10": len(l10),
                "n_l20": len(l20),
                "mean_occupancy_l10": float(l10.mean()),
                "mean_occupancy_l20": float(l20.mean()),
                "delta_l20_minus_l10": delta,
            }
        )
    adjusted = float(np.mean(deltas)) if deltas else float("nan")
    agreement = sign_agreement(deltas, overall_delta)
    return rows, adjusted, agreement

def evidence_classification(
    delta: float,
    low: float,
    high: float,
    direction_probability: float,
    loo_agreement: float,
    direction_agreement: float,
) -> tuple[str, bool]:
    ci_excludes_zero = low > 0 or high < 0
    robust_direction = bool(
        ci_excludes_zero
        and direction_probability >= 0.95
        and loo_agreement >= 0.80
        and direction_agreement >= 2.0 / 3.0
    )
    if robust_direction and delta > 0:
        return "resolved_L20_increase", ci_excludes_zero
    if robust_direction and delta < 0:
        return "resolved_L20_decrease", ci_excludes_zero
    if direction_probability >= 0.80:
        return "directional_but_unresolved", ci_excludes_zero
    return "no_resolved_hbond_shift", ci_excludes_zero

def analyze_shifts(case_values: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, object]] = []
    direction_records: list[dict[str, object]] = []
    for keys, group in case_values.groupby(
        ["protein_id", "protein_name", "residue_number", "residue_label"], sort=False
    ):
        record: dict[str, object] = {
            "protein_id": keys[0],
            "protein_name": keys[1],
            "residue_number": keys[2],
            "residue_label": keys[3],
            "spearman_heavy_contact_vs_any_hbond": float(
                group["heavy_atom_contact_fraction"].corr(
                    group["any_hbond_occupancy"], method="spearman"
                )
            ),
            "spearman_heavy_contact_vs_sidechain_hbond": float(
                group["heavy_atom_contact_fraction"].corr(
                    group["sidechain_hbond_occupancy"], method="spearman"
                )
            ),
        }
        for metric, output_token, direction_label in (
            ("any_hbond_occupancy", "", "any_residue_atom"),
            ("sidechain_hbond_occupancy", "sidechain_", "sidechain_only"),
        ):
            l10 = group.loc[group["pet_kind"].eq("PET_L10"), metric].to_numpy(float)
            l20 = group.loc[group["pet_kind"].eq("PET_L20"), metric].to_numpy(float)
            rng = np.random.default_rng(
                stable_seed(str(keys[0]), int(keys[2]), f"{direction_label}|delta")
            )
            delta_boot = bootstrap_means(l20, rng) - bootstrap_means(l10, rng)
            delta = float(l20.mean() - l10.mean())
            low, high = ci95(delta_boot)
            probability_positive = float(np.mean(delta_boot > 0))
            direction_probability = (
                probability_positive if delta >= 0 else 1.0 - probability_positive
            )
            loo_agreement = sign_agreement(leave_one_out_deltas(l10, l20), delta)
            direction_rows, adjusted, direction_agreement = direction_metrics(
                group, delta, metric
            )
            for row in direction_rows:
                direction_records.append(
                    {
                        "protein_id": keys[0],
                        "protein_name": keys[1],
                        "residue_number": keys[2],
                        "residue_label": keys[3],
                        "metric": direction_label,
                        **row,
                    }
                )
            evidence_class, ci_excludes_zero = evidence_classification(
                delta,
                low,
                high,
                direction_probability,
                loo_agreement,
                direction_agreement,
            )
            occupancy_token = f"{output_token}hbond_occupancy"
            record.update(
                {
                    "n_l10": len(l10),
                    "n_l20": len(l20),
                    f"mean_{occupancy_token}_l10": float(l10.mean()),
                    f"mean_{occupancy_token}_l20": float(l20.mean()),
                    f"median_{occupancy_token}_l10": float(np.median(l10)),
                    f"median_{occupancy_token}_l20": float(np.median(l20)),
                    f"{output_token}delta_l20_minus_l10": delta,
                    f"{output_token}delta_ci95_low": low,
                    f"{output_token}delta_ci95_high": high,
                    f"{output_token}bootstrap_probability_delta_positive": probability_positive,
                    f"{output_token}bootstrap_direction_probability": direction_probability,
                    f"{output_token}cliffs_delta_l20_vs_l10": cliffs_delta(l20, l10),
                    f"{output_token}loo_sign_agreement": loo_agreement,
                    f"{output_token}direction_adjusted_delta": adjusted,
                    f"{output_token}direction_sign_agreement": direction_agreement,
                    f"{output_token}n_common_directions": len(direction_rows),
                    f"{output_token}ci_excludes_zero": ci_excludes_zero,
                    f"{output_token}absolute_delta_ge_0_05": abs(delta) >= 0.05,
                    f"{output_token}evidence_class": evidence_class,
                }
            )
        records.append(record)
    return (
        pd.DataFrame(records).sort_values(["protein_id", "residue_number"]),
        pd.DataFrame(direction_records).sort_values(
            ["protein_id", "residue_number", "metric", "direction"]
        ),
    )
