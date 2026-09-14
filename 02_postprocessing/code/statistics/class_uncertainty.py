"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import math

import numpy as np

import pandas as pd

MIN_CONTACT_FRACTION = 0.20

N_BOOTSTRAP = 20_000

BASE_SEED = 20_260_723

PET_FROM = "PET_L10"

PET_TO = "PET_L20"

PROTEIN_ORDER = ["Pro00083", "Pro00057", "Pro00062", "Pro00075", "Pro00121", "Pro00137"]

PROTEIN_NAMES = {
    "Pro00083": "IsPETase",
    "Pro00057": "TfCut1",
    "Pro00062": "LCC",
    "Pro00075": "FoCut5a",
    "Pro00121": "HiC",
    "Pro00137": "PHL7",
}

CLASS_ORDER = [
    "Catalytic/oxyanion",
    "W-loop/flexible or rigid pair",
    "Cleft/rim patch",
    "Flap/binding loop",
    "Subsite/hotspot/stability",
    "Other/candidate",
]

STABILITY_WARNING_PROTEINS = {"Pro00083", "Pro00062", "Pro00075"}

def build_case_class_mass(
    contacts: pd.DataFrame,
    manifest: pd.DataFrame,
    annotation: pd.DataFrame,
) -> pd.DataFrame:
    if contacts.duplicated(["case", "residue"]).any():
        raise ValueError("Duplicate case/residue rows in contact table")
    if manifest.duplicated("case").any():
        raise ValueError("Duplicate cases in retained manifest")

    working = contacts.merge(manifest, on="case", how="inner", validate="many_to_one")
    working = working.merge(annotation, on=["protein_id", "residue"], how="left", validate="many_to_one")
    working["landmark_class"] = working["landmark_class"].fillna("Other/candidate")
    working["filtered_contact_mass"] = working["contact_fraction"].where(
        working["contact_fraction"] >= MIN_CONTACT_FRACTION,
        0.0,
    )

    grouped = (
        working.groupby(
            [
                "case",
                "protein_id",
                "protein_label",
                "pet_kind",
                "direction",
                "replica_id",
                "stability_pass_bool",
                "landmark_class",
            ],
            as_index=False,
        )["filtered_contact_mass"]
        .sum()
    )
    index_columns = [
        "case",
        "protein_id",
        "protein_label",
        "pet_kind",
        "direction",
        "replica_id",
        "stability_pass_bool",
    ]
    wide = (
        grouped.pivot(index=index_columns, columns="landmark_class", values="filtered_contact_mass")
        .reindex(columns=CLASS_ORDER, fill_value=0.0)
        .fillna(0.0)
        .reset_index()
    )
    wide["total_filtered_contact_mass"] = wide[CLASS_ORDER].sum(axis=1)
    if (wide["total_filtered_contact_mass"] <= 0).any():
        cases = ", ".join(wide.loc[wide["total_filtered_contact_mass"] <= 0, "case"].astype(str))
        raise ValueError(f"Retained trajectories with zero filtered contact mass: {cases}")
    return wide

def pooled_share(matrix: np.ndarray) -> np.ndarray:
    totals = matrix.sum(axis=0)
    denominator = float(totals.sum())
    if denominator <= 0:
        raise ValueError("Cannot normalize zero contact mass")
    return totals / denominator

def bootstrap_pooled_shares(matrix: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    indices = rng.integers(0, len(matrix), size=(N_BOOTSTRAP, len(matrix)))
    totals = matrix[indices].sum(axis=1)
    denominators = totals.sum(axis=1, keepdims=True)
    if np.any(denominators <= 0):
        raise ValueError("Bootstrap resample produced zero contact mass")
    return totals / denominators

def sign_agreement(values: list[float], reference: float) -> float:
    if not values or math.isclose(reference, 0.0, abs_tol=1e-15):
        return float("nan")
    signs = [np.sign(value) for value in values if not math.isclose(value, 0.0, abs_tol=1e-15)]
    if not signs:
        return float("nan")
    return float(np.mean(np.asarray(signs) == np.sign(reference)))

def leave_one_out_metrics(from_matrix: np.ndarray, to_matrix: np.ndarray, class_index: int, delta: float) -> tuple[float, float, float]:
    values: list[float] = []
    if len(from_matrix) > 1:
        for index in range(len(from_matrix)):
            values.append(float(pooled_share(to_matrix)[class_index] - pooled_share(np.delete(from_matrix, index, axis=0))[class_index]))
    if len(to_matrix) > 1:
        for index in range(len(to_matrix)):
            values.append(float(pooled_share(np.delete(to_matrix, index, axis=0))[class_index] - pooled_share(from_matrix)[class_index]))
    return sign_agreement(values, delta), float(min(values)), float(max(values))

def direction_metrics(case_mass: pd.DataFrame, protein_id: str, class_index: int, delta: float) -> tuple[int, str, str, float, float]:
    protein = case_mass[case_mass["protein_id"] == protein_id]
    from_directions = set(protein.loc[protein["pet_kind"] == PET_FROM, "direction"])
    to_directions = set(protein.loc[protein["pet_kind"] == PET_TO, "direction"])
    common = sorted(from_directions & to_directions)
    deltas: list[float] = []
    labels: list[str] = []
    for direction in common:
        from_matrix = protein[(protein["pet_kind"] == PET_FROM) & (protein["direction"] == direction)][CLASS_ORDER].to_numpy(float)
        to_matrix = protein[(protein["pet_kind"] == PET_TO) & (protein["direction"] == direction)][CLASS_ORDER].to_numpy(float)
        direction_delta = float(pooled_share(to_matrix)[class_index] - pooled_share(from_matrix)[class_index])
        deltas.append(direction_delta)
        labels.append(f"{direction}:{100.0 * direction_delta:+.4f}")
    adjusted = float(np.mean(deltas)) if deltas else float("nan")
    return len(common), ";".join(common), ";".join(labels), adjusted, sign_agreement(deltas, delta)

def analyze(case_mass: pd.DataFrame, figure7: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    primary = figure7[figure7["comparison"] == "PET_L20_minus_PET_L10"].copy()

    for protein_index, protein_id in enumerate(PROTEIN_ORDER):
        protein = case_mass[case_mass["protein_id"] == protein_id]
        from_rows = protein[protein["pet_kind"] == PET_FROM]
        to_rows = protein[protein["pet_kind"] == PET_TO]
        from_matrix = from_rows[CLASS_ORDER].to_numpy(float)
        to_matrix = to_rows[CLASS_ORDER].to_numpy(float)
        if len(from_matrix) == 0 or len(to_matrix) == 0:
            raise ValueError(f"Missing L10/L20 group for {protein_id}")

        from_share = pooled_share(from_matrix)
        to_share = pooled_share(to_matrix)
        rng = np.random.Generator(np.random.PCG64(BASE_SEED + 1009 * protein_index))
        from_bootstrap = bootstrap_pooled_shares(from_matrix, rng)
        to_bootstrap = bootstrap_pooled_shares(to_matrix, rng)
        bootstrap_delta = to_bootstrap - from_bootstrap

        for class_index, landmark_class in enumerate(CLASS_ORDER):
            delta = float(to_share[class_index] - from_share[class_index])
            bootstrap_values = bootstrap_delta[:, class_index]
            ci_low, ci_high = np.quantile(bootstrap_values, [0.025, 0.975])
            all_zero = bool(
                np.allclose(from_matrix[:, class_index], 0.0)
                and np.allclose(to_matrix[:, class_index], 0.0)
            )
            probability_gt_zero = float("nan") if all_zero else float(np.mean(bootstrap_values > 0))
            directional_probability = float("nan") if all_zero else max(probability_gt_zero, 1.0 - probability_gt_zero)
            loo_agreement, loo_min, loo_max = leave_one_out_metrics(
                from_matrix,
                to_matrix,
                class_index,
                delta,
            )
            common_n, common_labels, direction_deltas, adjusted_delta, direction_agreement = direction_metrics(
                case_mass,
                protein_id,
                class_index,
                delta,
            )
            ci_excludes_zero = bool(ci_low > 0 or ci_high < 0)
            small_n = min(len(from_matrix), len(to_matrix)) < 3
            if all_zero:
                evidence_label = "structural_zero_in_both_groups"
            elif small_n:
                evidence_label = "descriptive_small_n"
            elif ci_excludes_zero and loo_agreement >= 0.8 and direction_agreement >= 2 / 3:
                evidence_label = "exploratory_directionally_consistent"
            else:
                evidence_label = "uncertain_across_trajectories"

            expected = primary[
                (primary["protein_short"] == PROTEIN_NAMES[protein_id])
                & (primary["landmark_class"] == landmark_class)
            ]
            if len(expected) != 1:
                raise ValueError(f"Missing Figure 7 point estimate for {protein_id} / {landmark_class}")
            figure7_delta = float(expected.iloc[0]["delta_share_pct_points"])

            rows.append(
                {
                    "protein_id": protein_id,
                    "protein_short": PROTEIN_NAMES[protein_id],
                    "landmark_class": landmark_class,
                    "n_l10": len(from_matrix),
                    "n_l20": len(to_matrix),
                    "l10_pooled_share": from_share[class_index],
                    "l20_pooled_share": to_share[class_index],
                    "delta_share": delta,
                    "delta_pct_points": 100.0 * delta,
                    "bootstrap_ci95_low_pct_points": 100.0 * float(ci_low),
                    "bootstrap_ci95_high_pct_points": 100.0 * float(ci_high),
                    "bootstrap_probability_gt_zero": probability_gt_zero,
                    "bootstrap_directional_probability": directional_probability,
                    "ci95_excludes_zero": ci_excludes_zero,
                    "leave_one_out_sign_agreement": loo_agreement,
                    "leave_one_out_min_pct_points": 100.0 * loo_min,
                    "leave_one_out_max_pct_points": 100.0 * loo_max,
                    "common_direction_count": common_n,
                    "common_directions": common_labels,
                    "direction_deltas_pct_points": direction_deltas,
                    "direction_adjusted_delta_pct_points": 100.0 * adjusted_delta,
                    "direction_sign_agreement": direction_agreement,
                    "small_n_flag": small_n,
                    "stability_warning": protein_id in STABILITY_WARNING_PROTEINS,
                    "all_zero_both_groups": all_zero,
                    "evidence_label": evidence_label,
                    "figure7_delta_pct_points": figure7_delta,
                    "figure7_match_error_pct_points": 100.0 * delta - figure7_delta,
                }
            )
    return pd.DataFrame(rows)
