"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



from typing import Iterable

import numpy as np

import pandas as pd

from sklearn.mixture import GaussianMixture

FRAME_STEP_PS = 10.0

REFERENCE_ENVELOPE_DEG = 30.0

DWELL_NS_VALUES = (0.2, 0.5, 1.0)

MIXTURE_SAMPLE_PS = (100, 200, 500)

def normalize_angle(angle: float | np.ndarray) -> float | np.ndarray:
    values = (np.asarray(angle, dtype=float) + 180.0) % 360.0 - 180.0
    if np.isscalar(angle):
        return float(values)
    return values

def circular_mean_deg(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=float)
    radians = np.deg2rad(array)
    vector = np.mean(np.exp(1j * radians))
    return float(normalize_angle(np.rad2deg(np.angle(vector))))

def circular_sd_deg(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=float)
    radians = np.deg2rad(array)
    resultant = float(np.abs(np.mean(np.exp(1j * radians))))
    resultant = min(max(resultant, np.finfo(float).tiny), 1.0)
    return float(np.rad2deg(np.sqrt(-2.0 * np.log(resultant))))

def persistent_core_transitions(
    states: Iterable[str], min_dwell_frames: int
) -> tuple[int, int, str]:
    values = list(states)
    if not values:
        return 0, 0, ""
    runs: list[tuple[str, int]] = []
    current = values[0]
    count = 1
    for value in values[1:]:
        if value == current:
            count += 1
        else:
            runs.append((current, count))
            current = value
            count = 1
    runs.append((current, count))
    persistent = [
        state
        for state, length in runs
        if state in {"A", "B", "C"} and length >= min_dwell_frames
    ]
    collapsed: list[str] = []
    for state in persistent:
        if not collapsed or collapsed[-1] != state:
            collapsed.append(state)
    transitions = sum(
        current_state != previous_state
        for previous_state, current_state in zip(collapsed, collapsed[1:])
    )
    return transitions, len(set(collapsed)), ";".join(collapsed)

def safe_conditional_mean(
    frame: pd.DataFrame, metric: str, indicator: str, value: int
) -> tuple[float, int]:
    subset = frame.loc[frame[indicator].eq(value), metric]
    if subset.empty:
        return float("nan"), 0
    return float(subset.mean()), int(len(subset))

def summarize_trajectory(group: pd.DataFrame) -> dict[str, object]:
    row: dict[str, object] = {
        "case": group["case"].iloc[0],
        "pet_kind": group["pet_kind"].iloc[0],
        "direction": group["direction"].iloc[0],
        "replica": int(group["replica"].iloc[0]),
        "n_frames": len(group),
        "chi1_circular_mean_deg": circular_mean_deg(group["chi1_deg"]),
        "chi1_circular_sd_deg": circular_sd_deg(group["chi1_deg"]),
        "chi2_circular_mean_deg": circular_mean_deg(group["chi2_deg"]),
        "chi2_circular_sd_deg": circular_sd_deg(group["chi2_deg"]),
        "chi2_q05_deg": float(group["chi2_deg"].quantile(0.05)),
        "chi2_q95_deg": float(group["chi2_deg"].quantile(0.95)),
        "chi2_q90_span_deg": float(
            group["chi2_deg"].quantile(0.95) - group["chi2_deg"].quantile(0.05)
        ),
        "nearest_phl7_reference_distance_mean_deg": float(
            group["nearest_phl7_reference_distance_deg"].mean()
        ),
        "nearest_phl7_reference_distance_q95_deg": float(
            group["nearest_phl7_reference_distance_deg"].quantile(0.95)
        ),
        "within_30deg_reference_envelope_fraction": float(
            group["within_30deg_reference_envelope"].mean()
        ),
        "chi1_trans_fraction": float(group["chi1_trans_fraction_indicator"].mean()),
        "w156_contact_fraction": float(group["w156_any_pet"].mean()),
        "same_unit_coengagement_fraction": float(
            group["w156_with_f63_or_i179_same_unit"].mean()
        ),
        "complete_same_unit_fraction": float(
            group["f63_w156_i179_same_unit"].mean()
        ),
    }
    for state in ("A", "B", "C"):
        row[f"nearest_{state}_fraction"] = float(
            group["nearest_ispetase_state"].eq(state).mean()
        )
        row[f"core_{state}_fraction"] = float(
            group["reference_core_state"].eq(state).mean()
        )
    row["core_intermediate_fraction"] = float(
        group["reference_core_state"].eq("intermediate").mean()
    )
    for dwell_ns in DWELL_NS_VALUES:
        frames = int(round(dwell_ns * 1000.0 / FRAME_STEP_PS))
        transitions, visited, sequence = persistent_core_transitions(
            group["reference_core_state"], frames
        )
        suffix = str(dwell_ns).replace(".", "_")
        row[f"persistent_transitions_{suffix}ns"] = transitions
        row[f"persistent_states_visited_{suffix}ns"] = visited
        row[f"persistent_state_sequence_{suffix}ns"] = sequence
        row[f"persistent_transition_rate_per_100ns_{suffix}ns"] = (
            transitions / 80.0 * 100.0
        )
    conditional_specs = {
        "w156_contact": "w156_any_pet",
        "same_unit": "w156_with_f63_or_i179_same_unit",
    }
    for prefix, indicator in conditional_specs.items():
        for value, label in ((1, "yes"), (0, "no")):
            mean, count = safe_conditional_mean(
                group, "chi2_deg", indicator, value
            )
            row[f"chi2_mean_{prefix}_{label}_deg"] = mean
            row[f"frames_{prefix}_{label}"] = count
            subset = group.loc[group[indicator].eq(value)]
            for state in ("A", "B", "C"):
                row[f"nearest_{state}_fraction_{prefix}_{label}"] = (
                    float(subset["nearest_ispetase_state"].eq(state).mean())
                    if len(subset)
                    else float("nan")
                )
        yes = row[f"chi2_mean_{prefix}_yes_deg"]
        no = row[f"chi2_mean_{prefix}_no_deg"]
        row[f"chi2_mean_{prefix}_yes_minus_no_deg"] = (
            float(yes) - float(no)
            if np.isfinite(yes) and np.isfinite(no)
            else float("nan")
        )
    return row

def trajectory_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = [
        summarize_trajectory(group)
        for _, group in frame.groupby("case", sort=True, observed=True)
    ]
    return pd.DataFrame(rows)

BOOTSTRAP_METRICS = (
    "chi2_circular_mean_deg",
    "chi2_circular_sd_deg",
    "chi2_q90_span_deg",
    "nearest_phl7_reference_distance_mean_deg",
    "nearest_A_fraction",
    "nearest_B_fraction",
    "nearest_C_fraction",
    "core_A_fraction",
    "core_B_fraction",
    "core_C_fraction",
    "core_intermediate_fraction",
    "persistent_transition_rate_per_100ns_0_5ns",
    "w156_contact_fraction",
    "same_unit_coengagement_fraction",
)

def bootstrap_group_and_difference(
    summary: pd.DataFrame,
    resamples: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    group_rows: list[dict[str, object]] = []
    difference_rows: list[dict[str, object]] = []
    by_length = {
        pet_kind: group.reset_index(drop=True)
        for pet_kind, group in summary.groupby("pet_kind", sort=True)
    }
    if set(by_length) != {"PET_L10", "PET_L20"}:
        raise ValueError(f"Unexpected PET groups: {sorted(by_length)}")
    for metric in BOOTSTRAP_METRICS:
        for pet_kind, group in by_length.items():
            values = group[metric].dropna().to_numpy(dtype=float)
            if not len(values):
                continue
            draws = rng.choice(values, size=(resamples, len(values)), replace=True)
            means = draws.mean(axis=1)
            group_rows.append(
                {
                    "metric": metric,
                    "pet_kind": pet_kind,
                    "n_trajectories": len(values),
                    "mean": float(values.mean()),
                    "bootstrap_95_low": float(np.quantile(means, 0.025)),
                    "bootstrap_95_high": float(np.quantile(means, 0.975)),
                    "resamples": resamples,
                    "seed": seed,
                }
            )
        l10 = by_length["PET_L10"][metric].dropna().to_numpy(dtype=float)
        l20 = by_length["PET_L20"][metric].dropna().to_numpy(dtype=float)
        draws_l10 = rng.choice(l10, size=(resamples, len(l10)), replace=True).mean(axis=1)
        draws_l20 = rng.choice(l20, size=(resamples, len(l20)), replace=True).mean(axis=1)
        differences = draws_l20 - draws_l10
        difference_rows.append(
            {
                "metric": metric,
                "n_l10": len(l10),
                "n_l20": len(l20),
                "l10_mean": float(l10.mean()),
                "l20_mean": float(l20.mean()),
                "l20_minus_l10": float(l20.mean() - l10.mean()),
                "bootstrap_95_low": float(np.quantile(differences, 0.025)),
                "bootstrap_95_high": float(np.quantile(differences, 0.975)),
                "resamples": resamples,
                "seed": seed,
            }
        )
    return pd.DataFrame(group_rows), pd.DataFrame(difference_rows)

def mixture_model_sensitivity(frame: pd.DataFrame, seed: int) -> pd.DataFrame:
    """Compare one to three descriptive chi1/chi2 Gaussian components by BIC."""
    rows: list[dict[str, object]] = []
    for case, group in frame.groupby("case", sort=True, observed=True):
        ordered = group.sort_values("time_ps")
        for sample_ps in MIXTURE_SAMPLE_PS:
            stride = int(round(sample_ps / FRAME_STEP_PS))
            sampled = ordered.iloc[::stride]
            features = sampled[
                ["chi1_delta_from_trans_deg", "chi2_deg"]
            ].to_numpy(dtype=float)
            case_rows: list[dict[str, object]] = []
            for components in (1, 2, 3):
                model = GaussianMixture(
                    n_components=components,
                    covariance_type="full",
                    n_init=10,
                    random_state=seed,
                )
                model.fit(features)
                order = np.argsort(model.means_[:, 1])
                means = model.means_[order]
                weights = model.weights_[order]
                component_separation = (
                    float(np.max(means[:, 1]) - np.min(means[:, 1]))
                    if components > 1
                    else 0.0
                )
                case_rows.append(
                    {
                        "case": case,
                        "pet_kind": sampled["pet_kind"].iloc[0],
                        "direction": sampled["direction"].iloc[0],
                        "replica": int(sampled["replica"].iloc[0]),
                        "sample_interval_ps": sample_ps,
                        "n_sampled_frames": len(sampled),
                        "components": components,
                        "bic": float(model.bic(features)),
                        "component_chi1_delta_means_deg": ";".join(
                            f"{value:.4f}" for value in means[:, 0]
                        ),
                        "component_chi2_means_deg": ";".join(
                            f"{value:.4f}" for value in means[:, 1]
                        ),
                        "component_weights": ";".join(
                            f"{value:.6f}" for value in weights
                        ),
                        "chi2_component_span_deg": component_separation,
                        "all_component_weights_ge_0_10": bool(
                            np.all(weights >= 0.10)
                        ),
                    }
                )
            selected_index = int(np.argmin([row["bic"] for row in case_rows]))
            selected_components = int(case_rows[selected_index]["components"])
            for index, row in enumerate(case_rows):
                row["bic_delta_from_selected"] = float(
                    row["bic"] - case_rows[selected_index]["bic"]
                )
                row["selected_by_bic"] = index == selected_index
                row["resolved_multibasin_by_30deg_rule"] = bool(
                    index == selected_index
                    and selected_components > 1
                    and row["chi2_component_span_deg"] >= REFERENCE_ENVELOPE_DEG
                    and row["all_component_weights_ge_0_10"]
                )
                rows.append(row)
    return pd.DataFrame(rows)
