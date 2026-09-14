"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import pandas as pd

PET_ORDER = ["PET_L4", "PET_L10", "PET_L20"]

PROTEIN_ORDER = ["Pro00057", "Pro00062", "Pro00075", "Pro00083", "Pro00121", "Pro00137"]

SHORT_LABELS = {
    "Pro00057": "CUT1",
    "Pro00062": "LCC",
    "Pro00075": "FoCut5a",
    "Pro00083": "IsPETase",
    "Pro00121": "HiC",
    "Pro00137": "PHL7",
}

CLASS_ORDER = [
    "Catalytic/oxyanion",
    "W-loop/flexible or rigid pair",
    "Flap/binding loop",
    "Cleft/rim patch",
    "Subsite/hotspot/stability",
    "Other/candidate",
]

def build_share_table(
    contacts: pd.DataFrame,
    manifest: pd.DataFrame,
    annotation: pd.DataFrame,
    *,
    min_contact_fraction: float | None,
    analysis_label: str,
) -> pd.DataFrame:
    class_map = annotation[["protein_id", "residue", "landmark_class"]].drop_duplicates()
    working = contacts.merge(
        manifest[["case", "protein_id", "protein_label", "pet_kind", "direction", "strict_bool", "stability_pass_bool"]],
        on="case",
        how="inner",
    )
    working = working.merge(class_map, on=["protein_id", "residue"], how="left")
    working["landmark_class"] = working["landmark_class"].fillna("Other/candidate")
    if min_contact_fraction is not None:
        working = working[working["contact_fraction"] >= min_contact_fraction].copy()

    grouped = (
        working.groupby(["protein_id", "protein_label", "pet_kind", "landmark_class"], as_index=False)
        .agg(
            contact_mass=("contact_fraction", "sum"),
            residue_case_signals=("contact_fraction", "count"),
            mean_contact_fraction_included=("contact_fraction", "mean"),
        )
    )
    totals = grouped.groupby(["protein_id", "pet_kind"])["contact_mass"].transform("sum")
    grouped["contact_mass_share"] = grouped["contact_mass"] / totals
    grouped["analysis"] = analysis_label
    grouped["min_contact_fraction"] = -1.0 if min_contact_fraction is None else min_contact_fraction
    grouped["protein_short"] = grouped["protein_id"].map(SHORT_LABELS).fillna(grouped["protein_id"])
    grouped["pet_kind"] = pd.Categorical(grouped["pet_kind"], PET_ORDER, ordered=True)
    grouped["protein_id"] = pd.Categorical(grouped["protein_id"], PROTEIN_ORDER, ordered=True)
    grouped["landmark_class"] = pd.Categorical(grouped["landmark_class"], CLASS_ORDER, ordered=True)
    return grouped.sort_values(["protein_id", "pet_kind", "landmark_class"]).reset_index(drop=True)
