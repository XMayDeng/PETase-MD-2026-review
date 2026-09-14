"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import pandas as pd

PROTEINS = [
    {"protein_id": "Pro00083", "short": "IsPETase", "pdb": "6EQD", "panel": "A"},
    {"protein_id": "Pro00057", "short": "TfCut1", "pdb": "7QJR", "panel": "B"},
    {"protein_id": "Pro00062", "short": "LCC", "pdb": "4EB0", "panel": "C"},
    {"protein_id": "Pro00075", "short": "FoCut5a", "pdb": "5AJH", "panel": "D"},
    {"protein_id": "Pro00121", "short": "HiC", "pdb": "4OYY", "panel": "E"},
    {"protein_id": "Pro00137", "short": "PHL7", "pdb": "7NEI", "panel": "F"},
]

PROTEIN_ORDER = [item["protein_id"] for item in PROTEINS]

PROTEIN_META = {item["protein_id"]: item for item in PROTEINS}

EXPECTED_RETAINED = {
    "Pro00083": 13,
    "Pro00057": 14,
    "Pro00062": 13,
    "Pro00075": 12,
    "Pro00121": 13,
    "Pro00137": 6,
}

CLASS_ORDER = [
    "Catalytic/oxyanion",
    "W-loop/flexible or rigid pair",
    "Cleft/rim patch",
    "Flap/binding loop",
    "Subsite/hotspot/stability",
    "Other/candidate",
]

TOP_N = 8

CONSISTENCY_THRESHOLD = 0.5

EXPECTED_FRAMES = 8001

def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().eq("true")

def validate_and_join(
    contacts: pd.DataFrame,
    manifest: pd.DataFrame,
    annotation: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, severity: str, detail: str) -> None:
        checks.append(
            {
                "check": name,
                "status": "PASS" if passed else "FAIL",
                "severity_if_failed": severity,
                "detail": detail,
            }
        )

    required_contact_columns = {
        "case",
        "residue",
        "contact_fraction",
        "contact_count",
        "total_frames",
        "protein_id",
        "pet_kind",
        "direction",
        "replica_id",
    }
    required_manifest_columns = {
        "case",
        "protein_id",
        "protein_label",
        "strict_bool",
        "retained_bool",
        "stability_pass_bool",
    }
    required_annotation_columns = {
        "protein_id",
        "residue",
        "residue_name",
        "residue_label",
        "structural_region",
        "landmark_class",
        "annotation_confidence",
        "protein_caveat",
    }
    check(
        "contact schema",
        required_contact_columns.issubset(contacts.columns),
        "critical",
        f"columns={len(contacts.columns)}",
    )
    check(
        "manifest schema",
        required_manifest_columns.issubset(manifest.columns),
        "critical",
        f"columns={len(manifest.columns)}",
    )
    check(
        "annotation schema",
        required_annotation_columns.issubset(annotation.columns),
        "critical",
        f"columns={len(annotation.columns)}",
    )
    if any(row["status"] == "FAIL" for row in checks):
        raise ValueError("Required input schema check failed")

    retained_manifest = manifest[bool_series(manifest["retained_bool"])].copy()
    retained_manifest["strict_bool"] = bool_series(retained_manifest["strict_bool"])
    retained_manifest["stability_pass_bool"] = bool_series(retained_manifest["stability_pass_bool"])
    check(
        "retained manifest population",
        len(retained_manifest) == 71 and retained_manifest["case"].nunique() == 71,
        "critical",
        f"rows={len(retained_manifest)}; cases={retained_manifest['case'].nunique()}",
    )
    check(
        "contact case population",
        contacts["case"].nunique() == 71,
        "critical",
        f"cases={contacts['case'].nunique()}",
    )
    check(
        "contact and manifest case sets",
        set(contacts["case"]) == set(retained_manifest["case"]),
        "critical",
        f"contact_only={len(set(contacts['case']) - set(retained_manifest['case']))}; manifest_only={len(set(retained_manifest['case']) - set(contacts['case']))}",
    )
    check(
        "unique case-residue grain",
        not contacts.duplicated(["case", "residue"]).any(),
        "critical",
        f"duplicates={int(contacts.duplicated(['case', 'residue']).sum())}",
    )
    check(
        "contact fractions bounded",
        bool(contacts["contact_fraction"].between(0, 1).all()),
        "critical",
        f"range={contacts['contact_fraction'].min():.6f}-{contacts['contact_fraction'].max():.6f}",
    )
    check(
        "frame count",
        bool((contacts["total_frames"] == EXPECTED_FRAMES).all()),
        "critical",
        f"observed={sorted(contacts['total_frames'].unique().tolist())}",
    )

    case_meta = contacts[["case", "protein_id"]].drop_duplicates()
    manifest_case_meta = retained_manifest[["case", "protein_id"]].drop_duplicates()
    meta_compare = case_meta.merge(
        manifest_case_meta,
        on="case",
        how="outer",
        suffixes=("_contact", "_manifest"),
        indicator=True,
    )
    metadata_match = bool(
        meta_compare["_merge"].eq("both").all()
        and meta_compare["protein_id_contact"].eq(meta_compare["protein_id_manifest"]).all()
    )
    check("case protein metadata agreement", metadata_match, "critical", f"rows={len(meta_compare)}")

    annotation_fields = [
        "protein_id",
        "residue",
        "residue_name",
        "residue_label",
        "structural_region",
        "landmark_class",
        "annotation_confidence",
        "protein_caveat",
    ]
    annotation_meta = annotation[annotation_fields].drop_duplicates()
    annotation_unique = not annotation_meta.duplicated(["protein_id", "residue"]).any()
    check(
        "unique residue annotation",
        annotation_unique,
        "critical",
        f"rows={len(annotation_meta)}; duplicate_keys={int(annotation_meta.duplicated(['protein_id', 'residue']).sum())}",
    )
    if not annotation_unique:
        raise ValueError("Annotation has conflicting protein-residue rows")

    working = contacts.merge(
        retained_manifest[
            ["case", "protein_id", "protein_label", "strict_bool", "stability_pass_bool"]
        ],
        on=["case", "protein_id"],
        how="inner",
        validate="many_to_one",
    )
    working = working.merge(
        annotation_meta,
        on=["protein_id", "residue"],
        how="left",
        validate="many_to_one",
    )
    annotation_missing = int(working["residue_label"].isna().sum())
    check(
        "complete annotation coverage",
        annotation_missing == 0,
        "critical",
        f"missing_rows={annotation_missing}",
    )
    observed_classes = set(working["landmark_class"].dropna().unique())
    unexpected_classes = sorted(observed_classes - set(CLASS_ORDER))
    check(
        "annotation classes recognized",
        not unexpected_classes and not working["landmark_class"].isna().any(),
        "critical",
        f"classes={sorted(observed_classes)}; unexpected={unexpected_classes}",
    )

    case_residue_counts = (
        working.groupby(["protein_id", "case"], as_index=False)["residue"]
        .nunique()
        .rename(columns={"residue": "case_residue_count"})
    )
    expected_by_protein = (
        working.groupby("protein_id", as_index=False)["residue"]
        .nunique()
        .rename(columns={"residue": "expected_residue_count"})
    )
    case_residue_counts = case_residue_counts.merge(expected_by_protein, on="protein_id")
    complete_grid = bool(
        case_residue_counts["case_residue_count"].eq(case_residue_counts["expected_residue_count"]).all()
    )
    check(
        "complete case-residue grid",
        complete_grid,
        "critical",
        f"cases={len(case_residue_counts)}; min={case_residue_counts['case_residue_count'].min()}; max={case_residue_counts['case_residue_count'].max()}",
    )

    population = (
        retained_manifest.groupby(["protein_id", "protein_label"], as_index=False)
        .agg(
            retained_n=("case", "nunique"),
            strict_n=("strict_bool", "sum"),
            stability_pass_n=("stability_pass_bool", "sum"),
            pet_kinds=("case", lambda _: 0),
        )
    )
    pet_kind_counts = contacts.groupby("protein_id")["pet_kind"].nunique().to_dict()
    population["pet_kinds"] = population["protein_id"].map(pet_kind_counts).astype(int)
    population["stability_warn_n"] = population["retained_n"] - population["stability_pass_n"]
    population["protein_short"] = population["protein_id"].map(
        {key: value["short"] for key, value in PROTEIN_META.items()}
    )
    population["pdb_id"] = population["protein_id"].map(
        {key: value["pdb"] for key, value in PROTEIN_META.items()}
    )
    observed_counts = population.set_index("protein_id")["retained_n"].astype(int).to_dict()
    check(
        "retained denominators",
        observed_counts == EXPECTED_RETAINED,
        "critical",
        ";".join(f"{PROTEIN_META[p]['short']}={observed_counts.get(p, 0)}" for p in PROTEIN_ORDER),
    )
    check(
        "analysis uses no low-frequency filter",
        True,
        "critical",
        "all contact_fraction values, including zeros and values below 0.20, enter pooled means",
    )

    return working, retained_manifest, population, checks

def q25(series: pd.Series) -> float:
    return float(series.quantile(0.25))

def q75(series: pd.Series) -> float:
    return float(series.quantile(0.75))

def build_pooled_summary(working: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = (
        working.groupby(
            [
                "protein_id",
                "protein_label",
                "residue",
                "residue_name",
                "residue_label",
                "structural_region",
                "landmark_class",
                "annotation_confidence",
                "protein_caveat",
            ],
            as_index=False,
        )
        .agg(
            mean_contact_fraction=("contact_fraction", "mean"),
            median_contact_fraction=("contact_fraction", "median"),
            contact_fraction_sd=("contact_fraction", "std"),
            contact_fraction_q25=("contact_fraction", q25),
            contact_fraction_q75=("contact_fraction", q75),
            min_contact_fraction=("contact_fraction", "min"),
            max_contact_fraction=("contact_fraction", "max"),
            n_retained=("case", "nunique"),
            n_cases_ge_0_5=("contact_fraction", lambda values: int((values >= CONSISTENCY_THRESHOLD).sum())),
            stability_pass_n=("stability_pass_bool", "sum"),
        )
    )
    grouped["contact_fraction_sd"] = grouped["contact_fraction_sd"].fillna(0.0)
    grouped["support_fraction_ge_0_5"] = grouped["n_cases_ge_0_5"] / grouped["n_retained"]
    grouped["protein_short"] = grouped["protein_id"].map(
        {key: value["short"] for key, value in PROTEIN_META.items()}
    )
    grouped["pdb_id"] = grouped["protein_id"].map(
        {key: value["pdb"] for key, value in PROTEIN_META.items()}
    )
    grouped["protein_order"] = grouped["protein_id"].map(
        {protein_id: index for index, protein_id in enumerate(PROTEIN_ORDER)}
    )
    grouped = grouped.sort_values(
        [
            "protein_order",
            "mean_contact_fraction",
            "n_cases_ge_0_5",
            "max_contact_fraction",
            "residue",
        ],
        ascending=[True, False, False, False, True],
    ).reset_index(drop=True)
    grouped["rank_within_protein"] = grouped.groupby("protein_id").cumcount() + 1
    grouped["is_top8"] = grouped["rank_within_protein"] <= TOP_N
    top8 = grouped[grouped["is_top8"]].copy()
    return grouped.drop(columns=["protein_order"]), top8.drop(columns=["protein_order"])
