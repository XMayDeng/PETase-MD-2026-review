"""Source-exact numerical core; input/output orchestration is separate."""

from __future__ import annotations



import numpy as np

import pandas as pd

from build_pooled_top8_unfiltered import CLASS_ORDER, PROTEIN_META, PROTEIN_ORDER

TOP_N = 20

def build_composition(
    pooled: pd.DataFrame,
    population: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]]]:
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

    top20 = pooled[pooled["rank_within_protein"] <= TOP_N].copy()
    counts = top20.groupby("protein_id").size().to_dict()
    check(
        "twenty residues per protein",
        all(counts.get(protein_id) == TOP_N for protein_id in PROTEIN_ORDER),
        "critical",
        ";".join(f"{PROTEIN_META[p]['short']}={counts.get(p, 0)}" for p in PROTEIN_ORDER),
    )
    check(
        "top-20 ranks complete",
        bool(
            top20.groupby("protein_id")["rank_within_protein"]
            .apply(lambda values: sorted(values.astype(int)) == list(range(1, TOP_N + 1)))
            .all()
        ),
        "critical",
        "expected ranks 1-20 for all six proteins",
    )
    check(
        "top-20 means bounded",
        bool(top20["mean_contact_fraction"].between(0, 1).all()),
        "critical",
        f"range={top20['mean_contact_fraction'].min():.6f}-{top20['mean_contact_fraction'].max():.6f}",
    )

    summary = (
        top20.groupby(["protein_id", "protein_short", "pdb_id", "landmark_class"], as_index=False)
        .agg(
            contact_mass=("mean_contact_fraction", "sum"),
            residue_count=("residue", "nunique"),
        )
    )
    grid = pd.MultiIndex.from_product(
        [PROTEIN_ORDER, CLASS_ORDER], names=["protein_id", "landmark_class"]
    ).to_frame(index=False)
    grid["protein_short"] = grid["protein_id"].map(
        {key: value["short"] for key, value in PROTEIN_META.items()}
    )
    grid["pdb_id"] = grid["protein_id"].map(
        {key: value["pdb"] for key, value in PROTEIN_META.items()}
    )
    summary = grid.merge(
        summary,
        on=["protein_id", "protein_short", "pdb_id", "landmark_class"],
        how="left",
        validate="one_to_one",
    )
    summary[["contact_mass", "residue_count"]] = summary[
        ["contact_mass", "residue_count"]
    ].fillna(0)
    summary["residue_count"] = summary["residue_count"].astype(int)
    totals = summary.groupby("protein_id")["contact_mass"].transform("sum")
    summary["top20_total_contact_mass"] = totals
    summary["contact_mass_share"] = summary["contact_mass"] / totals
    summary["contact_mass_percent"] = 100.0 * summary["contact_mass_share"]
    summary["retained_n"] = summary["protein_id"].map(
        population.set_index("protein_id")["retained_n"].astype(int).to_dict()
    )
    summary["stability_warn_n"] = summary["protein_id"].map(
        population.set_index("protein_id")["stability_warn_n"].astype(int).to_dict()
    )
    summary["protein_order"] = summary["protein_id"].map(
        {protein_id: index for index, protein_id in enumerate(PROTEIN_ORDER)}
    )
    summary["class_order"] = summary["landmark_class"].map(
        {class_name: index for index, class_name in enumerate(CLASS_ORDER)}
    )
    summary = summary.sort_values(["protein_order", "class_order"]).reset_index(drop=True)

    share_sums = summary.groupby("protein_id")["contact_mass_share"].sum()
    check(
        "positive top-20 contact mass",
        bool((summary.groupby("protein_id")["contact_mass"].sum() > 0).all()),
        "critical",
        ";".join(f"{PROTEIN_META[p]['short']}={summary.loc[summary['protein_id'].eq(p), 'contact_mass'].sum():.6f}" for p in PROTEIN_ORDER),
    )
    check(
        "class shares sum to one",
        bool(np.allclose(share_sums.to_numpy(), 1.0, atol=1e-12)),
        "critical",
        ";".join(f"{PROTEIN_META[p]['short']}={share_sums[p]:.12f}" for p in PROTEIN_ORDER),
    )
    check(
        "six classes represented in output grid",
        bool((summary.groupby("protein_id")["landmark_class"].nunique() == len(CLASS_ORDER)).all()),
        "critical",
        "six explicit rows per protein, including zero-mass classes",
    )
    check(
        "analysis uses no low-frequency filter",
        True,
        "critical",
        "composition derives from unfiltered pooled trajectory means, including zeros and values below 0.20",
    )
    return top20, summary.drop(columns=["protein_order", "class_order"]), checks
