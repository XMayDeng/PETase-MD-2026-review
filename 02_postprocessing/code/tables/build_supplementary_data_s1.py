#!/usr/bin/env python3
"""Export the executed residue-to-contact-class map as Supplementary Data S1."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[3]
SOURCE = REPO / "02_postprocessing/inputs/annotations/full_panel_contact_landmark_annotation.csv"
# The verification entry point redirects this output to a temporary directory.
OUTPUT = REPO / "02_postprocessing/results/tables/Supplementary_Data_S1_residue_class_assignments.csv"
EXPECTED_SOURCE_SHA256 = (
    "dd195bf6936675d448f241fc6b73ef9d398dc3288ab715183c37c35c35e8dfef"
)
EXPECTED_ROWS = 1532

PROTEIN_METADATA = {
    "Pro00057": ("TfCut1", "7QJR", 1),
    "Pro00062": ("LCC", "4EB0", 2),
    "Pro00075": ("FoCut5a", "5AJH", 3),
    "Pro00083": ("IsPETase", "6EQD", 0),
    "Pro00121": ("HiC", "4OYY", 4),
    "Pro00137": ("PHL7", "7NEI", 5),
}

CLASS_SET = {
    "Catalytic/oxyanion",
    "W-loop/flexible or rigid pair",
    "Cleft/rim patch",
    "Flap/binding loop",
    "Subsite/hotspot/stability",
    "Other/candidate",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(output_path=None, check=False) -> None:
    if output_path is None:
        output_path = OUTPUT
    source_sha = sha256_file(SOURCE)
    if source_sha != EXPECTED_SOURCE_SHA256:
        raise RuntimeError(
            f"Source hash changed: expected {EXPECTED_SOURCE_SHA256}, got {source_sha}"
        )

    source = pd.read_csv(SOURCE)
    required = {
        "protein_id",
        "residue",
        "residue_name",
        "residue_label",
        "structural_region",
        "landmark_class",
        "annotation_confidence",
        "annotation_evidence_basis",
        "nearby_landmark",
        "nearby_landmark_residue",
        "nearby_landmark_role",
        "nearby_landmark_class",
        "nearest_landmark_ca_distance_a",
    }
    missing = required.difference(source.columns)
    if missing:
        raise RuntimeError(f"Missing source columns: {sorted(missing)}")
    if len(source) != EXPECTED_ROWS:
        raise RuntimeError(f"Expected {EXPECTED_ROWS} rows, found {len(source)}")
    if source.duplicated(["protein_id", "residue"]).any():
        raise RuntimeError("Duplicate protein-residue assignments found")
    if set(source["protein_id"]) != set(PROTEIN_METADATA):
        raise RuntimeError("Protein set differs from the registered six-enzyme panel")
    if set(source["landmark_class"]) != CLASS_SET:
        raise RuntimeError("Contact-class set differs from the registered six classes")

    output = source[
        [
            "protein_id",
            "residue",
            "residue_name",
            "residue_label",
            "structural_region",
            "landmark_class",
            "annotation_confidence",
            "annotation_evidence_basis",
            "nearby_landmark",
            "nearby_landmark_residue",
            "nearby_landmark_role",
            "nearby_landmark_class",
            "nearest_landmark_ca_distance_a",
        ]
    ].copy()
    output.insert(1, "enzyme", output["protein_id"].map(lambda x: PROTEIN_METADATA[x][0]))
    output.insert(2, "pdb_reference", output["protein_id"].map(lambda x: PROTEIN_METADATA[x][1]))
    output["panel_order"] = output["protein_id"].map(lambda x: PROTEIN_METADATA[x][2])
    output = output.sort_values(["panel_order", "residue"]).drop(columns="panel_order")
    output = output.rename(
        columns={
            "residue": "simulation_residue_number",
            "landmark_class": "assigned_contact_class",
            "annotation_evidence_basis": "assignment_basis",
            "nearby_landmark": "nearest_configured_landmark",
            "nearby_landmark_residue": "nearest_configured_landmark_number",
            "nearby_landmark_role": "nearest_configured_landmark_role",
            "nearby_landmark_class": "nearest_configured_landmark_class",
            "nearest_landmark_ca_distance_a": "nearest_landmark_ca_distance_A",
        }
    )
    output["nearest_landmark_ca_distance_A"] = output[
        "nearest_landmark_ca_distance_A"
    ].round(6)
    content = output.to_csv(index=False).encode('utf-8')
    if output_path.exists():
        if output_path.read_bytes() != content:
            raise ValueError(f'Existing output differs: {output_path}')
    elif check:
        raise FileNotFoundError(output_path)
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open('xb') as handle:
            handle.write(content)

    print(f"source={SOURCE}")
    print(f"source_sha256={source_sha}")
    print(f"output={output_path}")
    print(f"output_sha256={sha256_file(output_path)}")
    print(f"rows={len(output)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--check', action='store_true')
    arguments = parser.parse_args()
    main(arguments.output, arguments.check)
