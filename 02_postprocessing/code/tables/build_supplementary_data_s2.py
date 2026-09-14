#!/usr/bin/env python3
"""Export all 95 CON-03 candidates without recalculating their statistics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
ANALYSIS_ROOT = REPO_ROOT / "02_postprocessing/results/aggregate_statistics"
SOURCE = ANALYSIS_ROOT / "trajectory_level_contact_shift_candidates.csv"
SUMMARY = ANALYSIS_ROOT / "manuscript_candidate_residues.csv"
TABLE = REPO_ROOT / "02_postprocessing/results/tables/Supplementary_Table_S2_residue_contact_shifts.tex"
OUTPUT = REPO_ROOT / "02_postprocessing/results/tables/Supplementary_Data_S2_residue_contact_shift_candidates.csv"
SOURCE_SHA256 = "4e4bcfb5a6e0b042f327efa6d996b3f99bf7ffb375db36ba61dd2fc8f63f6eda"
SUMMARY_SHA256 = "01d753098732516550f2f0112529545ecfd735258deb7ef013540779fa0b53ab"

PROTEINS = {
    "Pro00083": ("IsPETase", "6EQD", 0, 17),
    "Pro00057": ("TfCut1", "7QJR", 1, 14),
    "Pro00062": ("LCC", "4EB0", 2, 15),
    "Pro00075": ("FoCut5a", "5AJH", 3, 19),
    "Pro00121": ("HiC", "4OYY", 4, 16),
    "Pro00137": ("PHL7", "7NEI", 5, 14),
}
OMITTED = {"protein_label", "protein_caveat", "protein_short", "figure_candidate"}
RENAMED = {
    "residue": "simulation_residue_number",
    "landmark_class": "assigned_contact_class",
    "manuscript_status": "reporting_status",
}
SUMMARY_STATUSES = {"candidate_with_caveats", "descriptive_core"}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_verified(path: Path, expected: str) -> tuple[list[str], list[dict[str, str]]]:
    raw = path.read_bytes()
    actual = sha256(raw)
    if actual != expected:
        raise ValueError(f"Source digest mismatch for {path}: {actual}")
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8"), newline=""))
    rows = list(reader)
    if reader.fieldnames is None or any(None in row for row in rows):
        raise ValueError(f"Malformed CSV: {path}")
    return reader.fieldnames, rows


def row_key(row: dict[str, str]) -> tuple[str, str]:
    return row["protein_id"], row["residue"]


def table_class(row: dict[str, str], included: bool) -> str:
    if not included:
        return ""
    if row["manuscript_status"] == "candidate_with_caveats":
        if row["shift_class"] != "L20_enriched_supported_candidate":
            raise ValueError("Unexpected supported-shift classification")
        return "supported L20 increase"
    if row["manuscript_status"] == "descriptive_core" and row["stable_contact_core"] == "True":
        return "stable contact core"
    raise ValueError("Unexpected Table S2 membership")


def verify_table(rows: list[dict[str, str]], table_path: Path) -> None:
    expected = []
    for row in rows:
        signed = lambda name: f"{float(row[name]):+.3f}"
        expected.append([
            PROTEINS[row["protein_id"]][0], row["residue_label"], row["landmark_class"],
            f"{row['n_l10']}/{row['n_l20']}",
            f"{float(row['mean_contact_l10']):.3f}", f"{float(row['mean_contact_l20']):.3f}",
            signed("mean_delta_l20_minus_l10"),
            f"[{signed('delta_ci95_low')}, {signed('delta_ci95_high')}]",
            signed("direction_adjusted_delta"), table_class(row, True),
        ])
    enzyme_names = tuple(meta[0] + " & " for meta in PROTEINS.values())
    actual = [
        [cell.strip() for cell in line.removesuffix("\\\\").strip().split(" & ")]
        for line in table_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(enzyme_names)
    ]
    if actual != expected:
        raise ValueError("Table S2 cells do not match the registered 29-row summary")


def build(source_path: Path, summary_path: Path, table_path: Path) -> tuple[bytes, dict]:
    fields, rows = read_verified(source_path, SOURCE_SHA256)
    summary_fields, summary_rows = read_verified(summary_path, SUMMARY_SHA256)
    if fields != summary_fields or len(fields) != 67:
        raise ValueError("The registered full and summary schemas differ")
    if len(rows) != 95 or len(summary_rows) != 29:
        raise ValueError("Expected 95 complete candidates and 29 summary entries")
    full = {row_key(row): row for row in rows}
    subset = {row_key(row): row for row in summary_rows}
    if len(full) != 95 or len(subset) != 29:
        raise ValueError("Duplicate enzyme-residue pairs")
    if Counter(row["protein_id"] for row in rows) != Counter({k: v[3] for k, v in PROTEINS.items()}):
        raise ValueError("Unexpected per-enzyme candidate counts")
    if any(key not in full or row != full[key] for key, row in subset.items()):
        raise ValueError("Summary entries differ from their complete-source rows")
    selected = {key for key, row in full.items() if row["manuscript_status"] in SUMMARY_STATUSES}
    if selected != set(subset):
        raise ValueError("The source reporting-status rule does not reproduce Table S2 membership")
    if any(min(int(row["n_l10"]), int(row["n_l20"])) < 3 for row in summary_rows):
        raise ValueError("Unexpected underpowered Table S2 entry")
    classes = Counter(table_class(row, True) for row in summary_rows)
    if classes != Counter({"supported L20 increase": 8, "stable contact core": 21}):
        raise ValueError("Expected eight supported increases and 21 tabulated stable cores")
    if sum(row["stable_contact_core"] == "True" for row in rows) != 24:
        raise ValueError("Unexpected complete-source stable-core count")
    verify_table(summary_rows, table_path)

    kept_fields = [field for field in fields if field not in OMITTED]
    output_fields = [RENAMED.get(field, field) for field in kept_fields]
    output_fields[1:1] = ["enzyme", "pdb_reference"]
    output_fields += ["included_in_table_s3", "table_s3_classification"]
    output_rows = []
    for row in sorted(rows, key=lambda r: (PROTEINS[r["protein_id"]][2], int(r["residue"]))):
        record = {RENAMED.get(field, field): row[field] for field in kept_fields}
        record["enzyme"], record["pdb_reference"] = PROTEINS[row["protein_id"]][:2]
        included = row_key(row) in subset
        record["included_in_table_s3"] = str(included)
        record["table_s3_classification"] = table_class(row, included)
        # Read and write original strings: do not round or re-estimate source statistics.
        if any(record[RENAMED.get(field, field)] != row[field] for field in kept_fields):
            raise ValueError("A source field changed during projection")
        output_rows.append(record)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=output_fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_rows)
    content = buffer.getvalue().encode("utf-8")
    reread = list(csv.DictReader(io.StringIO(content.decode("utf-8"), newline="")))
    if reread != output_rows:
        raise ValueError("CSV round-trip changed exported values")
    return content, {
        "rows": len(rows), "columns": len(output_fields), "table_s3_rows": len(subset),
        "supported_increases": classes["supported L20 increase"],
        "tabulated_stable_cores": classes["stable contact core"],
        "other_candidates_retained_in_csv": len(rows) - len(subset),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--summary-source", type=Path, default=SUMMARY)
    parser.add_argument("--table-source", type=Path, default=TABLE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true", help="Verify an existing export without writing")
    args = parser.parse_args()
    content, summary = build(args.source, args.summary_source, args.table_source)
    if args.check:
        if not args.output.is_file() or args.output.read_bytes() != content:
            raise SystemExit("The existing Data S2 file differs from the verified export")
    elif args.output.exists():
        if args.output.read_bytes() != content:
            raise SystemExit("Refusing to overwrite a different Data S2 file; review the existing changes first")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(content)
    print(f"source_sha256={SOURCE_SHA256}")
    print(f"summary_source_sha256={SUMMARY_SHA256}")
    print(f"output={args.output}")
    print(f"output_sha256={sha256(content)}")
    for name, value in summary.items():
        print(f"{name}={value}")
    print("status=verified" if args.check else "status=exported_or_already_identical")


if __name__ == "__main__":
    main()
