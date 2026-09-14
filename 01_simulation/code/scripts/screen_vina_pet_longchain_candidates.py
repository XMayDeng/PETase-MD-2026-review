#!/usr/bin/env python3
"""Screen Vina PET-fragment poses for CHARMM long-chain PET placement.

This is a static pre-gate before expensive EM/NVT/NPT/short-MD checks:

  Vina PET_L4 poses x CHARMM PET_L10/L20 windows x forward/reverse directions

The script keeps Vina as a geometry source only. It writes transformed CHARMM
PET complex candidates and a ranked queue for downstream 5 ns retention gates.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

import yaml


SIM_ROOT_REL = "."
DEFAULT_SOURCE_PET_KIND = "PET_L4"
DEFAULT_TARGET_PET_KIND = "PET_L10"
DEFAULT_ACTIVE_SITE_REGISTRY = f"{SIM_ROOT_REL}/inputs/active_site_registry.yaml"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_repo_path(value: str | Path, root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def load_script_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    if spec.loader is None:
        raise ImportError(f"Could not load module from {path}")
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_yaml(path: Path) -> dict:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML is not a mapping: {path}")
    return data


def pet_kind_token(pet_kind: str) -> str:
    return pet_kind.lower()


def default_unit_map(root: Path, pet_kind: str) -> Path:
    token = pet_kind.lower()
    return root / SIM_ROOT_REL / "inputs" / "pets" / pet_kind / "alignment" / f"{token}_unit_map.yaml"


def default_template_pdb(root: Path, pet_kind: str) -> Path:
    token = pet_kind_token(pet_kind)
    return (
        root
        / SIM_ROOT_REL
        / "inputs"
        / "pets"
        / pet_kind
        / "templates"
        / f"{token}-template-001"
        / "input"
        / "pet_final_only_raw.pdb"
    )


def infer_run_config(docking_dir: Path) -> dict:
    path = docking_dir / "run_config.yaml"
    if path.exists():
        return load_yaml(path)
    return {}


def infer_docking_path(docking_dir: Path | None, relpath: str, explicit: Path | None, root: Path) -> Path | None:
    if explicit is not None:
        return resolve_repo_path(explicit, root)
    if docking_dir is not None:
        return docking_dir / relpath
    return None


def parse_window(value: str) -> list[int]:
    if "-" not in value:
        raise argparse.ArgumentTypeError("Window must be formatted like START-END")
    start, end = value.split("-", 1)
    return [int(start), int(end)]


def selected_windows(target_unit_map: dict, policy: str, overrides: list[list[int]] | None) -> list[list[int]]:
    if overrides:
        return overrides
    windows = target_unit_map.get("candidate_windows", {})
    if policy not in windows:
        raise KeyError(f"Window policy {policy!r} not found in target unit map")
    return [list(item) for item in windows[policy]]


def affinity_key(value: object) -> float:
    if value is None or value == "":
        return 999.0
    return float(value)


def bool_value(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in {"1", "true", "yes", "y"}
    return bool(value)


def float_or(value: object, default: float = 999.0) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def int_or(value: object, default: int = 0) -> int:
    if value is None or value == "":
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def qc_v2_sort_rank(row: dict) -> int:
    if bool_value(row.get("catalytic_pose_qc_v2_pass")):
        return 0
    verdict = str(row.get("catalytic_pose_qc_v2_verdict") or "").upper()
    if verdict == "WARN":
        return 1
    if verdict == "FAIL":
        return 2
    return 3


def candidate_sort_key(row: dict) -> tuple:
    eligible_rank = 0 if bool_value(row.get("static_eligible_for_md_gate")) else 1
    expected_backbone_n_count = int_or(row.get("qc_v2_expected_oxyanion_backbone_n_donor_hit_count"))
    nearest_expected_backbone_n = float_or(
        row.get("qc_v2_nearest_expected_oxyanion_backbone_n_distance_angstrom")
    )
    backbone_n_count = int_or(row.get("qc_v2_oxyanion_backbone_n_donor_count"))
    nearest_backbone_n = float_or(row.get("qc_v2_nearest_oxyanion_backbone_n_distance_angstrom"))
    attack_angle = float_or(row.get("qc_v2_ser_attack_angle_deg"))
    return (
        eligible_rank,
        qc_v2_sort_rank(row),
        -expected_backbone_n_count,
        nearest_expected_backbone_n,
        -backbone_n_count,
        nearest_backbone_n,
        abs(attack_angle - 105.0),
        int(row["contacts_lt_1a"]),
        int(row["contacts_lt_2a"]),
        float(row["motif_rmsd_angstrom"]),
        float_or(row.get("ser_carbonyl_min_angstrom")),
        affinity_key(row.get("affinity_kcal_mol")),
        int(row["pose"]),
        int(row["window_start"]),
        str(row["direction"]),
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "rank",
        "static_eligible_for_md_gate",
        "pose",
        "affinity_kcal_mol",
        "candidate",
        "window_start",
        "window_end",
        "direction",
        "motif_policy",
        "motif_rmsd_angstrom",
        "n_alignment_atoms",
        "contacts_lt_1a",
        "contacts_lt_2a",
        "contacts_lt_4a",
        "ser_carbonyl_min_angstrom",
        "catalytic_pose_qc_v2_verdict",
        "catalytic_pose_qc_v2_pass",
        "qc_v2_overall_verdict",
        "qc_v2_expected_carbonyl",
        "qc_v2_evaluated_carbonyl",
        "qc_v2_carbonyl_matches_motif",
        "qc_v2_best_carbonyl",
        "qc_v2_ser_carbonyl_c_angstrom",
        "qc_v2_ser_attack_angle_deg",
        "qc_v2_oxyanion_donor_count",
        "qc_v2_oxyanion_backbone_n_donor_count",
        "qc_v2_expected_oxyanion_backbone_n_donor_hit_count",
        "qc_v2_nearest_oxyanion_donor_label",
        "qc_v2_nearest_oxyanion_donor_distance_angstrom",
        "qc_v2_nearest_oxyanion_backbone_n_label",
        "qc_v2_nearest_oxyanion_backbone_n_distance_angstrom",
        "qc_v2_nearest_expected_oxyanion_backbone_n_label",
        "qc_v2_nearest_expected_oxyanion_backbone_n_distance_angstrom",
        "qc_v2_expected_oxyanion_backbone_n_donor_pass",
        "qc_v2_gate_fail_reasons",
        "static_gate_pass",
        "ser_carbonyl_gate_pass",
        "catalytic_pose_qc_v2_gate_pass",
        "complex_pdb",
        "alignment_report",
        "catalytic_pose_qc_v2_report",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_markdown(path: Path, rows: list[dict], *, queued_rows: list[dict], thresholds: dict) -> None:
    lines = [
        "# Vina PET Long-Chain Candidate Screen",
        "",
        "Vina poses are used only as geometry. Downstream MD must use CHARMM PET parameters.",
        "",
        "## Static Gate",
        "",
        f"- Max motif RMSD: `{thresholds['max_motif_rmsd_angstrom']} A`",
        f"- Max severe clashes `<1 A`: `{thresholds['max_severe_clashes_lt_1a']}`",
        f"- Max Ser-carbonyl distance: `{thresholds['max_ser_carbonyl_angstrom']} A`",
        f"- Require catalytic pose QC v2: `{thresholds['require_catalytic_pose_qc_v2']}`",
        "",
        "## Ranked Candidates",
        "",
        "| rank | queue | pose | affinity | candidate | QC v2 | RMSD | Ser-carbonyl | attack angle | expected N hits | nearest expected N | backbone N donors | nearest backbone N | <1A | <2A |",
        "|---:|:---:|---:|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    queued_ids = {row["rank"] for row in queued_rows}
    for row in rows[:30]:
        affinity = row["affinity_kcal_mol"] if row["affinity_kcal_mol"] is not None else "NA"
        ser = row["ser_carbonyl_min_angstrom"] if row["ser_carbonyl_min_angstrom"] is not None else "NA"
        qc_verdict = row.get("catalytic_pose_qc_v2_verdict") or "NA"
        evaluated = row.get("qc_v2_evaluated_carbonyl") or row.get("qc_v2_best_carbonyl") or "NA"
        attack_angle = row.get("qc_v2_ser_attack_angle_deg") or "NA"
        backbone_count = row.get("qc_v2_oxyanion_backbone_n_donor_count")
        backbone_count = backbone_count if backbone_count is not None else "NA"
        nearest_backbone = row.get("qc_v2_nearest_oxyanion_backbone_n_distance_angstrom") or "NA"
        expected_count = row.get("qc_v2_expected_oxyanion_backbone_n_donor_hit_count")
        expected_count = expected_count if expected_count is not None else "NA"
        nearest_expected = row.get("qc_v2_nearest_expected_oxyanion_backbone_n_distance_angstrom") or "NA"
        lines.append(
            f"| {row['rank']} | {'yes' if row['rank'] in queued_ids else 'no'} | "
            f"{row['pose']} | {affinity} | {row['candidate']} | "
            f"{qc_verdict} `{evaluated}` | {row['motif_rmsd_angstrom']} | {ser} | {attack_angle} | "
            f"{expected_count} | {nearest_expected} | "
            f"{backbone_count} | {nearest_backbone} | "
            f"{row['contacts_lt_1a']} | {row['contacts_lt_2a']} |"
        )
    lines.extend(
        [
            "",
            "## Next Step",
            "",
            "Run EM/NVT/NPT + 5 ns retention gates only for candidates in `short_md_gate_queue/`.",
            "Do not submit 100 ns production MD until a candidate passes the 5 ns retention gate.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def stage_short_md_queue(out_dir: Path, rows: list[dict], stage_top_n: int) -> list[dict]:
    queue_dir = out_dir / "short_md_gate_queue"
    if queue_dir.exists():
        shutil.rmtree(queue_dir)
    queue_dir.mkdir(parents=True, exist_ok=True)
    queued: list[dict] = []
    for row in [item for item in rows if bool_value(item["static_eligible_for_md_gate"])][:stage_top_n]:
        label = f"candidate_{len(queued) + 1:03d}_pose{int(row['pose']):03d}_{row['candidate']}"
        dest = queue_dir / label
        dest.mkdir(parents=True, exist_ok=True)
        shutil.copy2(Path(row["complex_pdb"]), dest / "protein_pet_complex.pdb")
        shutil.copy2(Path(row["alignment_report"]), dest / "alignment_report.yaml")
        qc_report = row.get("catalytic_pose_qc_v2_report")
        if qc_report:
            qc_dir = Path(qc_report).parent
            if qc_dir.exists():
                shutil.copytree(qc_dir, dest / "catalytic_pose_qc_v2", dirs_exist_ok=True)
        meta = dict(row)
        meta["queued_candidate_dir"] = str(dest)
        (dest / "screening_candidate.yaml").write_text(
            yaml.safe_dump(meta, sort_keys=False),
            encoding="utf-8",
        )
        queued.append(meta)
    return queued


def active_site_from_run_config(run_config: dict) -> dict[str, object]:
    active_site = run_config.get("active_site") if isinstance(run_config.get("active_site"), dict) else {}
    triad = active_site.get("triad") if isinstance(active_site.get("triad"), dict) else {}
    ser_atom = run_config.get("ser_atom") if isinstance(run_config.get("ser_atom"), dict) else {}
    oxyanion_hole = active_site.get("oxyanion_hole") if isinstance(active_site.get("oxyanion_hole"), dict) else {}
    return {
        "ser_resid": run_config.get("ser_resid") or triad.get("Ser") or triad.get("ser"),
        "ser_atom_name": ser_atom.get("name") or active_site.get("ser_atom") or "OG",
        "his_resid": triad.get("His") or triad.get("his"),
        "acid_resid": triad.get("Acid") or triad.get("acid") or triad.get("Asp") or triad.get("asp"),
        "expected_oxyanion_backbone_n_donors": oxyanion_hole.get(
            "expected_backbone_n_donors",
            [],
        ),
    }


def qc_v2_thresholds(args: argparse.Namespace) -> dict[str, float | int]:
    return {
        "max_ser_carbonyl_c_distance_angstrom": args.qc_v2_max_ser_carbonyl_c_distance,
        "min_ser_attack_angle_deg": args.qc_v2_min_ser_attack_angle,
        "max_ser_attack_angle_deg": args.qc_v2_max_ser_attack_angle,
        "max_oxyanion_donor_distance_angstrom": args.qc_v2_max_oxyanion_donor_distance,
        "min_oxyanion_donors": args.qc_v2_min_oxyanion_donors,
        "min_oxyanion_backbone_n_donors": args.qc_v2_min_oxyanion_backbone_n_donors,
        "min_expected_oxyanion_backbone_n_donors": args.qc_v2_min_expected_oxyanion_backbone_n_donors,
        "max_ser_his_distance_angstrom": args.qc_v2_max_ser_his_distance,
        "max_severe_clashes_lt_1a": args.max_severe_clashes,
    }


def effective_ser_carbonyl_threshold(args: argparse.Namespace, *, require_qc_v2: bool) -> float:
    configured = float(args.max_ser_carbonyl)
    if not require_qc_v2:
        return configured
    return max(configured, float(args.qc_v2_max_ser_carbonyl_c_distance))


def run_catalytic_pose_qc_v2(
    *,
    qc_mod,
    complex_pdb: Path,
    out_dir: Path,
    ser_resid: int | None,
    ser_atom_name: str,
    his_resid: int | None,
    acid_resid: int | None,
    thresholds: dict[str, float | int],
    expected_oxyanion_backbone_n_donors: list[dict[str, object]] | None = None,
) -> dict | None:
    if ser_resid is None:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    result = qc_mod.score_static_pose(
        atoms=qc_mod.parse_pdb_atoms(complex_pdb),
        ser_resid=int(ser_resid),
        ser_atom_name=str(ser_atom_name or "OG"),
        his_resid=int(his_resid) if his_resid is not None else None,
        acid_resid=int(acid_resid) if acid_resid is not None else None,
        pet_chain="B",
        pet_resnames=set(qc_mod.DEFAULT_PET_RESNAMES),
        thresholds=thresholds,
        expected_oxyanion_backbone_n_donors=expected_oxyanion_backbone_n_donors or [],
    )
    result["inputs"] = {
        "complex_pdb": str(complex_pdb),
        "source": "screen_vina_pet_longchain_candidates",
    }
    (out_dir / "catalytic_pose_qc_v2.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    (out_dir / "catalytic_pose_qc_v2.yaml").write_text(
        yaml.safe_dump(result, sort_keys=False),
        encoding="utf-8",
    )
    qc_mod.write_csv(out_dir / "catalytic_pose_qc_v2.csv", result)
    qc_mod.write_markdown(out_dir / "catalytic_pose_qc_v2.md", result)
    return result


def screen_pose_worker(payload: tuple[dict, dict]) -> list[dict]:
    """Screen one Vina pose against all requested long-chain windows.

    The worker writes only under one pose-specific alignment directory and
    returns row dictionaries for the parent process to rank and summarize.
    """
    pose_spec, common = payload
    root = repo_root()
    script_dir = root / SIM_ROOT_REL / "code" / "scripts"
    align_mod = load_script_module(
        script_dir / "align_charmm_pet_to_vina_pose.py",
        f"align_charmm_pet_to_vina_pose_worker_{os.getpid()}_{pose_spec['index']}",
    )
    qc_mod = load_script_module(
        script_dir / "score_catalytic_pose_qc.py",
        f"score_catalytic_pose_qc_worker_{os.getpid()}_{pose_spec['index']}",
    )

    pose_index = int(pose_spec["index"])
    pose_id = f"pose_{pose_index:03d}"
    pose_pdb = Path(pose_spec["pose_pdb"])
    ligand_atom_map = Path(common["ligand_atom_map"])
    source_unit_map = Path(common["source_unit_map"])
    target_pet_pdb = Path(common["target_pet_pdb"])
    target_unit_map_path = Path(common["target_unit_map_path"])
    protein_pdb = Path(common["protein_pdb"]) if common.get("protein_pdb") else None
    alignment_root = Path(common["alignment_root"])
    align_out = alignment_root / pose_id

    align_args = align_mod.build_arg_parser().parse_args(
        [
            "--vina-pose-pdb",
            str(pose_pdb),
            "--ligand-atom-map",
            str(ligand_atom_map),
            "--source-unit-map",
            str(source_unit_map),
            "--target-pet-pdb",
            str(target_pet_pdb),
            "--target-unit-map",
            str(target_unit_map_path),
            "--motif-policy",
            str(common["motif_policy"]),
            "--max-motif-rmsd",
            str(common["max_motif_rmsd"]),
            "--max-severe-clashes",
            str(common["max_severe_clashes"]),
            "--out-dir",
            str(align_out),
            "--force",
        ]
    )
    if protein_pdb is not None:
        align_args.protein_pdb = protein_pdb
    if common.get("ser_resid") is not None:
        align_args.ser_resid = int(common["ser_resid"])
        align_args.ser_atom_name = str(common["ser_atom_name"])
    align_args.window = common["windows"]

    rows: list[dict] = []
    summary = align_mod.run_alignment(align_args)
    for candidate in summary["candidates"]:
        contacts = candidate["contacts"]
        ser_carbonyl = candidate.get("ser_carbonyl_min_angstrom")
        require_ser_gate = bool(common["require_ser_gate"])
        max_ser_carbonyl = float(common["max_ser_carbonyl"])
        ser_gate_pass = (
            (ser_carbonyl is not None and float(ser_carbonyl) <= max_ser_carbonyl)
            or (ser_carbonyl is None and not require_ser_gate)
        )
        static_gate_pass = bool_value(candidate.get("static_gate_pass"))
        complex_pdb = align_out / candidate["complex_pdb"]
        alignment_report = align_out / "candidates" / candidate["candidate"] / "alignment_report.yaml"
        alignment_report_data = load_yaml(alignment_report)
        expected_carbonyl = expected_target_carbonyl_from_alignment_report(alignment_report_data)
        qc_dir = alignment_report.parent / "catalytic_pose_qc_v2"
        qc_result = run_catalytic_pose_qc_v2(
            qc_mod=qc_mod,
            complex_pdb=complex_pdb,
            out_dir=qc_dir,
            ser_resid=int(common["ser_resid"]) if common.get("ser_resid") is not None else None,
            ser_atom_name=str(common["ser_atom_name"]),
            his_resid=int(common["his_resid"]) if common.get("his_resid") is not None else None,
            acid_resid=int(common["acid_resid"]) if common.get("acid_resid") is not None else None,
            thresholds=common["qc_thresholds"],
            expected_oxyanion_backbone_n_donors=common["expected_oxyanion_backbone_n_donors"],
        )
        qc_report = qc_dir / "catalytic_pose_qc_v2.md" if qc_result else None
        qc_fields = qc_v2_row_fields(qc_result, qc_report, expected_carbonyl=expected_carbonyl)
        require_qc_v2 = bool(common["require_qc_v2"])
        qc_gate_pass = (
            bool_value(qc_fields["catalytic_pose_qc_v2_gate_pass"])
            if require_qc_v2
            else True
        )
        row = {
            "rank": None,
            "static_eligible_for_md_gate": bool(static_gate_pass and ser_gate_pass and qc_gate_pass),
            "pose": pose_index,
            "affinity_kcal_mol": pose_spec["affinity"],
            "candidate": candidate["candidate"],
            "window_start": int(candidate["window"][0]),
            "window_end": int(candidate["window"][1]),
            "direction": candidate["direction"],
            "motif_policy": candidate["motif_policy"],
            "motif_rmsd_angstrom": float(candidate["motif_rmsd_angstrom"]),
            "n_alignment_atoms": int(candidate["n_alignment_atoms"]),
            "contacts_lt_1a": int(contacts["lt_1a"]),
            "contacts_lt_2a": int(contacts["lt_2a"]),
            "contacts_lt_4a": int(contacts["lt_4a"]),
            "ser_carbonyl_min_angstrom": ser_carbonyl,
            "static_gate_pass": static_gate_pass,
            "ser_carbonyl_gate_pass": ser_gate_pass,
            "complex_pdb": str(complex_pdb),
            "alignment_report": str(alignment_report),
        }
        row.update(qc_fields)
        rows.append(row)
    return rows


def carbonyl_label(candidate: dict) -> str:
    return f"{candidate['carbonyl_unit']}:{candidate['carbonyl_atom']}"


def split_carbonyl_label(label: str) -> tuple[int, str]:
    unit, atom = label.split(":", 1)
    return int(unit), atom


def qc_v2_candidate_verdict(candidate: dict) -> str:
    flags = candidate["flags"]
    if bool_value(flags["static_pose_qc_v2_pass"]):
        return "PASS"
    if bool_value(flags["ser_carbonyl_distance_pass"]) and bool_value(flags["clash_pass"]):
        return "WARN"
    return "FAIL"


def expected_target_carbonyl_from_alignment_report(report: dict) -> str | None:
    motif = report.get("motif_metadata") if isinstance(report.get("motif_metadata"), dict) else {}
    source_unit = motif.get("carbonyl_unit")
    source_atom = motif.get("carbonyl_atom")
    if source_unit is None or source_atom is None:
        return None
    source_label = f"{int(source_unit)}:{source_atom}"
    for ref in report.get("alignment_atom_refs", []):
        if "->" not in str(ref):
            continue
        source_ref, target_ref = str(ref).split("->", 1)
        if source_ref == source_label:
            return target_ref
    return None


def find_qc_v2_candidate(result: dict, label: str) -> dict | None:
    unit, atom = split_carbonyl_label(label)
    for candidate in result.get("candidates", []):
        if int(candidate["carbonyl_unit"]) == unit and str(candidate["carbonyl_atom"]) == atom:
            return candidate
    return None


def qc_v2_row_fields(
    result: dict | None,
    report_path: Path | None,
    *,
    expected_carbonyl: str | None = None,
) -> dict:
    if not result:
        return {
            "catalytic_pose_qc_v2_verdict": None,
            "catalytic_pose_qc_v2_pass": None,
            "qc_v2_overall_verdict": None,
            "qc_v2_expected_carbonyl": expected_carbonyl,
            "qc_v2_evaluated_carbonyl": None,
            "qc_v2_carbonyl_matches_motif": None,
            "qc_v2_best_carbonyl": None,
            "qc_v2_ser_carbonyl_c_angstrom": None,
            "qc_v2_ser_attack_angle_deg": None,
            "qc_v2_oxyanion_donor_count": None,
            "qc_v2_oxyanion_backbone_n_donor_count": None,
            "qc_v2_expected_oxyanion_backbone_n_donor_hit_count": None,
            "qc_v2_nearest_oxyanion_donor_label": None,
            "qc_v2_nearest_oxyanion_donor_distance_angstrom": None,
            "qc_v2_nearest_oxyanion_backbone_n_label": None,
            "qc_v2_nearest_oxyanion_backbone_n_distance_angstrom": None,
            "qc_v2_nearest_expected_oxyanion_backbone_n_label": None,
            "qc_v2_nearest_expected_oxyanion_backbone_n_distance_angstrom": None,
            "qc_v2_expected_oxyanion_backbone_n_donor_pass": None,
            "qc_v2_gate_fail_reasons": None,
            "catalytic_pose_qc_v2_gate_pass": None,
            "catalytic_pose_qc_v2_report": None,
        }
    best = result["best_candidate"]
    evaluated = best
    matches_motif = None
    if expected_carbonyl:
        motif_candidate = find_qc_v2_candidate(result, expected_carbonyl)
        matches_motif = motif_candidate is not None
        evaluated = motif_candidate or best
    evaluated_pass = bool_value(evaluated["flags"]["static_pose_qc_v2_pass"])
    gate_pass = evaluated_pass and (matches_motif is not False)
    nearest_donor = evaluated.get("nearest_oxyanion_donor") or {}
    nearest_backbone = evaluated.get("nearest_oxyanion_backbone_n_donor") or {}
    nearest_expected_backbone = evaluated.get("nearest_expected_oxyanion_backbone_n_donor") or {}
    return {
        "catalytic_pose_qc_v2_verdict": qc_v2_candidate_verdict(evaluated),
        "catalytic_pose_qc_v2_pass": evaluated_pass,
        "qc_v2_overall_verdict": result["verdict"],
        "qc_v2_expected_carbonyl": expected_carbonyl,
        "qc_v2_evaluated_carbonyl": carbonyl_label(evaluated),
        "qc_v2_carbonyl_matches_motif": matches_motif,
        "qc_v2_best_carbonyl": carbonyl_label(best),
        "qc_v2_ser_carbonyl_c_angstrom": evaluated["ser_og_to_carbonyl_c_angstrom"],
        "qc_v2_ser_attack_angle_deg": evaluated["ser_og_carbonyl_c_o_angle_deg"],
        "qc_v2_oxyanion_donor_count": evaluated["oxyanion_donor_count"],
        "qc_v2_oxyanion_backbone_n_donor_count": evaluated["oxyanion_backbone_n_donor_count"],
        "qc_v2_expected_oxyanion_backbone_n_donor_hit_count": evaluated.get(
            "expected_oxyanion_backbone_n_donor_hit_count"
        ),
        "qc_v2_nearest_oxyanion_donor_label": nearest_donor.get("label"),
        "qc_v2_nearest_oxyanion_donor_distance_angstrom": nearest_donor.get("distance_angstrom"),
        "qc_v2_nearest_oxyanion_backbone_n_label": nearest_backbone.get("label"),
        "qc_v2_nearest_oxyanion_backbone_n_distance_angstrom": nearest_backbone.get("distance_angstrom"),
        "qc_v2_nearest_expected_oxyanion_backbone_n_label": nearest_expected_backbone.get("label"),
        "qc_v2_nearest_expected_oxyanion_backbone_n_distance_angstrom": nearest_expected_backbone.get("distance_angstrom"),
        "qc_v2_expected_oxyanion_backbone_n_donor_pass": evaluated["flags"].get(
            "expected_oxyanion_backbone_n_donor_pass"
        ),
        "qc_v2_gate_fail_reasons": ";".join(evaluated.get("gate_fail_reasons", [])),
        "catalytic_pose_qc_v2_gate_pass": gate_pass,
        "catalytic_pose_qc_v2_report": str(report_path) if report_path else None,
    }


def screen_candidates(args: argparse.Namespace) -> dict:
    root = repo_root()
    script_dir = root / SIM_ROOT_REL / "code" / "scripts"
    vina_mod = load_script_module(
        script_dir / "run_vina_pet_fragment_docking.py",
        "run_vina_pet_fragment_docking_for_screen",
    )
    align_mod = load_script_module(
        script_dir / "align_charmm_pet_to_vina_pose.py",
        "align_charmm_pet_to_vina_pose_for_screen",
    )
    qc_mod = load_script_module(
        script_dir / "score_catalytic_pose_qc.py",
        "score_catalytic_pose_qc_for_screen",
    )

    docking_dir = resolve_repo_path(args.docking_dir, root) if args.docking_dir else None
    run_config = infer_run_config(docking_dir) if docking_dir else {}

    poses_pdbqt = infer_docking_path(docking_dir, "outputs/vina_poses.pdbqt", args.poses_pdbqt, root)
    ligand_atom_map = infer_docking_path(docking_dir, "pdbqt/ligand_atom_map.tsv", args.ligand_atom_map, root)
    if poses_pdbqt is None or ligand_atom_map is None:
        raise ValueError("--docking-dir or both --poses-pdbqt and --ligand-atom-map are required")

    protein_pdb = (
        resolve_repo_path(args.protein_pdb, root)
        if args.protein_pdb
        else resolve_repo_path(run_config["protein_pdb"], root)
        if run_config.get("protein_pdb")
        else None
    )
    active_site = active_site_from_run_config(run_config)
    if not active_site.get("expected_oxyanion_backbone_n_donors") and run_config.get("protein_id"):
        registry_path = resolve_repo_path(args.active_site_registry, root)
        if registry_path.exists():
            try:
                registry_site = qc_mod.load_active_site_from_registry(
                    registry_path,
                    str(run_config["protein_id"]),
                )
                active_site["expected_oxyanion_backbone_n_donors"] = registry_site.get(
                    "expected_oxyanion_backbone_n_donors",
                    [],
                )
            except (KeyError, ValueError):
                pass
    ser_resid = args.ser_resid or active_site.get("ser_resid")
    ser_atom_name = args.ser_atom_name or active_site.get("ser_atom_name") or "OG"
    his_resid = args.his_resid or active_site.get("his_resid")
    acid_resid = args.acid_resid or active_site.get("acid_resid")
    expected_oxyanion_backbone_n_donors = active_site.get(
        "expected_oxyanion_backbone_n_donors",
        [],
    )

    source_unit_map = (
        resolve_repo_path(args.source_unit_map, root)
        if args.source_unit_map
        else default_unit_map(root, args.source_pet_kind)
    )
    target_unit_map_path = (
        resolve_repo_path(args.target_unit_map, root)
        if args.target_unit_map
        else default_unit_map(root, args.target_pet_kind)
    )
    target_pet_pdb = (
        resolve_repo_path(args.target_pet_pdb, root)
        if args.target_pet_pdb
        else default_template_pdb(root, args.target_pet_kind)
    )
    out_dir = resolve_repo_path(args.out_dir, root)
    if out_dir.exists():
        if not args.force:
            raise FileExistsError(f"Output directory already exists: {out_dir}")
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    target_unit_map = load_yaml(target_unit_map_path)
    windows = selected_windows(target_unit_map, args.window_policy, args.window)

    poses = vina_mod.parse_pdbqt_poses(poses_pdbqt)
    if not poses:
        raise RuntimeError(f"No Vina poses parsed from {poses_pdbqt}")
    poses = sorted(
        poses,
        key=lambda pose: (
            affinity_key(pose.affinity),
            int(pose.index),
        ),
    )[: args.max_poses]

    pose_dir = out_dir / "poses"
    alignment_root = out_dir / "alignments"
    rows: list[dict] = []
    require_ser_gate = bool(args.require_ser_carbonyl_gate or ser_resid is not None)
    require_qc_v2 = bool(args.require_catalytic_pose_qc_v2)
    max_ser_carbonyl = effective_ser_carbonyl_threshold(args, require_qc_v2=require_qc_v2)
    qc_thresholds = qc_v2_thresholds(args)

    pose_specs: list[dict] = []
    pose_dir.mkdir(parents=True, exist_ok=True)
    for pose in poses:
        pose_id = f"pose_{int(pose.index):03d}"
        pose_pdb = pose_dir / f"{pose_id}.pdb"
        vina_mod.write_pose_as_pdb(
            pose,
            pose_pdb,
            atom_map_path=Path("../") / ligand_atom_map.name,
        )
        pose_specs.append(
            {
                "index": int(pose.index),
                "affinity": pose.affinity,
                "pose_pdb": str(pose_pdb),
            }
        )

    common = {
        "alignment_root": str(alignment_root),
        "ligand_atom_map": str(ligand_atom_map),
        "source_unit_map": str(source_unit_map),
        "target_pet_pdb": str(target_pet_pdb),
        "target_unit_map_path": str(target_unit_map_path),
        "protein_pdb": str(protein_pdb) if protein_pdb else None,
        "windows": windows,
        "motif_policy": args.motif_policy,
        "max_motif_rmsd": args.max_motif_rmsd,
        "max_severe_clashes": args.max_severe_clashes,
        "max_ser_carbonyl": max_ser_carbonyl,
        "require_ser_gate": require_ser_gate,
        "require_qc_v2": require_qc_v2,
        "ser_resid": int(ser_resid) if ser_resid is not None else None,
        "ser_atom_name": str(ser_atom_name),
        "his_resid": int(his_resid) if his_resid is not None else None,
        "acid_resid": int(acid_resid) if acid_resid is not None else None,
        "qc_thresholds": qc_thresholds,
        "expected_oxyanion_backbone_n_donors": expected_oxyanion_backbone_n_donors,
    }
    workers = max(1, int(args.workers))
    pose_workers = min(workers, len(pose_specs))
    if pose_workers > 1:
        with ProcessPoolExecutor(max_workers=pose_workers) as executor:
            futures = [
                executor.submit(screen_pose_worker, (pose_spec, common))
                for pose_spec in pose_specs
            ]
            for future in as_completed(futures):
                rows.extend(future.result())
    else:
        for pose_spec in pose_specs:
            rows.extend(screen_pose_worker((pose_spec, common)))

    rows = sorted(rows, key=candidate_sort_key)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank

    queued_rows = stage_short_md_queue(out_dir, rows, args.stage_top_n)
    thresholds = {
        "max_motif_rmsd_angstrom": args.max_motif_rmsd,
        "max_severe_clashes_lt_1a": args.max_severe_clashes,
        "max_ser_carbonyl_angstrom": max_ser_carbonyl,
        "configured_max_ser_carbonyl_angstrom": args.max_ser_carbonyl,
        "require_catalytic_pose_qc_v2": require_qc_v2,
        "catalytic_pose_qc_v2": qc_thresholds,
    }
    summary = {
        "status": "PASS" if queued_rows else "WARN",
        "note": (
            "PASS means at least one static candidate is queued for 5 ns retention gate; "
            "it does not mean production MD is approved."
        ),
        "n_poses_screened": len(poses),
        "n_candidates_screened": len(rows),
        "n_static_eligible_for_md_gate": len(
            [row for row in rows if bool_value(row["static_eligible_for_md_gate"])]
        ),
        "n_queued_for_short_md_gate": len(queued_rows),
        "inputs": {
            "docking_dir": str(docking_dir) if docking_dir else None,
            "poses_pdbqt": str(poses_pdbqt),
            "ligand_atom_map": str(ligand_atom_map),
            "protein_pdb": str(protein_pdb) if protein_pdb else None,
            "source_unit_map": str(source_unit_map),
            "target_pet_pdb": str(target_pet_pdb),
            "target_unit_map": str(target_unit_map_path),
            "active_site": {
                "ser_resid": int(ser_resid) if ser_resid is not None else None,
                "ser_atom_name": str(ser_atom_name) if ser_atom_name is not None else None,
                "his_resid": int(his_resid) if his_resid is not None else None,
                "acid_resid": int(acid_resid) if acid_resid is not None else None,
                "expected_oxyanion_backbone_n_donors": expected_oxyanion_backbone_n_donors,
            },
        },
        "screening": {
            "source_pet_kind": args.source_pet_kind,
            "target_pet_kind": args.target_pet_kind,
            "motif_policy": args.motif_policy,
            "window_policy": args.window_policy,
            "windows": windows,
            "max_poses": args.max_poses,
            "thresholds": thresholds,
            "require_catalytic_pose_qc_v2": require_qc_v2,
        },
        "top_candidates": rows[: min(20, len(rows))],
        "queued_candidates": queued_rows,
    }
    write_csv(out_dir / "screening_summary.csv", rows)
    (out_dir / "screening_summary.yaml").write_text(
        yaml.safe_dump(summary, sort_keys=False),
        encoding="utf-8",
    )
    (out_dir / "screening_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    write_markdown(out_dir / "screening_summary.md", rows, queued_rows=queued_rows, thresholds=thresholds)
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Screen top-N Vina PET_L4 poses against CHARMM PET_L10/L20 placement candidates."
    )
    parser.add_argument("--docking-dir", type=Path)
    parser.add_argument("--poses-pdbqt", type=Path)
    parser.add_argument("--ligand-atom-map", type=Path)
    parser.add_argument("--protein-pdb", type=Path)
    parser.add_argument("--active-site-registry", default=DEFAULT_ACTIVE_SITE_REGISTRY, type=Path)
    parser.add_argument("--ser-resid", type=int)
    parser.add_argument("--ser-atom-name", default=None)
    parser.add_argument("--his-resid", type=int)
    parser.add_argument("--acid-resid", type=int)
    parser.add_argument("--source-pet-kind", default=DEFAULT_SOURCE_PET_KIND)
    parser.add_argument("--target-pet-kind", default=DEFAULT_TARGET_PET_KIND)
    parser.add_argument("--source-unit-map", type=Path)
    parser.add_argument("--target-pet-pdb", type=Path)
    parser.add_argument("--target-unit-map", type=Path)
    parser.add_argument("--motif-policy", default="catalytic_ester_auto")
    parser.add_argument(
        "--window-policy",
        choices=["limited_first_pass", "all_windows"],
        default="all_windows",
    )
    parser.add_argument("--window", action="append", type=parse_window)
    parser.add_argument("--max-poses", type=int, default=10)
    parser.add_argument("--stage-top-n", type=int, default=8)
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("AUTOPROPET_SCREEN_POSE_WORKERS", "1")),
        help="Parallel pose workers within one conformer screen.",
    )
    parser.add_argument("--max-motif-rmsd", type=float, default=3.0)
    parser.add_argument("--max-severe-clashes", type=int, default=0)
    parser.add_argument("--max-ser-carbonyl", type=float, default=4.0)
    parser.add_argument("--qc-v2-max-ser-carbonyl-c-distance", type=float, default=5.0)
    parser.add_argument("--qc-v2-min-ser-attack-angle", type=float, default=80.0)
    parser.add_argument("--qc-v2-max-ser-attack-angle", type=float, default=130.0)
    parser.add_argument("--qc-v2-max-oxyanion-donor-distance", type=float, default=3.5)
    parser.add_argument("--qc-v2-min-oxyanion-donors", type=int, default=1)
    parser.add_argument("--qc-v2-min-oxyanion-backbone-n-donors", type=int, default=1)
    parser.add_argument("--qc-v2-min-expected-oxyanion-backbone-n-donors", type=int, default=1)
    parser.add_argument("--qc-v2-max-ser-his-distance", type=float, default=4.0)
    parser.add_argument(
        "--require-ser-carbonyl-gate",
        action="store_true",
        help="Require Ser-carbonyl gate even if no serine residue was provided.",
    )
    parser.add_argument(
        "--require-catalytic-pose-qc-v2",
        action="store_true",
        help="Require static catalytic pose QC v2 PASS before staging a candidate for short MD.",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = screen_candidates(args)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
