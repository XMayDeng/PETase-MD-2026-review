#!/usr/bin/env python3
"""Run AutoDock Vina PET-fragment docking and geometry QC.

This script is intentionally a geometry probe:
  - AutoDock Vina proposes short PET-fragment poses.
  - Geometry QC checks catalytic distance and severe clashes.
  - Vina PDBQT charges, atom types, and scores are not MD force-field inputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from itertools import zip_longest
from pathlib import Path
from typing import Iterable

import yaml


DEFAULT_GROOVE_CONFIG = (
    "CodeData/ML/exploration/models/b3_interaction_mode/groove_config.json"
)
DEFAULT_ACTIVE_SITE_REGISTRY = (
    "./inputs/active_site_registry.yaml"
)
DEFAULT_DOCKING_PROFILE = "vina-local"
DEFAULT_VINA_ENV_REL = "extra_env/autodock_vina_2026.04/vina"
SIM_ROOT_REL = "."


@dataclass
class Atom:
    serial: int | None
    name: str
    resname: str
    chain: str
    resid: int | None
    x: float
    y: float
    z: float
    element: str
    atom_type: str | None = None
    charge: float | None = None

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_heavy(self) -> bool:
        element = (self.element or "").upper()
        atom_type = (self.atom_type or "").upper()
        name = self.name.strip().upper()
        return not (element == "H" or atom_type.startswith("H") or name.startswith("H"))


@dataclass
class Pose:
    index: int
    affinity: float | None
    atoms: list[Atom]


@dataclass
class PoseMetrics:
    pose: int
    affinity_kcal_mol: float | None
    heavy_atoms: int
    centroid_to_ser_og_angstrom: float
    min_atom_to_ser_og_angstrom: float
    min_carbonyl_c_to_ser_og_angstrom: float | None
    nearest_carbonyl_atom: str | None
    contacts_lt_4a: int
    close_contacts_lt_2a: int
    severe_clashes_lt_1a: int
    catalytic_distance_pass: bool
    clash_pass: bool
    recommended: bool


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_repo_path(value: str | Path, root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def profile_path(profile: str | None, root: Path) -> Path | None:
    if not profile:
        return None
    path = Path(profile)
    if path.suffix or "/" in profile:
        return resolve_repo_path(path, root)
    return root / SIM_ROOT_REL / "code" / "profiles" / "docking" / f"{profile}.yaml"


def load_profile(profile: str | None, root: Path) -> dict:
    path = profile_path(profile, root)
    if path is None:
        return {}
    if not path.exists():
        raise FileNotFoundError(f"Docking profile not found: {path}")
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Docking profile is not a mapping: {path}")
    data["_profile_path"] = str(path)
    return data


def resolve_env_path(profile: dict, root: Path) -> Path | None:
    value = os.environ.get("AUTOPROPET_VINA_ENV") or profile.get("vina_env")
    return resolve_repo_path(value, root) if value else None


def tool_command(tool_path: Path, *, env_path: Path | None = None) -> list[str]:
    """Return an executable command for a resolved Vina helper.

    Conda/venv entrypoint wrappers in `extra_env/autodock_vina_2026.04/vina/bin`
    may contain an absolute shebang from the machine where the environment was
    created. Such wrappers may point to an unavailable machine, so run
    Meeko entrypoints through the environment's own Python interpreter.
    """
    if env_path is not None and tool_path.name in {
        "mk_prepare_ligand.py",
        "mk_prepare_receptor.py",
    }:
        env_python = env_path / "bin" / "python"
        if env_python.exists() and os.access(env_python, os.X_OK):
            return [str(env_python), str(tool_path)]
    return [str(tool_path)]


def configure_babel_libdir_from_obabel(obabel_path: Path) -> None:
    if os.environ.get("BABEL_LIBDIR"):
        return
    env_dir = obabel_path.parent.parent
    libdir_parent = env_dir / "lib" / "openbabel"
    if not libdir_parent.is_dir():
        return
    versions = sorted(
        [p for p in libdir_parent.iterdir() if p.is_dir()],
        reverse=True,
    )
    if versions:
        os.environ["BABEL_LIBDIR"] = str(versions[0])


def resolve_tool(
    *,
    tool_name: str,
    explicit: str | None,
    profile: dict,
    profile_key: str,
    env_var: str,
    root: Path,
    required: bool = True,
) -> Path | None:
    candidates: list[Path] = []
    for value in (explicit, os.environ.get(env_var), profile.get(profile_key)):
        if value:
            candidates.append(resolve_repo_path(value, root))

    env_path = resolve_env_path(profile, root)
    if env_path is not None:
        candidates.append(env_path / "bin" / tool_name)

    found = shutil.which(tool_name)
    if found:
        candidates.append(Path(found))

    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate

    if required:
        tried = ", ".join(str(p) for p in candidates) or "PATH"
        raise FileNotFoundError(f"Could not find executable {tool_name}; tried {tried}")
    return None


def parse_pdb_atom_line(line: str) -> Atom | None:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        serial = int(line[6:11])
    except ValueError:
        serial = None
    name = line[12:16].strip()
    resname = line[17:20].strip()
    chain = line[21].strip() or ""
    try:
        resid = int(line[22:26])
    except ValueError:
        resid = None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    element = line[76:78].strip() if len(line) >= 78 else ""
    if not element:
        element = re.sub(r"[^A-Za-z]", "", name)[:1].upper()
    return Atom(serial, name, resname, chain, resid, x, y, z, element)


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        atom = parse_pdb_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records parsed from {path}")
    return atoms


def parse_pdbqt_atom_line(line: str) -> Atom | None:
    atom = parse_pdb_atom_line(line)
    if atom is None:
        return None
    parts = line.split()
    charge = None
    atom_type = None
    if len(parts) >= 2:
        atom_type = parts[-1]
    if len(parts) >= 3:
        try:
            charge = float(parts[-2])
        except ValueError:
            charge = None
    element = atom_type or atom.element
    return Atom(
        atom.serial,
        atom.name,
        atom.resname,
        atom.chain,
        atom.resid,
        atom.x,
        atom.y,
        atom.z,
        element,
        atom_type=atom_type,
        charge=charge,
    )


def parse_pdbqt_poses(path: Path) -> list[Pose]:
    poses: list[Pose] = []
    atoms: list[Atom] = []
    current_index: int | None = None
    current_affinity: float | None = None
    next_index = 1

    def flush() -> None:
        nonlocal atoms, current_index, current_affinity, next_index
        if atoms:
            index = current_index if current_index is not None else next_index
            poses.append(Pose(index, current_affinity, atoms))
            next_index = max(next_index, index + 1)
        atoms = []
        current_index = None
        current_affinity = None

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MODEL"):
            flush()
            parts = line.split()
            current_index = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else next_index
            continue
        if line.startswith("ENDMDL"):
            flush()
            continue
        if line.startswith("REMARK VINA RESULT"):
            match = re.search(r"RESULT:\s+([-+]?\d+(?:\.\d+)?)", line)
            if match:
                current_affinity = float(match.group(1))
            continue
        atom = parse_pdbqt_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    flush()
    if not poses:
        raise ValueError(f"No docking poses parsed from {path}")
    return poses


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def centroid(atoms: Iterable[Atom]) -> tuple[float, float, float]:
    items = list(atoms)
    if not items:
        raise ValueError("Cannot calculate centroid of empty atom list")
    return (
        sum(a.x for a in items) / len(items),
        sum(a.y for a in items) / len(items),
        sum(a.z for a in items) / len(items),
    )


def find_ser_atom(
    protein_atoms: list[Atom],
    ser_resid: int,
    ser_atom_name: str,
    ser_chain: str | None,
) -> Atom:
    candidates = []
    for atom in protein_atoms:
        if atom.resid != ser_resid:
            continue
        if atom.name.strip().upper() != ser_atom_name.upper():
            continue
        if ser_chain and atom.chain != ser_chain:
            continue
        candidates.append(atom)
    if candidates:
        return candidates[0]
    chain_note = f" chain {ser_chain}" if ser_chain else ""
    raise ValueError(f"Cannot find Ser{ser_resid}{chain_note} atom {ser_atom_name}")


def protein_id_lookup_keys(protein_id: str) -> list[str]:
    keys = [str(protein_id)]
    try:
        numeric = int(str(protein_id))
    except ValueError:
        numeric = None
    if numeric is not None:
        keys.extend([f"{numeric:05d}", str(numeric)])
    unique: list[str] = []
    for key in keys:
        if key not in unique:
            unique.append(key)
    return unique


def get_protein_entry(mapping: dict, protein_id: str) -> dict | None:
    for key in protein_id_lookup_keys(protein_id):
        entry = mapping.get(key)
        if isinstance(entry, dict):
            return entry
    return None


def load_ser_resid_from_groove_config(path: Path, protein_id: str) -> int:
    data = json.loads(path.read_text(encoding="utf-8"))
    proteins = data.get("proteins", {})
    entry = get_protein_entry(proteins, protein_id)
    if not entry:
        raise KeyError(f"Protein {protein_id} not found in groove config: {path}")
    triad = entry.get("catalytic_triad", {})
    ser = triad.get("Ser") or triad.get("SER") or triad.get("ser")
    if ser is None:
        raise KeyError(f"Protein {protein_id} has no catalytic Ser in {path}")
    return int(ser)


def load_triad_from_groove_config(path: Path, protein_id: str) -> dict[str, int]:
    data = json.loads(path.read_text(encoding="utf-8"))
    proteins = data.get("proteins", {})
    entry = get_protein_entry(proteins, protein_id)
    if not entry:
        raise KeyError(f"Protein {protein_id} not found in groove config: {path}")
    triad = entry.get("catalytic_triad", {})
    ser = triad.get("Ser") or triad.get("SER") or triad.get("ser")
    asp = triad.get("Asp") or triad.get("ASP") or triad.get("asp") or triad.get("Glu") or triad.get("GLU") or triad.get("glu")
    his = triad.get("His") or triad.get("HIS") or triad.get("his")
    if ser is None:
        raise KeyError(f"Protein {protein_id} has no catalytic Ser in {path}")
    if asp is None or his is None:
        raise KeyError(f"Protein {protein_id} has incomplete catalytic triad in {path}")
    return {"Ser": int(ser), "Acid": int(asp), "His": int(his)}


def load_active_site_from_registry(path: Path, protein_id: str) -> dict[str, object]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Active-site registry is not a mapping: {path}")
    proteins = data.get("proteins", {})
    if not isinstance(proteins, dict):
        raise ValueError(f"Active-site registry has no proteins mapping: {path}")
    entry = get_protein_entry(proteins, protein_id)
    if not entry:
        raise KeyError(f"Protein {protein_id} not found in active-site registry: {path}")

    triad = entry.get("catalytic_triad", {})
    if not isinstance(triad, dict):
        raise ValueError(f"Protein {protein_id} catalytic_triad is not a mapping in {path}")
    ser = triad.get("ser") or triad.get("Ser") or triad.get("SER")
    acid = (
        triad.get("acid")
        or triad.get("Acid")
        or triad.get("asp")
        or triad.get("Asp")
        or triad.get("ASP")
        or triad.get("glu")
        or triad.get("Glu")
        or triad.get("GLU")
    )
    his = triad.get("his") or triad.get("His") or triad.get("HIS")
    if ser is None or acid is None or his is None:
        raise KeyError(f"Protein {protein_id} has incomplete catalytic triad in {path}")

    return {
        "source": "active_site_registry",
        "registry_path": str(path),
        "protein_id": protein_id,
        "status": str(entry.get("status", "candidate_needs_review")),
        "triad": {"Ser": int(ser), "Acid": int(acid), "His": int(his)},
        "ser_atom": entry.get("ser_atom", "OG"),
        "acid_resname": triad.get("acid_resname"),
        "triad_spatial_qc": entry.get("triad_spatial_qc"),
        "oxyanion_hole": entry.get("oxyanion_hole"),
        "evidence": entry.get("evidence", []),
        "notes": entry.get("notes"),
    }


def catalytic_site_qc(protein_atoms: list[Atom], triad: dict[str, int]) -> dict[str, object]:
    ser = find_ser_atom(protein_atoms, triad["Ser"], "OG", None)
    his_atoms = [
        atom
        for atom in protein_atoms
        if atom.resid == triad["His"] and atom.resname.upper() == "HIS" and atom.name.upper() in {"ND1", "NE2"}
    ]
    acid_atoms = [
        atom
        for atom in protein_atoms
        if atom.resid == triad["Acid"]
        and atom.resname.upper() in {"ASP", "GLU"}
        and atom.name.upper() in {"OD1", "OD2", "OE1", "OE2"}
    ]
    if not his_atoms or not acid_atoms:
        return {
            "pass": False,
            "triad": triad,
            "reason": "Missing His or acidic side-chain atoms for catalytic triad QC.",
        }
    ser_his = min(distance(ser.xyz, atom.xyz) for atom in his_atoms)
    his_acid = min(distance(his.xyz, acid.xyz) for his in his_atoms for acid in acid_atoms)
    ser_acid = min(distance(ser.xyz, acid.xyz) for acid in acid_atoms)
    passed = ser_his <= 6.0 and his_acid <= 6.0 and ser_acid <= 12.0
    return {
        "pass": passed,
        "triad": triad,
        "ser_his_min_angstrom": round(ser_his, 3),
        "his_acid_min_angstrom": round(his_acid, 3),
        "ser_acid_min_angstrom": round(ser_acid, 3),
        "thresholds_angstrom": {
            "ser_his_max": 6.0,
            "his_acid_max": 6.0,
            "ser_acid_max": 12.0,
        },
    }


def carbonyl_c_candidates(atoms: list[Atom]) -> list[Atom]:
    candidates = [
        atom
        for atom in atoms
        if atom.is_heavy
        and (atom.atom_type or "").upper() == "C"
        and atom.charge is not None
        and atom.charge >= 0.30
    ]
    if candidates:
        return candidates
    return [
        atom
        for atom in atoms
        if atom.is_heavy
        and atom.name.strip().upper().startswith("C")
        and (atom.atom_type or "").upper() == "C"
    ]


def pair_counts(
    ligand_atoms: list[Atom],
    receptor_atoms: list[Atom],
    cutoffs: tuple[float, ...] = (1.0, 2.0, 4.0),
) -> dict[float, int]:
    counts = {cutoff: 0 for cutoff in cutoffs}
    for lig_atom in ligand_atoms:
        if not lig_atom.is_heavy:
            continue
        for rec_atom in receptor_atoms:
            if not rec_atom.is_heavy:
                continue
            d = distance(lig_atom.xyz, rec_atom.xyz)
            for cutoff in cutoffs:
                if d < cutoff:
                    counts[cutoff] += 1
    return counts


def summarize_poses(
    poses: list[Pose],
    protein_atoms: list[Atom],
    ser_atom: Atom,
    catalytic_cutoff: float,
    max_severe_clashes: int,
) -> list[PoseMetrics]:
    metrics: list[PoseMetrics] = []
    for pose in poses:
        heavy_atoms = [atom for atom in pose.atoms if atom.is_heavy]
        if not heavy_atoms:
            continue
        cent = centroid(heavy_atoms)
        min_atom_ser = min(distance(atom.xyz, ser_atom.xyz) for atom in heavy_atoms)
        carbonyl_atoms = carbonyl_c_candidates(heavy_atoms)
        nearest_carbonyl = None
        min_carbonyl_ser = None
        if carbonyl_atoms:
            nearest_carbonyl = min(carbonyl_atoms, key=lambda atom: distance(atom.xyz, ser_atom.xyz))
            min_carbonyl_ser = distance(nearest_carbonyl.xyz, ser_atom.xyz)
        counts = pair_counts(heavy_atoms, protein_atoms)
        catalytic_pass = (
            min_carbonyl_ser is not None and min_carbonyl_ser <= catalytic_cutoff
        )
        clash_pass = counts[1.0] <= max_severe_clashes
        metrics.append(
            PoseMetrics(
                pose=pose.index,
                affinity_kcal_mol=pose.affinity,
                heavy_atoms=len(heavy_atoms),
                centroid_to_ser_og_angstrom=distance(cent, ser_atom.xyz),
                min_atom_to_ser_og_angstrom=min_atom_ser,
                min_carbonyl_c_to_ser_og_angstrom=min_carbonyl_ser,
                nearest_carbonyl_atom=nearest_carbonyl.name if nearest_carbonyl else None,
                contacts_lt_4a=counts[4.0],
                close_contacts_lt_2a=counts[2.0],
                severe_clashes_lt_1a=counts[1.0],
                catalytic_distance_pass=catalytic_pass,
                clash_pass=clash_pass,
                recommended=False,
            )
        )

    passing = [m for m in metrics if m.catalytic_distance_pass and m.clash_pass]
    if passing:
        # Vina scores are not absolute PET binding energies, but they are useful
        # for ranking poses from the same receptor/ligand/grid after geometry gating.
        recommended = min(
            passing,
            key=lambda m: (
                m.affinity_kcal_mol if m.affinity_kcal_mol is not None else float("inf"),
                m.min_carbonyl_c_to_ser_og_angstrom or float("inf"),
                m.close_contacts_lt_2a,
            ),
        )
        recommended.recommended = True
    return metrics


def metrics_to_row(metric: PoseMetrics) -> dict[str, object]:
    return {
        "pose": metric.pose,
        "affinity_kcal_mol": metric.affinity_kcal_mol,
        "heavy_atoms": metric.heavy_atoms,
        "centroid_to_ser_og_angstrom": f"{metric.centroid_to_ser_og_angstrom:.3f}",
        "min_atom_to_ser_og_angstrom": f"{metric.min_atom_to_ser_og_angstrom:.3f}",
        "min_carbonyl_c_to_ser_og_angstrom": (
            f"{metric.min_carbonyl_c_to_ser_og_angstrom:.3f}"
            if metric.min_carbonyl_c_to_ser_og_angstrom is not None
            else ""
        ),
        "nearest_carbonyl_atom": metric.nearest_carbonyl_atom or "",
        "contacts_lt_4a": metric.contacts_lt_4a,
        "close_contacts_lt_2a": metric.close_contacts_lt_2a,
        "severe_clashes_lt_1a": metric.severe_clashes_lt_1a,
        "catalytic_distance_pass": metric.catalytic_distance_pass,
        "clash_pass": metric.clash_pass,
        "recommended": metric.recommended,
    }


def write_metrics(metrics: list[PoseMetrics], out_dir: Path) -> None:
    csv_path = out_dir / "pose_geometry_summary.csv"
    fieldnames = list(metrics_to_row(metrics[0]).keys()) if metrics else [
        "pose",
        "affinity_kcal_mol",
        "heavy_atoms",
        "centroid_to_ser_og_angstrom",
        "min_atom_to_ser_og_angstrom",
        "min_carbonyl_c_to_ser_og_angstrom",
        "nearest_carbonyl_atom",
        "contacts_lt_4a",
        "close_contacts_lt_2a",
        "severe_clashes_lt_1a",
        "catalytic_distance_pass",
        "clash_pass",
        "recommended",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for metric in metrics:
            writer.writerow(metrics_to_row(metric))

    md_lines = [
        "# Vina PET Fragment Geometry Summary",
        "",
        "Vina scores are used only for same-setup ranking. MD initialization must use CHARMM PET parameters.",
        "",
        "| pose | affinity | SerOG-carbonyl C | severe clashes | contacts <4A | recommended |",
        "|---:|---:|---:|---:|---:|:---:|",
    ]
    for metric in metrics:
        carbonyl = (
            f"{metric.min_carbonyl_c_to_ser_og_angstrom:.3f}"
            if metric.min_carbonyl_c_to_ser_og_angstrom is not None
            else "NA"
        )
        affinity = (
            f"{metric.affinity_kcal_mol:.3f}"
            if metric.affinity_kcal_mol is not None
            else "NA"
        )
        md_lines.append(
            f"| {metric.pose} | {affinity} | {carbonyl} | "
            f"{metric.severe_clashes_lt_1a} | {metric.contacts_lt_4a} | "
            f"{'yes' if metric.recommended else 'no'} |"
        )
    (out_dir / "pose_geometry_summary.md").write_text("\n".join(md_lines) + "\n", encoding="utf-8")


def pdb_atom_record(atom: Atom, serial: int, *, resname: str = "UNL", chain: str = "B", resid: int = 1) -> str:
    name = (atom.name or atom.element or "X")[:4]
    element = (atom.element or name[:1] or "X")[:2].upper()
    return (
        f"HETATM{serial:5d} {name:<4s} {resname:>3s} {chain:1s}{resid:4d}    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def write_pose_as_pdb(
    pose: Pose,
    out_path: Path,
    *,
    atom_map_path: Path | None = None,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "REMARK Generated from AutoDock Vina PDBQT pose coordinates.\n",
        "REMARK Coordinates only; use CHARMM PET assets for MD topology and charges.\n",
        f"REMARK VINA_POSE_INDEX {pose.index}\n",
    ]
    if pose.affinity is not None:
        lines.append(f"REMARK VINA_AFFINITY_KCAL_MOL {pose.affinity:.3f}\n")
    if atom_map_path is not None:
        lines.append(f"REMARK ATOM_IDENTITY_MAP {atom_map_path.as_posix()}\n")
    for serial, atom in enumerate(pose.atoms, start=1):
        lines.append(pdb_atom_record(atom, serial))
    lines.extend(["TER\n", "END\n"])
    out_path.write_text("".join(lines), encoding="utf-8")


def write_ligand_atom_map(
    *,
    charmm_atoms: list[Atom],
    pdbqt_atoms: list[Atom],
    out_tsv: Path,
    summary_yaml: Path,
) -> dict[str, object]:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    charmm_by_serial = {
        atom.serial: atom
        for atom in charmm_atoms
        if atom.serial is not None
    }

    def original_serial_from_pdbqt_atom(atom: Atom) -> int | None:
        match = re.search(r"(\d+)$", atom.name.strip())
        if not match:
            return None
        serial = int(match.group(1))
        return serial if serial in charmm_by_serial else None

    fieldnames = [
        "map_index",
        "charmm_serial",
        "charmm_resid",
        "charmm_resname",
        "charmm_atom",
        "charmm_element",
        "pdbqt_serial",
        "pdbqt_resid",
        "pdbqt_atom",
        "vina_atom_type",
        "pdbqt_charge",
        "mapping_method",
    ]
    rows = []
    mapped_from_pdbqt_names = 0
    for index, pdbqt_atom in enumerate(pdbqt_atoms, start=1):
        original_serial = original_serial_from_pdbqt_atom(pdbqt_atom)
        charmm_atom = charmm_by_serial.get(original_serial) if original_serial is not None else None
        if charmm_atom is not None:
            mapped_from_pdbqt_names += 1
        rows.append(
            {
                "map_index": index,
                "charmm_serial": charmm_atom.serial if charmm_atom else "",
                "charmm_resid": charmm_atom.resid if charmm_atom else "",
                "charmm_resname": charmm_atom.resname if charmm_atom else "",
                "charmm_atom": charmm_atom.name if charmm_atom else "",
                "charmm_element": charmm_atom.element if charmm_atom else "",
                "pdbqt_serial": pdbqt_atom.serial if pdbqt_atom else "",
                "pdbqt_resid": pdbqt_atom.resid if pdbqt_atom else "",
                "pdbqt_atom": pdbqt_atom.name if pdbqt_atom else "",
                "vina_atom_type": pdbqt_atom.atom_type if pdbqt_atom else "",
                "pdbqt_charge": pdbqt_atom.charge if pdbqt_atom and pdbqt_atom.charge is not None else "",
                "mapping_method": "pdbqt_atom_name_original_serial" if charmm_atom else "unmapped",
            }
        )
    if mapped_from_pdbqt_names == 0:
        rows = []
        for index, (charmm_atom, pdbqt_atom) in enumerate(
            zip_longest(charmm_atoms, pdbqt_atoms),
            start=1,
        ):
            rows.append(
                {
                    "map_index": index,
                    "charmm_serial": charmm_atom.serial if charmm_atom else "",
                    "charmm_resid": charmm_atom.resid if charmm_atom else "",
                    "charmm_resname": charmm_atom.resname if charmm_atom else "",
                    "charmm_atom": charmm_atom.name if charmm_atom else "",
                    "charmm_element": charmm_atom.element if charmm_atom else "",
                    "pdbqt_serial": pdbqt_atom.serial if pdbqt_atom else "",
                    "pdbqt_resid": pdbqt_atom.resid if pdbqt_atom else "",
                    "pdbqt_atom": pdbqt_atom.name if pdbqt_atom else "",
                    "vina_atom_type": pdbqt_atom.atom_type if pdbqt_atom else "",
                    "pdbqt_charge": pdbqt_atom.charge if pdbqt_atom and pdbqt_atom.charge is not None else "",
                    "mapping_method": "serial_order_fallback",
                }
            )
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    mapped_rows = [row for row in rows if row["charmm_serial"]]
    mapped_charmm_heavy = {
        int(row["charmm_serial"])
        for row in mapped_rows
        if charmm_by_serial.get(int(row["charmm_serial"]))
        and charmm_by_serial[int(row["charmm_serial"])].is_heavy
    }
    mapping_method = (
        "pdbqt_atom_name_original_serial"
        if mapped_from_pdbqt_names
        else "serial_order_fallback"
    )
    summary = {
        "mapping_method": mapping_method,
        "mapping_complete": len(mapped_rows) == len(pdbqt_atoms),
        "charmm_atom_count": len(charmm_atoms),
        "charmm_heavy_atom_count": sum(atom.is_heavy for atom in charmm_atoms),
        "pdbqt_atom_count": len(pdbqt_atoms),
        "pdbqt_heavy_atom_count": sum(atom.is_heavy for atom in pdbqt_atoms),
        "mapped_pdbqt_atom_count": len(mapped_rows),
        "mapped_charmm_heavy_atom_count": len(mapped_charmm_heavy),
        "note": (
            "This map preserves CHARMM input atom identity against the docking-only "
            "PDBQT atom order. Vina charges/types are not MD inputs."
        ),
    }
    summary_yaml.write_text(yaml.safe_dump(summary, sort_keys=False), encoding="utf-8")
    return summary


def run_command(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(command, capture_output=True, text=True)
    log_path.write_text(
        "$ " + " ".join(command) + "\n\n"
        + "## stdout\n"
        + result.stdout
        + "\n## stderr\n"
        + result.stderr
        + f"\n## returncode\n{result.returncode}\n",
        encoding="utf-8",
    )
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with return code {result.returncode}: {' '.join(command)}")


def write_run_config(
    path: Path,
    *,
    args: argparse.Namespace,
    profile: dict,
    tools: dict[str, Path],
    ser_atom: Atom,
    box_center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    active_site: dict[str, object] | None = None,
    catalytic_qc: dict[str, object] | None = None,
) -> None:
    data = {
        "protein_id": args.protein_id,
        "protein_pdb": str(args.protein_pdb),
        "ligand_pdb": str(args.ligand_pdb),
        "profile": args.profile,
        "profile_path": profile.get("_profile_path"),
        "tools": {key: str(value) for key, value in tools.items()},
        "ser_resid": args.ser_resid,
        "ser_chain": args.ser_chain,
        "ser_atom": {
            "name": ser_atom.name,
            "resid": ser_atom.resid,
            "chain": ser_atom.chain,
            "xyz_angstrom": [ser_atom.x, ser_atom.y, ser_atom.z],
        },
        "box_center_angstrom": list(box_center),
        "box_size_angstrom": list(box_size),
        "exhaustiveness": args.exhaustiveness,
        "num_modes": args.num_modes,
        "seed": args.seed,
        "cpu": args.cpu,
        "note": (
            "AutoDock Vina is used as a PET fragment geometry generator only; "
            "CHARMM PET parameters remain authoritative for MD."
        ),
    }
    if catalytic_qc is not None:
        data["catalytic_site_qc"] = catalytic_qc
    if active_site is not None:
        data["active_site"] = active_site
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def parse_box(values: list[float] | None, fallback: tuple[float, float, float]) -> tuple[float, float, float]:
    if not values:
        return fallback
    if len(values) != 3:
        raise ValueError("Box values must contain exactly 3 floats")
    return (float(values[0]), float(values[1]), float(values[2]))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run AutoDock Vina PET_L4 short-fragment docking and geometry QC."
    )
    parser.add_argument("--protein-id", help="Protein id for active-site lookup, e.g. 00059.")
    parser.add_argument("--protein-pdb", required=True, type=Path)
    parser.add_argument("--ligand-pdb", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--profile", default=DEFAULT_DOCKING_PROFILE)
    parser.add_argument("--active-site-registry", default=DEFAULT_ACTIVE_SITE_REGISTRY, type=Path)
    parser.add_argument(
        "--site-source",
        choices=["auto", "registry", "groove-config", "manual"],
        default="auto",
        help="Active-site source. auto uses the registry first, then groove_config as fallback.",
    )
    parser.add_argument("--groove-config", default=DEFAULT_GROOVE_CONFIG, type=Path)
    parser.add_argument("--ser-resid", type=int)
    parser.add_argument("--ser-chain")
    parser.add_argument("--ser-atom-name", default="OG")
    parser.add_argument("--box-center", nargs=3, type=float)
    parser.add_argument("--box-size", nargs=3, type=float)
    parser.add_argument("--exhaustiveness", type=int)
    parser.add_argument("--num-modes", type=int)
    parser.add_argument("--seed", type=int, default=20260413)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--catalytic-cutoff", type=float, default=4.0)
    parser.add_argument("--max-severe-clashes", type=int, default=0)
    parser.add_argument(
        "--allow-invalid-catalytic-site",
        action="store_true",
        help="Continue even if the selected catalytic triad is invalid or not spatially clustered.",
    )
    parser.add_argument(
        "--require-trusted-active-site",
        action="store_true",
        help="Reject registry entries whose status is not trusted.",
    )
    parser.add_argument("--vina-bin")
    parser.add_argument("--obabel-bin")
    parser.add_argument("--prepare-ligand-bin")
    parser.add_argument("--prepare-receptor-bin")
    parser.add_argument("--force", action="store_true", help="Replace an existing output directory.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    root = repo_root()
    profile = load_profile(args.profile, root)

    args.protein_pdb = resolve_repo_path(args.protein_pdb, root)
    args.ligand_pdb = resolve_repo_path(args.ligand_pdb, root)
    args.active_site_registry = resolve_repo_path(args.active_site_registry, root)
    args.groove_config = resolve_repo_path(args.groove_config, root)
    args.out_dir = resolve_repo_path(args.out_dir, root)

    if args.out_dir.exists():
        if not args.force:
            raise FileExistsError(f"Output directory already exists: {args.out_dir}")
        shutil.rmtree(args.out_dir)
    args.out_dir.mkdir(parents=True, exist_ok=False)
    inputs_dir = args.out_dir / "inputs"
    pdbqt_dir = args.out_dir / "pdbqt"
    outputs_dir = args.out_dir / "outputs"
    logs_dir = args.out_dir / "logs"
    for path in (inputs_dir, pdbqt_dir, outputs_dir, logs_dir):
        path.mkdir(parents=True, exist_ok=True)

    protein_copy = inputs_dir / args.protein_pdb.name
    ligand_copy = inputs_dir / args.ligand_pdb.name
    shutil.copy2(args.protein_pdb, protein_copy)
    shutil.copy2(args.ligand_pdb, ligand_copy)

    protein_atoms = parse_pdb_atoms(protein_copy)
    ser_resid = args.ser_resid
    active_site: dict[str, object] | None = None
    triad_for_qc = None
    if args.site_source == "manual" and ser_resid is None:
        raise ValueError("--site-source manual requires --ser-resid")
    if ser_resid is None:
        if not args.protein_id:
            raise ValueError("--protein-id or --ser-resid is required")
        if args.site_source in {"auto", "registry"}:
            if args.active_site_registry.exists():
                try:
                    active_site = load_active_site_from_registry(
                        args.active_site_registry,
                        args.protein_id,
                    )
                    triad_for_qc = active_site["triad"]  # type: ignore[assignment]
                    ser_resid = triad_for_qc["Ser"]
                    args.ser_resid = ser_resid
                    registry_ser_atom = active_site.get("ser_atom")
                    if registry_ser_atom and args.ser_atom_name == "OG":
                        args.ser_atom_name = str(registry_ser_atom)
                except KeyError:
                    if args.site_source == "registry":
                        raise
            elif args.site_source == "registry":
                raise FileNotFoundError(
                    f"Active-site registry not found: {args.active_site_registry}"
                )
        if ser_resid is None and args.site_source in {"auto", "groove-config"}:
            triad_for_qc = load_triad_from_groove_config(args.groove_config, args.protein_id)
            ser_resid = triad_for_qc["Ser"]
            args.ser_resid = ser_resid
            active_site = {
                "source": "groove_config",
                "groove_config_path": str(args.groove_config),
                "protein_id": args.protein_id,
                "status": "candidate_unverified",
                "triad": triad_for_qc,
            }
        if ser_resid is None:
            raise KeyError(
                f"Protein {args.protein_id} has no active-site entry in the selected source"
            )
    else:
        active_site = {
            "source": "manual_ser_resid",
            "protein_id": args.protein_id,
            "status": "manual_unverified",
            "ser_resid": ser_resid,
            "ser_atom": args.ser_atom_name,
        }
    ser_atom = find_ser_atom(protein_atoms, ser_resid, args.ser_atom_name, args.ser_chain)
    catalytic_qc = catalytic_site_qc(protein_atoms, triad_for_qc) if triad_for_qc else None
    if active_site is not None and catalytic_qc is not None:
        active_site["runtime_spatial_qc"] = catalytic_qc
    site_status = str(active_site.get("status")) if active_site else ""
    invalid_site_status = site_status == "invalid"
    failed_spatial_qc = catalytic_qc is not None and not catalytic_qc.get("pass")
    untrusted_required = (
        args.require_trusted_active_site
        and active_site is not None
        and site_status != "trusted"
    )
    if (invalid_site_status or failed_spatial_qc or untrusted_required) and not args.allow_invalid_catalytic_site:
        report = {
            "active_site": active_site,
            "catalytic_site_qc": catalytic_qc,
            "require_trusted_active_site": args.require_trusted_active_site,
        }
        (args.out_dir / "active_site_qc.yaml").write_text(
            yaml.safe_dump(report, sort_keys=False),
            encoding="utf-8",
        )
        raise ValueError(
            "Selected catalytic site is invalid, untrusted, or not spatially clustered. "
            "Use the active-site registry, a corrected --ser-resid, or rerun with "
            "--allow-invalid-catalytic-site only for exploratory diagnostics."
        )

    default_box_size = tuple(profile.get("default_box_size_angstrom", [25.0, 25.0, 25.0]))
    box_size = parse_box(args.box_size, default_box_size)  # type: ignore[arg-type]
    box_center = parse_box(args.box_center, ser_atom.xyz)
    exhaustiveness = args.exhaustiveness or int(profile.get("default_exhaustiveness", 8))
    num_modes = args.num_modes or int(profile.get("default_num_modes", 20))
    cpu = args.cpu or int(profile.get("default_cpu", 4))
    args.exhaustiveness = exhaustiveness
    args.num_modes = num_modes
    args.cpu = cpu

    ligand_sdf = pdbqt_dir / "ligand.sdf"
    ligand_pdbqt = pdbqt_dir / "ligand.pdbqt"
    ligand_atom_map = pdbqt_dir / "ligand_atom_map.tsv"
    ligand_atom_map_summary = pdbqt_dir / "ligand_atom_map_summary.yaml"
    receptor_prefix = pdbqt_dir / "receptor"
    receptor_pdbqt = pdbqt_dir / "receptor.pdbqt"
    box_txt = pdbqt_dir / "box.txt"
    poses_pdbqt = outputs_dir / "vina_poses.pdbqt"
    recommended_pose_pdb = outputs_dir / "recommended_pose.pdb"

    env_path = resolve_env_path(profile, root)
    tools = {
        "vina": resolve_tool(
            tool_name="vina",
            explicit=args.vina_bin,
            profile=profile,
            profile_key="vina_bin",
            env_var="AUTOPROPET_VINA_BIN",
            root=root,
        ),
        "obabel": resolve_tool(
            tool_name="obabel",
            explicit=args.obabel_bin,
            profile=profile,
            profile_key="obabel_bin",
            env_var="AUTOPROPET_OBABEL_BIN",
            root=root,
        ),
        "prepare_ligand": resolve_tool(
            tool_name="mk_prepare_ligand.py",
            explicit=args.prepare_ligand_bin,
            profile=profile,
            profile_key="prepare_ligand_bin",
            env_var="AUTOPROPET_PREPARE_LIGAND_BIN",
            root=root,
        ),
        "prepare_receptor": resolve_tool(
            tool_name="mk_prepare_receptor.py",
            explicit=args.prepare_receptor_bin,
            profile=profile,
            profile_key="prepare_receptor_bin",
            env_var="AUTOPROPET_PREPARE_RECEPTOR_BIN",
            root=root,
        ),
    }

    configure_babel_libdir_from_obabel(Path(tools["obabel"]))

    write_run_config(
        args.out_dir / "run_config.yaml",
        args=args,
        profile=profile,
        tools=tools,  # type: ignore[arg-type]
        ser_atom=ser_atom,
        box_center=box_center,
        box_size=box_size,
        active_site=active_site,
        catalytic_qc=catalytic_qc,
    )

    run_command(
        [
            str(tools["obabel"]),
            "-ipdb",
            str(ligand_copy),
            "-osdf",
            "-O",
            str(ligand_sdf),
        ],
        logs_dir / "01_obabel_ligand.log",
    )
    run_command(
        tool_command(tools["prepare_ligand"], env_path=env_path)
        + [
            "-i",
            str(ligand_sdf),
            "-o",
            str(ligand_pdbqt),
            "--charge_model",
            "gasteiger",
            "--rename_atoms",
        ],
        logs_dir / "02_prepare_ligand.log",
    )
    ligand_pdbqt_poses = parse_pdbqt_poses(ligand_pdbqt)
    ligand_mapping_summary = write_ligand_atom_map(
        charmm_atoms=parse_pdb_atoms(ligand_copy),
        pdbqt_atoms=ligand_pdbqt_poses[0].atoms,
        out_tsv=ligand_atom_map,
        summary_yaml=ligand_atom_map_summary,
    )
    run_command(
        tool_command(tools["prepare_receptor"], env_path=env_path)
        + [
            "--read_pdb",
            str(protein_copy),
            "-o",
            str(receptor_prefix),
            "-p",
            str(receptor_pdbqt),
            "-v",
            str(box_txt),
            "--box_center",
            f"{box_center[0]:.3f}",
            f"{box_center[1]:.3f}",
            f"{box_center[2]:.3f}",
            "--box_size",
            f"{box_size[0]:.3f}",
            f"{box_size[1]:.3f}",
            f"{box_size[2]:.3f}",
        ],
        logs_dir / "03_prepare_receptor.log",
    )
    run_command(
        [
            str(tools["vina"]),
            "--receptor",
            str(receptor_pdbqt),
            "--ligand",
            str(ligand_pdbqt),
            "--center_x",
            f"{box_center[0]:.3f}",
            "--center_y",
            f"{box_center[1]:.3f}",
            "--center_z",
            f"{box_center[2]:.3f}",
            "--size_x",
            f"{box_size[0]:.3f}",
            "--size_y",
            f"{box_size[1]:.3f}",
            "--size_z",
            f"{box_size[2]:.3f}",
            "--exhaustiveness",
            str(exhaustiveness),
            "--num_modes",
            str(num_modes),
            "--seed",
            str(args.seed),
            "--cpu",
            str(cpu),
            "--out",
            str(poses_pdbqt),
        ],
        logs_dir / "04_vina.log",
    )

    poses = parse_pdbqt_poses(poses_pdbqt)
    metrics = summarize_poses(
        poses,
        [atom for atom in protein_atoms if atom.is_heavy],
        ser_atom,
        args.catalytic_cutoff,
        args.max_severe_clashes,
    )
    if not metrics:
        raise RuntimeError("No pose metrics were produced")
    write_metrics(metrics, outputs_dir)

    recommended = [metric for metric in metrics if metric.recommended]
    if recommended:
        recommended_pose = next(
            pose for pose in poses if pose.index == recommended[0].pose
        )
        write_pose_as_pdb(
            recommended_pose,
            recommended_pose_pdb,
            atom_map_path=ligand_atom_map.relative_to(args.out_dir),
        )
    status = "PASS" if recommended else "WARN"
    summary = {
        "status": status,
        "out_dir": str(args.out_dir),
        "protein_id": args.protein_id,
        "ser_resid": ser_resid,
        "active_site_source": active_site.get("source") if active_site else None,
        "active_site_status": active_site.get("status") if active_site else None,
        "box_center_angstrom": list(box_center),
        "box_size_angstrom": list(box_size),
        "recommended_pose": recommended[0].pose if recommended else None,
        "recommended_pose_pdb": str(recommended_pose_pdb) if recommended else None,
        "ligand_atom_map": str(ligand_atom_map),
        "ligand_atom_map_complete": ligand_mapping_summary["mapping_complete"],
        "n_poses": len(metrics),
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
