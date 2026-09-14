#!/usr/bin/env python3
"""Align CHARMM PET_L10/PET_L20 windows to a Vina PET_L4 pose.

This is the first static geometry stage after short-fragment docking. It uses
Vina coordinates only as geometry and writes transformed CHARMM PET coordinates
for downstream GROMACS building.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import yaml


MOTIF_POLICIES = [
    "full_l4_core",
    "central_2unit_core",
    "unit2_core",
    "unit3_core",
    "catalytic_ester_left",
    "catalytic_ester_right",
    "catalytic_ester_auto",
    "carbonyl_pair_only",
]


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resid: int
    x: float
    y: float
    z: float
    element: str
    record: str = "ATOM"

    @property
    def xyz(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)

    @property
    def is_heavy(self) -> bool:
        return self.element.upper() != "H" and not self.name.upper().startswith("H")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_repo_path(value: str | Path, root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def parse_pdb_atom_line(line: str) -> Atom | None:
    if not line.startswith(("ATOM", "HETATM")):
        return None
    try:
        serial = int(line[6:11])
        resid = int(line[22:26])
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None
    name = line[12:16].strip()
    element = line[76:78].strip() if len(line) >= 78 else ""
    if not element:
        element = "".join(ch for ch in name if ch.isalpha())[:1].upper()
    return Atom(
        serial=serial,
        name=name,
        resname=line[17:21].strip(),
        chain=line[21:22].strip(),
        resid=resid,
        x=x,
        y=y,
        z=z,
        element=element,
        record=line[0:6].strip() or "ATOM",
    )


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        atom = parse_pdb_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records parsed from {path}")
    return atoms


def atom_key(resid: int, name: str) -> str:
    return f"{resid}:{name}"


def parse_atom_ref(ref: str) -> tuple[int, str]:
    resid, name = ref.split(":", 1)
    return int(resid), name


def load_yaml(path: Path) -> dict:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML is not a mapping: {path}")
    return data


def load_ligand_atom_map(path: Path) -> dict[tuple[int, str], int]:
    mapping: dict[tuple[int, str], int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row.get("charmm_resid") or not row.get("charmm_atom"):
                continue
            mapping[(int(row["charmm_resid"]), row["charmm_atom"])] = int(row["map_index"])
    if not mapping:
        raise ValueError(f"No CHARMM-to-PDBQT atom identity rows parsed from {path}")
    return mapping


def atoms_by_serial(atoms: list[Atom]) -> dict[int, Atom]:
    return {atom.serial: atom for atom in atoms}


def atoms_by_ref(atoms: list[Atom]) -> dict[str, Atom]:
    return {atom_key(atom.resid, atom.name): atom for atom in atoms}


def refs_for_units(unit_map: dict, start: int, end: int) -> list[list[str]]:
    units = unit_map.get("repeat_units", [])
    by_unit = {int(item["unit"]): item for item in units}
    refs = []
    for unit in range(start, end + 1):
        if unit not in by_unit:
            raise KeyError(f"Unit {unit} not found in unit map")
        refs.append(list(by_unit[unit]["core_heavy_atoms"]))
    return refs


def unit_entries(unit_map: dict) -> dict[int, dict]:
    return {int(item["unit"]): item for item in unit_map.get("repeat_units", [])}


def names_from_refs(refs: list[str]) -> list[str]:
    return [parse_atom_ref(ref)[1] for ref in refs]


def core_atom_names(unit_map: dict, unit: int) -> list[str]:
    entries = unit_entries(unit_map)
    if unit not in entries:
        raise KeyError(f"Unit {unit} not found in unit map")
    return names_from_refs(entries[unit]["core_heavy_atoms"])


def existing_names(unit_map: dict, unit: int, names: list[str]) -> list[str]:
    available = set(core_atom_names(unit_map, unit))
    return [name for name in names if name in available]


def central_units(n_units: int, width: int = 2) -> list[int]:
    if width >= n_units:
        return list(range(1, n_units + 1))
    start = (n_units - width) // 2 + 1
    return list(range(start, start + width))


def nearest_carbonyl_to_ser(
    *,
    source_unit_map: dict,
    ligand_map: dict[tuple[int, str], int],
    vina_by_serial: dict[int, Atom],
    ser_atom: Atom | None,
) -> tuple[int, str]:
    if ser_atom is None:
        raise ValueError("--motif-policy catalytic_ester_auto requires --protein-pdb and --ser-resid")
    candidates: list[tuple[float, int, str]] = []
    for unit in range(1, int(source_unit_map["n_units"]) + 1):
        for carbonyl_name in ("C1", "C8"):
            map_index = ligand_map.get((unit, carbonyl_name))
            if map_index is None or map_index not in vina_by_serial:
                continue
            candidates.append((distance(ser_atom.xyz, vina_by_serial[map_index].xyz), unit, carbonyl_name))
    if not candidates:
        raise ValueError("No mapped PET carbonyl carbon atoms found for catalytic ester auto policy")
    _, unit, carbonyl_name = min(candidates)
    return unit, carbonyl_name


def catalytic_ester_atom_names(unit: int, carbonyl_name: str, n_units: int) -> dict[int, list[str]]:
    if carbonyl_name == "C1":
        motif = {
            unit: ["O1", "C1", "O2", "C2", "C3", "C4", "C5", "C6", "C7"],
        }
        if unit > 1:
            motif[unit - 1] = ["C9", "C10"]
        return motif
    if carbonyl_name == "C8":
        motif = {
            unit: ["C2", "C3", "C4", "C5", "C6", "C7", "C8", "O3", "O4", "C9", "C10"],
        }
        if unit < n_units:
            motif[unit + 1] = ["O1", "C1"]
        return motif
    raise ValueError(f"Unsupported catalytic carbonyl name: {carbonyl_name}")


def motif_source_atom_names(
    *,
    policy: str,
    source_unit_map: dict,
    ligand_map: dict[tuple[int, str], int],
    vina_by_serial: dict[int, Atom],
    ser_atom: Atom | None,
) -> tuple[dict[int, list[str]], dict[str, object]]:
    n_units = int(source_unit_map["n_units"])
    if policy == "full_l4_core":
        return (
            {unit: core_atom_names(source_unit_map, unit) for unit in range(1, n_units + 1)},
            {"policy": policy},
        )
    if policy == "central_2unit_core":
        units = central_units(n_units, 2)
        return (
            {unit: core_atom_names(source_unit_map, unit) for unit in units},
            {"policy": policy, "source_units": units},
        )
    if policy == "unit2_core":
        unit = min(2, n_units)
        return ({unit: core_atom_names(source_unit_map, unit)}, {"policy": policy, "source_units": [unit]})
    if policy == "unit3_core":
        unit = min(3, n_units)
        return ({unit: core_atom_names(source_unit_map, unit)}, {"policy": policy, "source_units": [unit]})
    if policy == "catalytic_ester_left":
        unit = max(1, n_units // 2)
        motif = catalytic_ester_atom_names(unit, "C8", n_units)
        return (motif, {"policy": policy, "carbonyl_unit": unit, "carbonyl_atom": "C8"})
    if policy == "catalytic_ester_right":
        unit = min(n_units, n_units // 2 + 1)
        motif = catalytic_ester_atom_names(unit, "C1", n_units)
        return (motif, {"policy": policy, "carbonyl_unit": unit, "carbonyl_atom": "C1"})
    if policy == "catalytic_ester_auto":
        unit, carbonyl_name = nearest_carbonyl_to_ser(
            source_unit_map=source_unit_map,
            ligand_map=ligand_map,
            vina_by_serial=vina_by_serial,
            ser_atom=ser_atom,
        )
        motif = catalytic_ester_atom_names(unit, carbonyl_name, n_units)
        return (motif, {"policy": policy, "carbonyl_unit": unit, "carbonyl_atom": carbonyl_name})
    if policy == "carbonyl_pair_only":
        left_unit = max(1, n_units // 2)
        right_unit = min(n_units, left_unit + 1)
        motif = {
            left_unit: ["C6", "C7", "C8", "O3", "O4"],
            right_unit: ["O1", "C1", "O2", "C2", "C3"],
        }
        return (
            motif,
            {
                "policy": policy,
                "carbonyl_pair": [f"{left_unit}:C8", f"{right_unit}:C1"],
            },
        )
    raise ValueError(f"Unsupported motif policy: {policy}")


def source_to_target_unit(source_unit: int, window: tuple[int, int], direction: str, n_source_units: int) -> int:
    target_units = list(range(window[0], window[1] + 1))
    if len(target_units) != n_source_units:
        raise ValueError(
            f"Window {window} has {len(target_units)} units, expected {n_source_units}"
        )
    if direction == "reverse":
        target_units = list(reversed(target_units))
    elif direction != "forward":
        raise ValueError(f"Unsupported direction: {direction}")
    return target_units[source_unit - 1]


def collect_alignment_points(
    *,
    vina_pose_atoms: list[Atom],
    ligand_map: dict[tuple[int, str], int],
    source_unit_map: dict,
    target_pet_atoms: list[Atom],
    target_unit_map: dict,
    window: tuple[int, int],
    direction: str,
    motif_policy: str,
    ser_atom: Atom | None = None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    vina_by_serial = atoms_by_serial(vina_pose_atoms)
    target_by_ref = atoms_by_ref(target_pet_atoms)
    n_source_units = int(source_unit_map["n_units"])
    motif_names, _ = motif_source_atom_names(
        policy=motif_policy,
        source_unit_map=source_unit_map,
        ligand_map=ligand_map,
        vina_by_serial=vina_by_serial,
        ser_atom=ser_atom,
    )
    source_coords = []
    target_coords = []
    atom_refs = []
    for source_unit, atom_names in motif_names.items():
        target_unit = source_to_target_unit(source_unit, window, direction, n_source_units)
        names = existing_names(source_unit_map, source_unit, atom_names)
        for source_atom_name in names:
            target_ref = atom_key(target_unit, source_atom_name)
            if target_ref not in target_by_ref:
                continue
            map_index = ligand_map.get((source_unit, source_atom_name))
            if map_index is None:
                continue
            vina_atom = vina_by_serial.get(map_index)
            target_atom = target_by_ref.get(target_ref)
            if vina_atom is None or target_atom is None:
                continue
            source_coords.append(vina_atom.xyz)
            target_coords.append(target_atom.xyz)
            atom_refs.append(f"{source_unit}:{source_atom_name}->{target_ref}")
    if len(source_coords) < 6:
        raise ValueError(
            f"Too few common atoms for {motif_policy} alignment: {len(source_coords)} "
            f"for window {window} {direction}"
        )
    return np.vstack(target_coords), np.vstack(source_coords), atom_refs


def kabsch_transform(source: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return R, t that map source row vectors onto target row vectors."""
    source_centroid = source.mean(axis=0)
    target_centroid = target.mean(axis=0)
    source_centered = source - source_centroid
    target_centered = target - target_centroid
    covariance = source_centered.T @ target_centered
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    translation = target_centroid - source_centroid @ rotation
    aligned = source @ rotation + translation
    rmsd = math.sqrt(float(np.mean(np.sum((aligned - target) ** 2, axis=1))))
    return rotation, translation, rmsd


def transform_atoms(atoms: list[Atom], rotation: np.ndarray, translation: np.ndarray) -> list[Atom]:
    transformed = []
    for atom in atoms:
        xyz = atom.xyz @ rotation + translation
        transformed.append(
            replace(atom, x=float(xyz[0]), y=float(xyz[1]), z=float(xyz[2]))
        )
    return transformed


def distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def find_ser_atom(protein_atoms: list[Atom], ser_resid: int, ser_atom_name: str) -> Atom | None:
    for atom in protein_atoms:
        if atom.resid == ser_resid and atom.name.upper() == ser_atom_name.upper():
            return atom
    return None


def carbonyl_c_atoms(pet_atoms: list[Atom]) -> list[Atom]:
    return [atom for atom in pet_atoms if atom.name in {"C1", "C8"} and atom.is_heavy]


def contact_counts(protein_atoms: list[Atom], pet_atoms: list[Atom]) -> dict[str, int]:
    counts = {"lt_1a": 0, "lt_2a": 0, "lt_4a": 0}
    protein_heavy = [atom for atom in protein_atoms if atom.is_heavy]
    pet_heavy = [atom for atom in pet_atoms if atom.is_heavy]
    for pet_atom in pet_heavy:
        for protein_atom in protein_heavy:
            d = distance(pet_atom.xyz, protein_atom.xyz)
            if d < 1.0:
                counts["lt_1a"] += 1
            if d < 2.0:
                counts["lt_2a"] += 1
            if d < 4.0:
                counts["lt_4a"] += 1
    return counts


def pet_resname_for_unit(unit: int, n_units: int) -> str:
    if unit == 1:
        return "PEA"
    if unit == n_units:
        return "PEB"
    return "PEM"


def pdb_line(atom: Atom, serial: int, *, resname: str | None = None, chain: str | None = None) -> str:
    name = atom.name[:4]
    output_resname = (resname or atom.resname or "UNK")[:3]
    output_chain = chain if chain is not None else (atom.chain or " ")
    element = (atom.element or name[:1] or "X")[:2].upper()
    return (
        f"{atom.record:<6}{serial:5d} {name:<4s} {output_resname:>3s} {output_chain:1s}{atom.resid:4d}    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def write_complex_pdb(
    path: Path,
    *,
    protein_atoms: list[Atom],
    pet_atoms: list[Atom],
    n_pet_units: int,
) -> None:
    lines = []
    serial = 1
    for atom in protein_atoms:
        lines.append(pdb_line(atom, serial, chain=atom.chain or "A"))
        serial += 1
    if protein_atoms:
        lines.append("TER\n")
    for atom in pet_atoms:
        lines.append(
            pdb_line(
                atom,
                serial,
                resname=pet_resname_for_unit(atom.resid, n_pet_units),
                chain="B",
            )
        )
        serial += 1
    lines.extend(["TER\n", "END\n"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="utf-8")


def candidate_name(window: tuple[int, int], direction: str) -> str:
    start, end = window
    if start == 1:
        prefix = "terminal"
    else:
        prefix = "internal_window"
    return f"{prefix}_{start}_{end}_{direction}"


def score_candidate(report: dict) -> tuple[float, float, float, float]:
    return (
        float(report["contacts"]["lt_1a"]),
        float(report["contacts"]["lt_2a"]),
        float(report["motif_rmsd_angstrom"]),
        float(report.get("ser_carbonyl_min_angstrom") or 999.0),
    )


def static_gate_passes(report: dict, *, max_motif_rmsd: float, max_severe_clashes: int) -> bool:
    return (
        float(report["motif_rmsd_angstrom"]) <= max_motif_rmsd
        and int(report["contacts"]["lt_1a"]) <= max_severe_clashes
    )


def run_alignment(args: argparse.Namespace) -> dict:
    root = repo_root()
    out_dir = resolve_repo_path(args.out_dir, root)
    if out_dir.exists():
        if not args.force:
            raise FileExistsError(f"Output directory already exists: {out_dir}")
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    vina_pose_pdb = resolve_repo_path(args.vina_pose_pdb, root)
    ligand_atom_map = resolve_repo_path(args.ligand_atom_map, root)
    source_unit_map_path = resolve_repo_path(args.source_unit_map, root)
    target_pet_pdb = resolve_repo_path(args.target_pet_pdb, root)
    target_unit_map_path = resolve_repo_path(args.target_unit_map, root)
    protein_pdb = resolve_repo_path(args.protein_pdb, root) if args.protein_pdb else None

    vina_pose_atoms = parse_pdb_atoms(vina_pose_pdb)
    target_pet_atoms = parse_pdb_atoms(target_pet_pdb)
    protein_atoms = parse_pdb_atoms(protein_pdb) if protein_pdb else []
    ligand_map = load_ligand_atom_map(ligand_atom_map)
    source_unit_map = load_yaml(source_unit_map_path)
    target_unit_map = load_yaml(target_unit_map_path)
    windows = args.window or target_unit_map["candidate_windows"]["limited_first_pass"]
    windows = [tuple(int(value) for value in item) for item in windows]

    ser_atom = (
        find_ser_atom(protein_atoms, args.ser_resid, args.ser_atom_name)
        if protein_atoms and args.ser_resid
        else None
    )
    vina_by_serial = atoms_by_serial(vina_pose_atoms)
    _, motif_metadata = motif_source_atom_names(
        policy=args.motif_policy,
        source_unit_map=source_unit_map,
        ligand_map=ligand_map,
        vina_by_serial=vina_by_serial,
        ser_atom=ser_atom,
    )
    reports = []
    for window in windows:
        for direction in ("forward", "reverse"):
            target_points, source_points, atom_refs = collect_alignment_points(
                vina_pose_atoms=vina_pose_atoms,
                ligand_map=ligand_map,
                source_unit_map=source_unit_map,
                target_pet_atoms=target_pet_atoms,
                target_unit_map=target_unit_map,
                window=window,
                direction=direction,
                motif_policy=args.motif_policy,
                ser_atom=ser_atom,
            )
            rotation, translation, rmsd = kabsch_transform(target_points, source_points)
            transformed_pet = transform_atoms(target_pet_atoms, rotation, translation)
            contacts = contact_counts(protein_atoms, transformed_pet) if protein_atoms else {
                "lt_1a": 0,
                "lt_2a": 0,
                "lt_4a": 0,
            }
            ser_carbonyl = None
            if ser_atom is not None:
                carbonyl_atoms = carbonyl_c_atoms(transformed_pet)
                if carbonyl_atoms:
                    ser_carbonyl = min(distance(ser_atom.xyz, atom.xyz) for atom in carbonyl_atoms)
            name = candidate_name(window, direction)
            candidate_dir = out_dir / "candidates" / name
            complex_pdb = candidate_dir / "protein_pet_complex.pdb"
            report_path = candidate_dir / "alignment_report.yaml"
            write_complex_pdb(
                complex_pdb,
                protein_atoms=protein_atoms,
                pet_atoms=transformed_pet,
                n_pet_units=int(target_unit_map["n_units"]),
            )
            report = {
                "candidate": name,
                "window": list(window),
                "direction": direction,
                "motif_policy": args.motif_policy,
                "motif_metadata": motif_metadata,
                "motif_rmsd_angstrom": round(rmsd, 4),
                "n_alignment_atoms": len(atom_refs),
                "contacts": contacts,
                "ser_resid": args.ser_resid,
                "ser_atom": args.ser_atom_name,
                "ser_carbonyl_min_angstrom": round(ser_carbonyl, 4) if ser_carbonyl is not None else None,
                "complex_pdb": str(complex_pdb.relative_to(out_dir)),
                "alignment_atom_refs": atom_refs,
            }
            report["static_gate_pass"] = static_gate_passes(
                report,
                max_motif_rmsd=args.max_motif_rmsd,
                max_severe_clashes=args.max_severe_clashes,
            )
            report_path.write_text(yaml.safe_dump(report, sort_keys=False), encoding="utf-8")
            reports.append(report)

    selected = min(reports, key=score_candidate)
    selected_passes = static_gate_passes(
        selected,
        max_motif_rmsd=args.max_motif_rmsd,
        max_severe_clashes=args.max_severe_clashes,
    )
    selected_dir = out_dir / "selected"
    selected_dir.mkdir(parents=True, exist_ok=True)
    selected_src = out_dir / selected["complex_pdb"]
    shutil.copy2(selected_src, selected_dir / "protein_pet_complex.pdb")
    (selected_dir / "alignment_report.yaml").write_text(
        yaml.safe_dump(selected, sort_keys=False),
        encoding="utf-8",
    )
    summary = {
        "status": "PASS" if selected_passes else "WARN",
        "selected_candidate": selected["candidate"],
        "selected_window": selected["window"],
        "selected_direction": selected["direction"],
        "selected_static_gate_pass": selected_passes,
        "n_candidates": len(reports),
        "motif_policy": args.motif_policy,
        "motif_metadata": motif_metadata,
        "ranking_policy": "lt_1a, lt_2a, motif_rmsd, ser_carbonyl_min",
        "static_gate": {
            "max_motif_rmsd_angstrom": args.max_motif_rmsd,
            "max_severe_clashes_lt_1a": args.max_severe_clashes,
        },
        "inputs": {
            "vina_pose_pdb": str(vina_pose_pdb),
            "ligand_atom_map": str(ligand_atom_map),
            "source_unit_map": str(source_unit_map_path),
            "target_pet_pdb": str(target_pet_pdb),
            "target_unit_map": str(target_unit_map_path),
            "protein_pdb": str(protein_pdb) if protein_pdb else None,
        },
        "candidates": [
            {
                key: report[key]
                for key in (
                    "candidate",
                    "window",
                    "direction",
                    "motif_policy",
                    "motif_rmsd_angstrom",
                    "n_alignment_atoms",
                    "contacts",
                    "ser_carbonyl_min_angstrom",
                    "static_gate_pass",
                    "complex_pdb",
                )
            }
            for report in sorted(reports, key=score_candidate)
        ],
    }
    (out_dir / "summary.yaml").write_text(yaml.safe_dump(summary, sort_keys=False), encoding="utf-8")
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def parse_window(value: str) -> list[int]:
    if "-" not in value:
        raise argparse.ArgumentTypeError("Window must be formatted like START-END")
    start, end = value.split("-", 1)
    return [int(start), int(end)]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Align CHARMM PET windows to a Vina PET_L4 recommended pose."
    )
    parser.add_argument("--vina-pose-pdb", required=True, type=Path)
    parser.add_argument("--ligand-atom-map", required=True, type=Path)
    parser.add_argument("--source-unit-map", required=True, type=Path)
    parser.add_argument("--target-pet-pdb", required=True, type=Path)
    parser.add_argument("--target-unit-map", required=True, type=Path)
    parser.add_argument("--protein-pdb", type=Path)
    parser.add_argument("--ser-resid", type=int)
    parser.add_argument("--ser-atom-name", default="OG")
    parser.add_argument(
        "--motif-policy",
        choices=MOTIF_POLICIES,
        default="catalytic_ester_auto",
        help=(
            "Subset of PET_L4 atoms used to place the long-chain CHARMM PET. "
            "Default catalytic_ester_auto focuses on the ester motif closest "
            "to the active-site serine; full_l4_core preserves the original "
            "full-fragment benchmark behavior."
        ),
    )
    parser.add_argument("--max-motif-rmsd", type=float, default=3.0)
    parser.add_argument("--max-severe-clashes", type=int, default=0)
    parser.add_argument("--window", action="append", type=parse_window)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = run_alignment(args)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
