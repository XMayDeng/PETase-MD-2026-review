#!/usr/bin/env python3
"""Build PET repeat-unit atom maps for Vina-to-CHARMM alignment.

The map is intentionally force-field-oriented: it records CHARMM atom identity
from the frozen PET asset so Vina pose coordinates can later be mapped back to
CHARMM PET_L10/PET_L20 topology without trusting PDBQT atom names alone.
"""
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import yaml


SIM_ROOT_REL = "."
DEFAULT_RAW_PDB_REL = "raw/charmm-gui-001/gromacs/step3_input.pdb"

CORE_HEAVY_ORDER = [
    "O1",
    "C1",
    "O2",
    "C2",
    "C3",
    "C4",
    "C5",
    "C6",
    "C7",
    "C8",
    "O3",
    "O4",
    "C9",
    "C10",
]
AROMATIC_RING_ATOMS = ["C2", "C3", "C4", "C5", "C6", "C7"]
ESTER_CARBONYL_ATOMS = ["C1", "O2", "C8", "O4"]
GLYCOL_LINKER_ATOMS = ["O1", "O3", "C9", "C10"]
TERMINAL_CAP_ATOMS = {"H11", "OZ", "HZ1"}


@dataclass(frozen=True)
class PdbAtom:
    serial: int
    name: str
    resname: str
    chain: str
    resid: int
    x: float
    y: float
    z: float
    element: str

    @property
    def is_hydrogen(self) -> bool:
        return self.element.upper() == "H" or self.name.upper().startswith("H")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_repo_path(value: str | Path, root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def parse_pdb_atom_line(line: str) -> PdbAtom | None:
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
        element = re.sub(r"[^A-Za-z]", "", name)[:1].upper()
    return PdbAtom(
        serial=serial,
        name=name,
        resname=line[17:21].strip(),
        chain=line[21:22].strip(),
        resid=resid,
        x=x,
        y=y,
        z=z,
        element=element,
    )


def parse_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        atom = parse_pdb_atom_line(line)
        if atom is not None:
            atoms.append(atom)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records parsed from {path}")
    return atoms


def pet_kind_unit_count(pet_kind: str) -> int:
    match = re.fullmatch(r"PET_?L?(\d+)", pet_kind.upper())
    if not match:
        raise ValueError(f"Cannot infer PET unit count from kind: {pet_kind}")
    return int(match.group(1))


def atom_ref(atom: PdbAtom) -> str:
    return f"{atom.resid}:{atom.name}"


def atoms_by_residue(atoms: list[PdbAtom]) -> dict[int, list[PdbAtom]]:
    grouped: dict[int, list[PdbAtom]] = {}
    for atom in atoms:
        grouped.setdefault(atom.resid, []).append(atom)
    return dict(sorted(grouped.items()))


def pick_names(unit_atoms: list[PdbAtom], ordered_names: list[str]) -> list[str]:
    available = {atom.name for atom in unit_atoms}
    return [name for name in ordered_names if name in available]


def build_default_windows(n_units: int, window_size: int = 4) -> dict[str, list[list[int]]]:
    if n_units < window_size:
        return {"limited_first_pass": [], "all_windows": []}
    all_windows = [[start, start + window_size - 1] for start in range(1, n_units - window_size + 2)]
    limited: list[list[int]] = []
    if n_units == window_size:
        limited = [[1, n_units]]
    else:
        starts = [1, max(1, (n_units - window_size) // 2 + 1), n_units - window_size + 1]
        for start in starts:
            item = [start, start + window_size - 1]
            if item not in limited:
                limited.append(item)
    return {"limited_first_pass": limited, "all_windows": all_windows}


def build_unit_map(
    atoms: list[PdbAtom],
    *,
    pet_kind: str,
    source_pdb: Path | None = None,
    asset_root: Path | None = None,
) -> dict:
    expected_units = pet_kind_unit_count(pet_kind)
    grouped = atoms_by_residue(atoms)
    residues = list(grouped)
    if len(residues) != expected_units:
        raise ValueError(
            f"{pet_kind} expected {expected_units} PETL residues, found {len(residues)}"
        )

    units = []
    missing_by_unit: dict[int, list[str]] = {}
    for idx, resid in enumerate(residues, start=1):
        unit_atoms = grouped[resid]
        by_name = {atom.name: atom for atom in unit_atoms}
        core_names = pick_names(unit_atoms, CORE_HEAVY_ORDER)
        missing = [name for name in CORE_HEAVY_ORDER if name not in by_name]
        if missing:
            missing_by_unit[idx] = missing
        terminal_cap_names = [
            atom.name
            for atom in unit_atoms
            if atom.name in TERMINAL_CAP_ATOMS
        ]
        core_refs = [
            atom_ref(by_name[name])
            for name in core_names
            if name in by_name and not by_name[name].is_hydrogen and name not in TERMINAL_CAP_ATOMS
        ]
        terminal_position = "internal"
        if idx == 1:
            terminal_position = "first"
        elif idx == expected_units:
            terminal_position = "last"
        units.append(
            {
                "unit": idx,
                "resid": resid,
                "resname": unit_atoms[0].resname,
                "terminal_position": terminal_position,
                "all_atom_refs": [atom_ref(atom) for atom in unit_atoms],
                "core_heavy_atoms": core_refs,
                "aromatic_ring_atoms": [
                    atom_ref(by_name[name])
                    for name in AROMATIC_RING_ATOMS
                    if name in by_name
                ],
                "ester_carbonyl_atoms": [
                    atom_ref(by_name[name])
                    for name in ESTER_CARBONYL_ATOMS
                    if name in by_name
                ],
                "glycol_linker_atoms": [
                    atom_ref(by_name[name])
                    for name in GLYCOL_LINKER_ATOMS
                    if name in by_name
                ],
                "terminal_cap_atoms": [
                    atom_ref(by_name[name])
                    for name in terminal_cap_names
                    if name in by_name
                ],
            }
        )

    if missing_by_unit:
        details = "; ".join(
            f"unit {unit}: {','.join(names)}" for unit, names in missing_by_unit.items()
        )
        raise ValueError(f"Missing required PET core heavy atom names: {details}")

    data = {
        "schema_version": 1,
        "pet_kind": pet_kind,
        "n_units": expected_units,
        "source_pdb": str(source_pdb) if source_pdb is not None else None,
        "asset_root": str(asset_root) if asset_root is not None else None,
        "residue_name": atoms[0].resname,
        "atom_name_policy": {
            "identity": "resid:atom_name",
            "alignment_core": "common repeat-unit core heavy atoms",
            "exclude": ["terminal_cap_atoms", "hydrogens"],
        },
        "unit_atom_name_sets": {
            "core_heavy_order": CORE_HEAVY_ORDER,
            "aromatic_ring_atoms": AROMATIC_RING_ATOMS,
            "ester_carbonyl_atoms": ESTER_CARBONYL_ATOMS,
            "glycol_linker_atoms": GLYCOL_LINKER_ATOMS,
            "terminal_cap_atom_names": sorted(TERMINAL_CAP_ATOMS),
        },
        "candidate_windows": build_default_windows(expected_units),
        "repeat_units": units,
    }
    return data


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build PET repeat-unit maps for geometry-guided alignment."
    )
    parser.add_argument("--pet-kind", required=True, help="PET kind, e.g. PET_L4.")
    parser.add_argument(
        "--asset-root",
        type=Path,
        help="PET asset root, e.g. Simulation/inputs/pets/PET_L4.",
    )
    parser.add_argument("--source-pdb", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = repo_root()
    asset_root = (
        resolve_repo_path(args.asset_root, root)
        if args.asset_root
        else root / SIM_ROOT_REL / "inputs" / "pets" / args.pet_kind
    )
    source_pdb = (
        resolve_repo_path(args.source_pdb, root)
        if args.source_pdb
        else asset_root / DEFAULT_RAW_PDB_REL
    )
    out_path = (
        resolve_repo_path(args.out, root)
        if args.out
        else asset_root / "alignment" / f"{args.pet_kind.lower()}_unit_map.yaml"
    )
    if out_path.exists() and not args.force:
        raise FileExistsError(f"Output already exists: {out_path}")

    atoms = parse_pdb_atoms(source_pdb)
    try:
        source_rel = source_pdb.relative_to(asset_root)
    except ValueError:
        source_rel = source_pdb
    data = build_unit_map(
        atoms,
        pet_kind=args.pet_kind,
        source_pdb=source_rel,
        asset_root=asset_root.relative_to(root) if asset_root.is_relative_to(root) else asset_root,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    print(out_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
