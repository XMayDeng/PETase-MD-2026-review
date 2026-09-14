#!/usr/bin/env python3
"""
Sanity checks for system_build geometry/topology before long MD runs.

Default target:
  - code/system_build/conf_pet.gro
  - code/system_build/P1.itp

What this checks:
  1) Box size and atoms outside [0, box)
  2) Protein-ligand minimum distance under minimum-image convention
  3) Protein distance to box faces (wrapped coordinates)
  4) Ligand size vs box size (end-to-end and axis span)
  5) Ligand bond lengths from .itp against current .gro coordinates

This script is heuristic. It is designed to catch obvious setup risks early.
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path


WATER_RESNAMES = {"SOL", "WAT", "HOH"}
ION_RESNAMES = {"NA", "CL", "K", "CA", "MG", "ZN", "SOD", "CLA", "POT", "LIT", "CES", "RUB"}


@dataclass
class Atom:
    resnr: int
    resname: str
    atomname: str
    atomnr: int
    x: float
    y: float
    z: float


@dataclass
class CheckResult:
    level: str
    name: str
    detail: str


def read_gro(path: Path) -> tuple[list[Atom], tuple[float, float, float]]:
    if not path.exists():
        raise SystemExit(f"GRO file not found: {path}")

    atoms: list[Atom] = []
    with path.open("r", encoding="utf-8") as f:
        _ = f.readline()  # title
        natoms_line = f.readline().strip()
        try:
            natoms = int(natoms_line)
        except ValueError as exc:
            raise SystemExit(f"Invalid GRO atom count in {path}: {natoms_line!r}") from exc

        for _ in range(natoms):
            line = f.readline()
            if len(line) < 44:
                raise SystemExit(f"Invalid GRO atom line in {path}: {line!r}")
            atoms.append(
                Atom(
                    resnr=int(line[0:5]),
                    resname=line[5:10].strip(),
                    atomname=line[10:15].strip(),
                    atomnr=int(line[15:20]),
                    x=float(line[20:28]),
                    y=float(line[28:36]),
                    z=float(line[36:44]),
                )
            )

        box_line = f.readline().split()
        if len(box_line) < 3:
            raise SystemExit(f"Invalid GRO box line in {path}: {' '.join(box_line)!r}")
        box = (float(box_line[0]), float(box_line[1]), float(box_line[2]))

    if any(v <= 0.0 for v in box):
        raise SystemExit(f"Invalid box size in {path}: {box}")
    return atoms, box


def parse_itp_bonds(path: Path) -> tuple[int, list[tuple[int, int]]]:
    if not path.exists():
        raise SystemExit(f"ITP file not found: {path}")

    in_atoms = False
    in_bonds = False
    natoms = 0
    bonds: list[tuple[int, int]] = []

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            lower = line.lower()
            in_atoms = lower.startswith("[ atoms ]")
            in_bonds = lower.startswith("[ bonds ]")
            continue

        if in_atoms:
            parts = line.split()
            if parts and parts[0].isdigit():
                natoms = max(natoms, int(parts[0]))
            continue

        if in_bonds:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                # convert 1-based to 0-based
                bonds.append((int(parts[0]) - 1, int(parts[1]) - 1))

    if natoms <= 0:
        raise SystemExit(f"No [ atoms ] parsed from {path}")
    if not bonds:
        raise SystemExit(f"No [ bonds ] parsed from {path}")
    return natoms, bonds


def wrap(v: float, box: float) -> float:
    out = v % box
    if out < 0:
        out += box
    return out


def dist3(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def dist3_mic(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    box: tuple[float, float, float],
) -> float:
    lx, ly, lz = box
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    dx -= lx * round(dx / lx)
    dy -= ly * round(dy / ly)
    dz -= lz * round(dz / lz)
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def min_face_distance(coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> float:
    if not coords:
        return float("nan")
    lx, ly, lz = box
    best = float("inf")
    for x, y, z in coords:
        wx = wrap(x, lx)
        wy = wrap(y, ly)
        wz = wrap(z, lz)
        best = min(best, wx, lx - wx, wy, ly - wy, wz, lz - wz)
    return best


def count_outside(coords: list[tuple[float, float, float]], box: tuple[float, float, float]) -> int:
    lx, ly, lz = box
    n = 0
    for x, y, z in coords:
        if x < 0.0 or x >= lx or y < 0.0 or y >= ly or z < 0.0 or z >= lz:
            n += 1
    return n


def min_group_distance_mic(
    a_coords: list[tuple[float, float, float]],
    b_coords: list[tuple[float, float, float]],
    box: tuple[float, float, float],
) -> float:
    best = float("inf")
    for a in a_coords:
        for b in b_coords:
            d = dist3_mic(a, b, box)
            if d < best:
                best = d
    return best


def max_pair_distance(coords: list[tuple[float, float, float]]) -> float:
    best = 0.0
    n = len(coords)
    for i in range(n):
        for j in range(i + 1, n):
            d = dist3(coords[i], coords[j])
            if d > best:
                best = d
    return best


def axis_span(coords: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    return max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)


def parse_args() -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]  # generalGromacs_PL

    p = argparse.ArgumentParser()
    p.add_argument(
        "--gro",
        type=Path,
        default=root / "code/system_build/conf_pet.gro",
        help="Input GRO to validate",
    )
    p.add_argument(
        "--itp",
        type=Path,
        default=root / "code/system_build/P1.itp",
        help="Ligand ITP used for bond-length sanity checks",
    )
    p.add_argument("--ligand-resname", default="PETL", help="Ligand residue name in GRO")
    p.add_argument(
        "--protein-face-min",
        type=float,
        default=0.6,
        help="Warn if wrapped protein min distance to box faces is below this (nm)",
    )
    p.add_argument(
        "--pl-clash-min",
        type=float,
        default=0.12,
        help="Fail if protein-ligand minimum distance is below this (nm)",
    )
    p.add_argument(
        "--pl-contact-max",
        type=float,
        default=0.8,
        help="Warn if protein-ligand minimum distance is above this (nm)",
    )
    p.add_argument(
        "--bond-max",
        type=float,
        default=0.25,
        help="Fail if any ligand bond is longer than this (nm)",
    )
    p.add_argument(
        "--strict",
        action="store_true",
        help="Return exit code 1 when WARN/FAIL exists (default: only FAIL returns 1)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    atoms, box = read_gro(args.gro)

    lig_name = args.ligand_resname.upper()
    ligand_coords: list[tuple[float, float, float]] = []
    protein_coords: list[tuple[float, float, float]] = []
    solvent_coords: list[tuple[float, float, float]] = []

    for a in atoms:
        res = a.resname.upper()
        xyz = (a.x, a.y, a.z)
        if res == lig_name:
            ligand_coords.append(xyz)
        elif res in WATER_RESNAMES or res in ION_RESNAMES:
            solvent_coords.append(xyz)
        else:
            protein_coords.append(xyz)

    results: list[CheckResult] = []

    def add(level: str, name: str, detail: str) -> None:
        results.append(CheckResult(level=level, name=name, detail=detail))

    lx, ly, lz = box
    box_min = min(box)
    add("INFO", "Input", f"gro={args.gro}")
    add("INFO", "Box", f"Lx={lx:.3f} Ly={ly:.3f} Lz={lz:.3f} nm")
    add(
        "INFO",
        "Groups",
        f"protein_atoms={len(protein_coords)} ligand_atoms={len(ligand_coords)} solvent_or_ion_atoms={len(solvent_coords)}",
    )

    if not protein_coords:
        add("FAIL", "Protein group", "No protein-like atoms detected in GRO")
    if not ligand_coords:
        add("FAIL", "Ligand group", f"No atoms found with ligand resname={lig_name}")

    prot_outside = count_outside(protein_coords, box) if protein_coords else 0
    lig_outside = count_outside(ligand_coords, box) if ligand_coords else 0
    if prot_outside > 0:
        add("WARN", "Protein coordinates", f"{prot_outside} protein atoms are outside [0, box)")
    else:
        add("OK", "Protein coordinates", "All protein atoms are inside [0, box)")
    if lig_outside > 0:
        add("WARN", "Ligand coordinates", f"{lig_outside} ligand atoms are outside [0, box) (crossing PBC in coordinate file)")
    else:
        add("OK", "Ligand coordinates", "All ligand atoms are inside [0, box)")

    if protein_coords:
        prot_face = min_face_distance(protein_coords, box)
        if prot_face < args.protein_face_min:
            add(
                "WARN",
                "Protein-box clearance",
                f"wrapped minimum protein-to-face distance={prot_face:.3f} nm < {args.protein_face_min:.3f} nm",
            )
        else:
            add(
                "OK",
                "Protein-box clearance",
                f"wrapped minimum protein-to-face distance={prot_face:.3f} nm",
            )

    if protein_coords and ligand_coords:
        pl_min = min_group_distance_mic(protein_coords, ligand_coords, box)
        if pl_min < args.pl_clash_min:
            add(
                "FAIL",
                "Protein-ligand minimum distance",
                f"{pl_min:.3f} nm < clash threshold {args.pl_clash_min:.3f} nm",
            )
        elif pl_min > args.pl_contact_max:
            add(
                "WARN",
                "Protein-ligand minimum distance",
                f"{pl_min:.3f} nm > contact threshold {args.pl_contact_max:.3f} nm (initially weak/no contact)",
            )
        else:
            add("OK", "Protein-ligand minimum distance", f"{pl_min:.3f} nm")

    if ligand_coords:
        lig_e2e = dist3(ligand_coords[0], ligand_coords[-1])
        sx, sy, sz = axis_span(ligand_coords)
        span_max = max(sx, sy, sz)
        add("INFO", "Ligand geometry", f"end_to_end={lig_e2e:.3f} nm axis_span=({sx:.3f}, {sy:.3f}, {sz:.3f}) nm")

        if lig_e2e > box_min:
            add(
                "WARN",
                "Ligand length vs box",
                f"ligand end-to-end {lig_e2e:.3f} nm > min box length {box_min:.3f} nm",
            )
        else:
            add("OK", "Ligand length vs box", f"ligand end-to-end {lig_e2e:.3f} nm <= min box length {box_min:.3f} nm")

        if span_max > box_min:
            add(
                "WARN",
                "Ligand span vs box",
                f"ligand max coordinate span {span_max:.3f} nm > min box length {box_min:.3f} nm",
            )
        else:
            add("OK", "Ligand span vs box", f"ligand max coordinate span {span_max:.3f} nm <= min box length {box_min:.3f} nm")

    # Bond sanity from ITP (if available)
    if args.itp.exists() and ligand_coords:
        natoms_itp, bonds = parse_itp_bonds(args.itp)
        if natoms_itp != len(ligand_coords):
            add(
                "WARN",
                "Ligand atom count",
                f"ITP atoms={natoms_itp}, ligand atoms in GRO={len(ligand_coords)} (bond check skipped)",
            )
        else:
            bond_lengths = [dist3(ligand_coords[i], ligand_coords[j]) for i, j in bonds]
            max_bond = max(bond_lengths)
            mean_bond = sum(bond_lengths) / len(bond_lengths)
            if max_bond > args.bond_max:
                add(
                    "FAIL",
                    "Ligand bond length",
                    f"max={max_bond:.3f} nm > threshold {args.bond_max:.3f} nm (mean={mean_bond:.3f} nm)",
                )
            else:
                add(
                    "OK",
                    "Ligand bond length",
                    f"max={max_bond:.3f} nm mean={mean_bond:.3f} nm",
                )
    elif not args.itp.exists():
        add("WARN", "Ligand bond length", f"ITP not found, skip bond check: {args.itp}")

    fail_count = sum(1 for r in results if r.level == "FAIL")
    warn_count = sum(1 for r in results if r.level == "WARN")

    print("=== system_build sanity report ===")
    for r in results:
        print(f"[{r.level:<4}] {r.name}: {r.detail}")
    print("----------------------------------")
    print(f"Summary: FAIL={fail_count} WARN={warn_count}")

    if fail_count > 0:
        return 1
    if args.strict and warn_count > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
