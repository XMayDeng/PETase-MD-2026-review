#!/usr/bin/env python3
"""Validate a CHARMM complex .gro box against the chain-length-aware margin.

Reads a GROMACS .gro file, computes the bounding box of the solute (heavy
atoms only, water and ions excluded), and checks that every margin
``box_dim − bounding`` is at least ``--margin`` nm (default 1.5).

Exit codes
----------
* 0 — PASS (margin satisfied on every axis)
* 1 — FAIL (any margin below threshold, or cubic-strict check failed)
* 2 — bad input (missing file, malformed .gro, unsupported chain length)

Stdout is a single-line JSON verdict so it can be consumed by the R1
Phase B fallback retry loop::

    {"verdict": "PASS|FAIL",
     "box_dim_nm": [13.0, 13.0, 13.0],
     "bounding_per_dim_nm": [11.20, 10.80, 11.00],
     "margin_per_dim_nm":   [1.80, 2.20, 2.00],
     "margin_threshold_nm": 1.5,
     "pet_kind": "PET_L10",
     "n_solute_heavy_atoms": 223,
     "cubic_check": "PASS|N/A"}

Usage::

    python validate_complex_box.py \\
        --gro <complex>.gro \\
        --pet-kind PET_L10 \\
        [--margin 1.5] \\
        [--exclude-resnames TIP3,SOD,CLA,...] \\
        [--strict-cubic]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = REPO_ROOT / "."


def _load_box_protocol():
    """Dynamic load to avoid forcing this script to be a package import."""
    path = SIM_ROOT / "code/pipelines/box_protocol.py"
    spec = importlib.util.spec_from_file_location("box_protocol", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


BP = _load_box_protocol()


# ---------------------------------------------------------------------------
# Pure helpers (unit-test target)
# ---------------------------------------------------------------------------

def parse_gro_box(gro_text: str) -> tuple[float, float, float]:
    """Return ``(Lx, Ly, Lz)`` in nm from a GROMACS .gro box line (last line).

    Triclinic boxes have 9 values: ``v1(x) v2(y) v3(z) v1(y) v1(z) v2(x)``
    etc.  The diagonal (Lx, Ly, Lz) is the first three; for cubic and
    orthorhombic boxes the off-diagonal values are zero.  We only need
    the diagonal for the margin check.
    """
    lines = [ln for ln in gro_text.rstrip("\n").splitlines() if ln]
    if len(lines) < 3:
        raise ValueError("Truncated .gro file")
    box_line = lines[-1]
    parts = box_line.split()
    if len(parts) < 3:
        raise ValueError(f"Unparseable box line: {box_line!r}")
    return float(parts[0]), float(parts[1]), float(parts[2])


def parse_gro_atoms(gro_text: str,
                    excluded_resnames: frozenset[str] | set[str]
                    ) -> list[tuple[str, str, float, float, float]]:
    """Return heavy solute atoms as ``(resname, atom_name, x, y, z)`` tuples.

    Excludes any residue whose resname is in ``excluded_resnames``, and any
    atom whose name starts with ``H`` (hydrogens are filtered to match the
    heavy-atom envelope definition).
    """
    lines = gro_text.splitlines()
    if len(lines) < 3:
        raise ValueError("Truncated .gro file")
    try:
        n_atoms = int(lines[1].strip())
    except ValueError as exc:
        raise ValueError(f"Bad atom count line: {lines[1]!r}") from exc
    out: list[tuple[str, str, float, float, float]] = []
    for raw in lines[2:2 + n_atoms]:
        if len(raw) < 44:
            continue
        resname = raw[5:10].strip().upper()
        atom_name = raw[10:15].strip()
        if resname in excluded_resnames:
            continue
        if atom_name.upper().startswith("H"):
            continue
        try:
            x = float(raw[20:28])
            y = float(raw[28:36])
            z = float(raw[36:44])
        except ValueError:
            continue
        out.append((resname, atom_name, x, y, z))
    return out


def compute_bounding(coords: list[tuple[float, float, float]]
                     ) -> tuple[float, float, float]:
    if not coords:
        raise ValueError("Empty solute heavy-atom set; cannot compute bounding")
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def is_box_cubic(box: tuple[float, float, float], rtol: float = 1e-3) -> bool:
    """``True`` if the three box edges agree within ``rtol`` (default 0.1%)."""
    lx, ly, lz = box
    mean = (lx + ly + lz) / 3.0
    if mean <= 0:
        return False
    return all(abs(edge - mean) / mean <= rtol for edge in box)


def evaluate_margins(box: tuple[float, float, float],
                     bounding: tuple[float, float, float],
                     threshold_nm: float
                     ) -> tuple[str, list[float]]:
    margins = [b - x for b, x in zip(box, bounding)]
    verdict = "PASS" if all(m >= threshold_nm for m in margins) else "FAIL"
    return verdict, margins


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gro", required=True, type=Path,
                   help="Path to the complex .gro file (post-solvate, "
                        "post-genion).")
    p.add_argument("--pet-kind", required=True,
                   choices=sorted(BP.SUPPORTED_PET_KINDS),
                   help="PET chain length; required so the verdict can be "
                        "audited even when --margin is overridden.")
    p.add_argument("--margin", type=float, default=BP.DEFAULT_MARGIN_NM,
                   help=f"Minimum solute-to-face margin in nm "
                        f"(default {BP.DEFAULT_MARGIN_NM}).")
    p.add_argument("--exclude-resnames", type=str, default=None,
                   help="Comma-separated resname list to treat as solvent/ions "
                        "(case-insensitive). Overrides the built-in default set.")
    p.add_argument("--strict-cubic", action="store_true",
                   help="Also assert the box is cubic (all three edges equal "
                        "within 0.1 %%); fail verdict if not.")
    p.add_argument("--out", type=Path, default=None,
                   help="Optional path to also write the JSON verdict.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.gro.exists():
        print(f"validate_complex_box: missing .gro: {args.gro}", file=sys.stderr)
        return 2
    try:
        gro_text = args.gro.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        print(f"validate_complex_box: cannot decode .gro: {exc}", file=sys.stderr)
        return 2

    if args.exclude_resnames is None:
        excluded = BP.SOLVENT_AND_ION_RESNAMES
    else:
        excluded = frozenset(
            r.strip().upper() for r in args.exclude_resnames.split(",")
            if r.strip()
        )

    try:
        box = parse_gro_box(gro_text)
        atoms = parse_gro_atoms(gro_text, excluded)
    except ValueError as exc:
        print(f"validate_complex_box: {exc}", file=sys.stderr)
        return 2

    if not atoms:
        print(f"validate_complex_box: no solute heavy atoms found in {args.gro} "
              f"(every residue matched the excluded set?)", file=sys.stderr)
        return 2

    coords = [(a[2], a[3], a[4]) for a in atoms]
    bounding = compute_bounding(coords)
    verdict, margins = evaluate_margins(box, bounding, args.margin)

    cubic_check = "N/A"
    if args.strict_cubic:
        cubic_pass = is_box_cubic(box)
        cubic_check = "PASS" if cubic_pass else "FAIL"
        if not cubic_pass:
            verdict = "FAIL"

    report = {
        "verdict": verdict,
        "pet_kind": args.pet_kind,
        "box_dim_nm": [round(v, 4) for v in box],
        "bounding_per_dim_nm": [round(v, 4) for v in bounding],
        "margin_per_dim_nm": [round(v, 4) for v in margins],
        "margin_threshold_nm": args.margin,
        "expected_cubic_dim_nm": BP.box_dim_for(args.pet_kind),
        "n_solute_heavy_atoms": len(atoms),
        "cubic_check": cubic_check,
    }
    print(json.dumps(report, separators=(",", ":")))
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
