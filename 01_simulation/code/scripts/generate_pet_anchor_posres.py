#!/usr/bin/env python3
"""Generate POSRES include for the L4 anchor of a long-chain PET case.

Spec: ``docs/plans/2026-05-03_md_workflow_v1_binding_attribution_spec.md``
§4.1 Phase C + PA-B (2026-05-13).

In R1 Phase C (``04_chain_relax`` stage, L10/L20 only) the protein backbone
and the 4 PETL residues forming the docked L4 anchor are position-
restrained while the chain-extension residues relax freely. The anchor
atom set is direction-dependent (``forward`` / ``reverse``) because the
Kabsch alignment can flip the residue order along the docked motif.

This script reads:

  * the long-chain target ``unit_map.yaml`` for the per-residue atom
    name list (``repeat_units[].all_atom_refs``), and
  * the CHARMM-GUI ``P1.itp`` topology for the ``[atoms]`` section that
    maps ``(resid, atom_name)`` to a 1-indexed atom serial within the
    PET molecule.

It emits ``posre_pet_anchor.itp`` with one position-restraint line per
selected atom at the requested force constant (default K=500
kJ·mol⁻¹·nm⁻², spec §4.1).

Usage::

    python generate_pet_anchor_posres.py \\
        --unit-map     inputs/pets/PET_L10/alignment/pet_l10_unit_map.yaml \\
        --p1-itp       inputs/pets/PET_L10/raw/charmm-gui-001/gromacs/toppar/P1.itp \\
        --window       1 4 \\
        --direction    forward \\
        --output       <build-dir>/dir_tail-side/outputs/posre_pet_anchor.itp \\
        [--posres-k 500]
"""
from __future__ import annotations

import argparse
import datetime
import re
import sys
from pathlib import Path

import yaml

DEFAULT_POSRES_K = 500
SUPPORTED_DIRECTIONS = ("forward", "reverse")


# ---------------------------------------------------------------------------
# Pure helpers (unit-test target)
# ---------------------------------------------------------------------------

def parse_p1_itp_atoms(itp_text: str) -> dict[tuple[int, str], int]:
    """Parse the ``[ atoms ]`` section of a CHARMM-GUI P1.itp.

    Returns ``{(resid, atom_name): nr}``. ``nr`` is the 1-indexed atom
    serial GROMACS uses inside POSRES ``[ position_restraints ]``.
    """
    out: dict[tuple[int, str], int] = {}
    in_atoms = False
    for raw in itp_text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            in_atoms = (line.strip("[] ").strip() == "atoms")
            continue
        if not in_atoms:
            continue
        # Layout: nr  type  resnr  residu  atom  cgnr  charge  mass [; ...]
        body = line.split(";", 1)[0]
        parts = body.split()
        if len(parts) < 5:
            continue
        try:
            nr = int(parts[0])
            resnr = int(parts[2])
        except ValueError:
            continue
        atom_name = parts[4]
        out[(resnr, atom_name)] = nr
    return out


_ATOM_REF_RE = re.compile(r"^\s*(\d+)\s*:\s*(\S+)\s*$")


def parse_atom_ref(ref: str) -> tuple[int, str]:
    """Parse a ``"<resid>:<atom_name>"`` ref from ``unit_map.yaml``."""
    m = _ATOM_REF_RE.match(ref)
    if not m:
        raise ValueError(f"Bad atom ref {ref!r}; expected '<resid>:<atom_name>'")
    return int(m.group(1)), m.group(2)


def unit_atom_refs_for_residues(unit_map: dict, residues: list[int]) -> list[str]:
    """Return concatenated ``all_atom_refs`` from ``unit_map`` repeat_units
    for the residue ids in ``residues``.

    Raises ``ValueError`` if any requested residue id is missing.
    """
    repeat_units = unit_map.get("repeat_units")
    if not isinstance(repeat_units, list):
        raise ValueError("unit_map.repeat_units missing or malformed")
    by_resid: dict[int, dict] = {}
    for entry in repeat_units:
        rid = entry.get("resid")
        if rid is not None:
            by_resid[int(rid)] = entry
    refs: list[str] = []
    for rid in residues:
        if rid not in by_resid:
            raise ValueError(
                f"unit_map has no repeat_unit for resid={rid}; "
                f"available = {sorted(by_resid)}"
            )
        atom_refs = by_resid[rid].get("all_atom_refs")
        if not isinstance(atom_refs, list):
            raise ValueError(f"resid={rid} all_atom_refs missing")
        refs.extend(str(r) for r in atom_refs)
    return refs


def resolve_anchor_atom_indices(
    *,
    unit_map: dict,
    p1_atoms_index: dict[tuple[int, str], int],
    window: tuple[int, int],
    direction: str,
) -> list[int]:
    """Return the sorted 1-indexed POSRES atom list for one (window, direction).

    Implements spec §4.1 PA-B pseudo-code: build the residue list, reverse
    when ``direction == "reverse"``, gather every atom ref (heavy + H), look
    up each in the P1.itp [atoms] map, sort ascending.
    """
    if direction not in SUPPORTED_DIRECTIONS:
        raise ValueError(
            f"direction={direction!r}; expected one of {SUPPORTED_DIRECTIONS}"
        )
    start, end = int(window[0]), int(window[1])
    if start > end:
        raise ValueError(f"Bad window {window!r}: start > end")
    residues = list(range(start, end + 1))
    if direction == "reverse":
        residues = list(reversed(residues))

    refs = unit_atom_refs_for_residues(unit_map, residues)
    serials: list[int] = []
    missing: list[str] = []
    for ref in refs:
        resid, atom_name = parse_atom_ref(ref)
        nr = p1_atoms_index.get((resid, atom_name))
        if nr is None:
            missing.append(ref)
            continue
        serials.append(nr)
    if missing:
        raise ValueError(
            f"P1.itp [atoms] is missing {len(missing)} entries for window "
            f"{window} direction={direction}: e.g. {missing[:5]}"
        )
    return sorted(set(serials))


def render_posre_itp(
    *,
    atom_indices: list[int],
    posres_k: int,
    case_id: str | None,
    direction: str,
    window: tuple[int, int],
) -> str:
    """Render a GROMACS ``[ position_restraints ]`` itp body.

    Format-stable; reproducible across runs (deterministic from inputs)."""
    lines: list[str] = []
    lines.append("; posre_pet_anchor.itp - L4 anchor POSRES")
    lines.append("; generated by generate_pet_anchor_posres.py "
                 "(spec section 4.1 / PA-B)")
    if case_id:
        lines.append(f"; case_id = {case_id}")
    lines.append(f"; direction = {direction}")
    lines.append(f"; window = [{window[0]}, {window[1]}]")
    lines.append(f"; n_atoms = {len(atom_indices)}")
    lines.append(f"; force_constant_kJ_per_mol_per_nm2 = {posres_k}")
    lines.append(f"; generated_at_utc = "
                 f"{datetime.datetime.now(datetime.timezone.utc).isoformat()}")
    lines.append("")
    lines.append("[ position_restraints ]")
    lines.append(";  atom_id    type        fx        fy        fz")
    for nr in atom_indices:
        lines.append(f"{nr:>10d} {1:>10d} {posres_k:>10d} {posres_k:>10d} {posres_k:>10d}")
    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--unit-map", required=True, type=Path,
                   help="Long-chain target unit_map.yaml")
    p.add_argument("--p1-itp", required=True, type=Path,
                   help="CHARMM-GUI P1.itp (toppar/P1.itp) for the long-chain target")
    p.add_argument("--window", nargs=2, required=True, type=int,
                   metavar=("START", "END"),
                   help="L4 anchor window (1-indexed inclusive residue range)")
    p.add_argument("--direction", required=True, choices=SUPPORTED_DIRECTIONS,
                   help="Kabsch atom-order direction; reverses residue order")
    p.add_argument("--output", required=True, type=Path,
                   help="Destination .itp path")
    p.add_argument("--posres-k", type=int, default=DEFAULT_POSRES_K,
                   help=f"POSRES force constant on x/y/z (default {DEFAULT_POSRES_K} "
                        f"kJ/mol/nm^2 per spec section 4.1)")
    p.add_argument("--case-id", default=None,
                   help="Optional case identifier embedded as a comment for audit")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.unit_map.exists():
        print(f"generate_pet_anchor_posres: missing unit_map: {args.unit_map}",
              file=sys.stderr)
        return 2
    if not args.p1_itp.exists():
        print(f"generate_pet_anchor_posres: missing P1.itp: {args.p1_itp}",
              file=sys.stderr)
        return 2

    unit_map = yaml.safe_load(args.unit_map.read_text(encoding="utf-8"))
    p1_atoms = parse_p1_itp_atoms(args.p1_itp.read_text(encoding="utf-8"))
    if not p1_atoms:
        print(f"generate_pet_anchor_posres: P1.itp has no parseable "
              f"[ atoms ] section: {args.p1_itp}", file=sys.stderr)
        return 2

    try:
        indices = resolve_anchor_atom_indices(
            unit_map=unit_map,
            p1_atoms_index=p1_atoms,
            window=tuple(args.window),  # type: ignore[arg-type]
            direction=args.direction,
        )
    except ValueError as exc:
        print(f"generate_pet_anchor_posres: {exc}", file=sys.stderr)
        return 2

    body = render_posre_itp(
        atom_indices=indices,
        posres_k=args.posres_k,
        case_id=args.case_id,
        direction=args.direction,
        window=tuple(args.window),  # type: ignore[arg-type]
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(body, encoding="utf-8")
    print(f"[done] wrote {len(indices)} POSRES atoms (K={args.posres_k}) → "
          f"{args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
