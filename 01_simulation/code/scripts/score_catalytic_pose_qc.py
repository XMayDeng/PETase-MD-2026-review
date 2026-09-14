#!/usr/bin/env python3
"""Score static PET catalytic-pose geometry before short MD gates.

This is a stricter explanatory QC layer than the original Ser-carbonyl distance
screen. It does not use Vina charges or scores. It only checks whether a
force-field PET pose has a plausible catalytic geometry in the receptor.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import yaml


DEFAULT_ACTIVE_SITE_REGISTRY = (
    "./inputs/active_site_registry.yaml"
)
DEFAULT_PET_RESNAMES = {"PET", "PETL", "PEA", "PEM", "PEB"}
DEFAULT_DONOR_NAMES = {
    "N",
    "ND1",
    "NE2",
    "NE1",
    "NE",
    "NH1",
    "NH2",
    "NZ",
    "OG",
    "OG1",
    "OH",
}


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resid: int
    xyz: tuple[float, float, float]
    element: str

    @property
    def is_heavy(self) -> bool:
        element = self.element.upper()
        name = self.name.upper()
        return element != "H" and not name.startswith("H")

    @property
    def label(self) -> str:
        chain = self.chain or "-"
        return f"{self.resname}{self.resid}:{self.name}:{chain}"


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
        xyz = (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
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
        xyz=xyz,
        element=element,
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


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def angle_degrees(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    ab = tuple(x - y for x, y in zip(a, b))
    cb = tuple(x - y for x, y in zip(c, b))
    denom = math.sqrt(sum(x * x for x in ab)) * math.sqrt(sum(x * x for x in cb))
    if denom == 0:
        raise ValueError("Cannot calculate angle with zero-length vector")
    cos_value = sum(x * y for x, y in zip(ab, cb)) / denom
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_value))))


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


def load_active_site_from_registry(path: Path, protein_id: str) -> dict[str, object]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    proteins = data.get("proteins", {})
    entry = None
    for key in protein_id_lookup_keys(protein_id):
        candidate = proteins.get(key)
        if isinstance(candidate, dict):
            entry = candidate
            break
    if not entry:
        raise KeyError(f"Protein {protein_id} not found in active-site registry: {path}")

    triad = entry.get("catalytic_triad", {})
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
    if ser is None:
        raise KeyError(f"Protein {protein_id} has no catalytic Ser in {path}")
    oxyanion_hole = entry.get("oxyanion_hole") if isinstance(entry.get("oxyanion_hole"), dict) else {}
    return {
        "source": "active_site_registry",
        "protein_id": protein_id,
        "status": entry.get("status"),
        "ser_resid": int(ser),
        "ser_atom_name": entry.get("ser_atom", "OG"),
        "his_resid": int(his) if his is not None else None,
        "acid_resid": int(acid) if acid is not None else None,
        "expected_oxyanion_backbone_n_donors": normalize_expected_backbone_n_donors(
            oxyanion_hole.get("expected_backbone_n_donors", [])
        ),
    }


def is_pet_atom(atom: Atom, pet_chain: str | None, pet_resnames: set[str]) -> bool:
    if pet_chain and atom.chain == pet_chain:
        return True
    resname = atom.resname.upper()
    return resname in pet_resnames or resname.startswith("PET")


def find_atom(
    atoms: Iterable[Atom],
    *,
    resid: int,
    name: str,
    is_pet: bool,
    pet_chain: str | None,
    pet_resnames: set[str],
) -> Atom:
    name_upper = name.upper()
    for atom in atoms:
        if is_pet_atom(atom, pet_chain, pet_resnames) != is_pet:
            continue
        if atom.resid == resid and atom.name.upper() == name_upper:
            return atom
    kind = "PET" if is_pet else "protein"
    raise ValueError(f"Could not find {kind} atom {resid}:{name}")


def optional_find_atom(
    atoms: Iterable[Atom],
    *,
    resid: int,
    name: str,
    is_pet: bool,
    pet_chain: str | None,
    pet_resnames: set[str],
) -> Atom | None:
    try:
        return find_atom(
            atoms,
            resid=resid,
            name=name,
            is_pet=is_pet,
            pet_chain=pet_chain,
            pet_resnames=pet_resnames,
        )
    except ValueError:
        return None


def carbonyl_partner_names(carbonyl_atom_name: str) -> tuple[str, str]:
    name = carbonyl_atom_name.upper()
    if name == "C1":
        return "O2", "O1"
    if name == "C8":
        return "O4", "O3"
    raise ValueError(f"Unsupported PET carbonyl atom name: {carbonyl_atom_name}")


def pair_count_metrics(
    ligand_atoms: list[Atom],
    receptor_atoms: list[Atom],
    cutoffs: tuple[float, ...] = (1.0, 2.0, 4.0),
) -> dict[str, float | int | None]:
    counts = {cutoff: 0 for cutoff in cutoffs}
    min_dist: float | None = None
    for lig_atom in ligand_atoms:
        if not lig_atom.is_heavy:
            continue
        for rec_atom in receptor_atoms:
            if not rec_atom.is_heavy:
                continue
            d = distance(lig_atom.xyz, rec_atom.xyz)
            min_dist = d if min_dist is None else min(min_dist, d)
            for cutoff in cutoffs:
                if d < cutoff:
                    counts[cutoff] += 1
    return {
        "min_heavy_atom_distance_angstrom": round(min_dist, 3) if min_dist is not None else None,
        "heavy_pairs_lt_1a": counts[1.0],
        "heavy_pairs_lt_2a": counts[2.0],
        "heavy_pairs_lt_4a": counts[4.0],
    }


def donor_candidates(
    protein_atoms: list[Atom],
    carbonyl_oxygen: Atom,
    *,
    max_distance: float,
    donor_names: set[str],
) -> list[dict[str, object]]:
    donors = []
    for atom in protein_atoms:
        if not atom.is_heavy:
            continue
        if atom.name.upper() not in donor_names:
            continue
        d = distance(atom.xyz, carbonyl_oxygen.xyz)
        if d <= max_distance:
            donors.append(
                {
                    "resname": atom.resname,
                    "resid": atom.resid,
                    "atom": atom.name,
                    "chain": atom.chain,
                    "distance_angstrom": round(d, 3),
                    "label": atom.label,
                }
            )
    return sorted(donors, key=lambda item: (float(item["distance_angstrom"]), str(item["label"])))


def backbone_n_donors(donors: list[dict[str, object]]) -> list[dict[str, object]]:
    return [donor for donor in donors if str(donor["atom"]).upper() == "N"]


def normalize_expected_backbone_n_donors(items: object) -> list[dict[str, object]]:
    if not items:
        return []
    if not isinstance(items, list):
        raise ValueError("expected_backbone_n_donors must be a list")
    donors: list[dict[str, object]] = []
    for item in items:
        if not isinstance(item, dict):
            raise ValueError("Each expected backbone-N donor must be a mapping")
        if item.get("resid") is None:
            raise ValueError("Expected backbone-N donor is missing resid")
        donors.append(
            {
                "resid": int(item["resid"]),
                "resname": str(item["resname"]).upper() if item.get("resname") else None,
                "atom": str(item.get("atom") or "N").upper(),
                "chain": str(item["chain"]) if item.get("chain") is not None else None,
                "maps_to_positive_control": item.get("maps_to_positive_control"),
            }
        )
    return donors


def expected_donor_label(spec: dict[str, object]) -> str:
    resname = spec.get("resname") or "RES"
    chain = spec.get("chain")
    label = f"{resname}{spec['resid']}:{spec.get('atom') or 'N'}"
    if chain:
        label = f"{label}:{chain}"
    return label


def resolve_expected_donor_atoms(
    protein_atoms: list[Atom],
    expected_donors: list[dict[str, object]],
) -> tuple[list[tuple[dict[str, object], Atom]], list[dict[str, object]]]:
    resolved: list[tuple[dict[str, object], Atom]] = []
    missing: list[dict[str, object]] = []
    for spec in expected_donors:
        resid = int(spec["resid"])
        atom_name = str(spec.get("atom") or "N").upper()
        resname = str(spec["resname"]).upper() if spec.get("resname") else None
        chain = str(spec["chain"]) if spec.get("chain") is not None else None
        matches = [
            atom
            for atom in protein_atoms
            if atom.resid == resid
            and atom.name.upper() == atom_name
            and (resname is None or atom.resname.upper() == resname)
            and (chain is None or atom.chain == chain)
        ]
        if matches:
            resolved.append((spec, matches[0]))
        else:
            missing.append(spec)
    return resolved, missing


def expected_backbone_n_donor_hits(
    resolved_expected_donors: list[tuple[dict[str, object], Atom]],
    carbonyl_oxygen: Atom,
    *,
    max_distance: float,
) -> list[dict[str, object]]:
    hits = []
    for spec, atom in resolved_expected_donors:
        d = distance(atom.xyz, carbonyl_oxygen.xyz)
        if d <= max_distance:
            hits.append(
                {
                    "resname": atom.resname,
                    "resid": atom.resid,
                    "atom": atom.name,
                    "chain": atom.chain,
                    "distance_angstrom": round(d, 3),
                    "label": atom.label,
                    "expected_label": expected_donor_label(spec),
                    "maps_to_positive_control": spec.get("maps_to_positive_control"),
                }
            )
    return sorted(hits, key=lambda item: (float(item["distance_angstrom"]), str(item["label"])))


def nearest_expected_backbone_n_donor(
    resolved_expected_donors: list[tuple[dict[str, object], Atom]],
    carbonyl_oxygen: Atom,
) -> dict[str, object] | None:
    if not resolved_expected_donors:
        return None
    spec, atom = min(
        resolved_expected_donors,
        key=lambda item: distance(item[1].xyz, carbonyl_oxygen.xyz),
    )
    return {
        "resname": atom.resname,
        "resid": atom.resid,
        "atom": atom.name,
        "chain": atom.chain,
        "distance_angstrom": round(distance(atom.xyz, carbonyl_oxygen.xyz), 3),
        "label": atom.label,
        "expected_label": expected_donor_label(spec),
        "maps_to_positive_control": spec.get("maps_to_positive_control"),
    }


def nearest_donor_candidate(
    protein_atoms: list[Atom],
    carbonyl_oxygen: Atom,
    *,
    donor_names: set[str],
) -> dict[str, object] | None:
    donor_atoms = [
        atom
        for atom in protein_atoms
        if atom.is_heavy and atom.name.upper() in donor_names
    ]
    if not donor_atoms:
        return None
    nearest = min(donor_atoms, key=lambda atom: distance(atom.xyz, carbonyl_oxygen.xyz))
    return {
        "resname": nearest.resname,
        "resid": nearest.resid,
        "atom": nearest.name,
        "chain": nearest.chain,
        "distance_angstrom": round(distance(nearest.xyz, carbonyl_oxygen.xyz), 3),
        "label": nearest.label,
    }


def gate_fail_reasons(flags: dict[str, bool]) -> list[str]:
    reasons = []
    if not flags["ser_carbonyl_distance_pass"]:
        reasons.append("ser_carbonyl_distance")
    if not flags["ser_attack_angle_pass"]:
        reasons.append("ser_attack_angle")
    if not flags["oxyanion_donor_pass"]:
        reasons.append("oxyanion_donor")
    if not flags["oxyanion_backbone_n_donor_pass"]:
        reasons.append("oxyanion_backbone_n_donor")
    if not flags["expected_oxyanion_backbone_n_donor_pass"]:
        reasons.append("expected_oxyanion_backbone_n_donor")
    if not flags["clash_pass"]:
        reasons.append("severe_clash")
    if not flags["ser_his_pass"]:
        reasons.append("ser_his_distance")
    return reasons


def min_distance_to_atoms(source: Atom, atoms: list[Atom]) -> dict[str, object] | None:
    if not atoms:
        return None
    nearest = min(atoms, key=lambda atom: distance(source.xyz, atom.xyz))
    return {
        "distance_angstrom": round(distance(source.xyz, nearest.xyz), 3),
        "atom": nearest.name,
        "resname": nearest.resname,
        "resid": nearest.resid,
        "chain": nearest.chain,
        "label": nearest.label,
    }


def atom_records_for_resid(atoms: list[Atom], resid: int) -> list[Atom]:
    return [atom for atom in atoms if atom.resid == resid]


def score_static_pose(
    *,
    atoms: list[Atom],
    ser_resid: int,
    ser_atom_name: str,
    his_resid: int | None,
    acid_resid: int | None,
    pet_chain: str | None,
    pet_resnames: set[str],
    thresholds: dict[str, float | int],
    expected_oxyanion_backbone_n_donors: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    protein_atoms = [atom for atom in atoms if not is_pet_atom(atom, pet_chain, pet_resnames)]
    pet_atoms = [atom for atom in atoms if is_pet_atom(atom, pet_chain, pet_resnames)]
    if not pet_atoms:
        raise ValueError("No PET atoms found in complex")
    expected_donors = normalize_expected_backbone_n_donors(
        expected_oxyanion_backbone_n_donors or []
    )
    resolved_expected_donors, missing_expected_donors = resolve_expected_donor_atoms(
        protein_atoms,
        expected_donors,
    )

    ser_atom = find_atom(
        atoms,
        resid=ser_resid,
        name=ser_atom_name,
        is_pet=False,
        pet_chain=pet_chain,
        pet_resnames=pet_resnames,
    )
    his_atoms = [
        atom
        for atom in protein_atoms
        if his_resid is not None
        and atom.resid == his_resid
        and atom.name.upper() in {"ND1", "NE2"}
    ]
    acid_atoms = [
        atom
        for atom in protein_atoms
        if acid_resid is not None
        and atom.resid == acid_resid
        and atom.name.upper() in {"OD1", "OD2", "OE1", "OE2"}
    ]
    whole_pet_contacts = pair_count_metrics(pet_atoms, protein_atoms)
    candidates = []
    for carbonyl in pet_atoms:
        if not carbonyl.is_heavy or carbonyl.name.upper() not in {"C1", "C8"}:
            continue
        carbonyl_o_name, ester_o_name = carbonyl_partner_names(carbonyl.name)
        carbonyl_oxygen = optional_find_atom(
            atoms,
            resid=carbonyl.resid,
            name=carbonyl_o_name,
            is_pet=True,
            pet_chain=pet_chain,
            pet_resnames=pet_resnames,
        )
        ester_oxygen = optional_find_atom(
            atoms,
            resid=carbonyl.resid,
            name=ester_o_name,
            is_pet=True,
            pet_chain=pet_chain,
            pet_resnames=pet_resnames,
        )
        if carbonyl_oxygen is None:
            continue
        ser_carbonyl_c = distance(ser_atom.xyz, carbonyl.xyz)
        ser_carbonyl_o = distance(ser_atom.xyz, carbonyl_oxygen.xyz)
        attack_angle = angle_degrees(ser_atom.xyz, carbonyl.xyz, carbonyl_oxygen.xyz)
        leaving_angle = (
            angle_degrees(ser_atom.xyz, carbonyl.xyz, ester_oxygen.xyz)
            if ester_oxygen is not None
            else None
        )
        oxy_donors = donor_candidates(
            protein_atoms,
            carbonyl_oxygen,
            max_distance=float(thresholds["max_oxyanion_donor_distance_angstrom"]),
            donor_names=DEFAULT_DONOR_NAMES,
        )
        oxy_backbone_n_donors = backbone_n_donors(oxy_donors)
        nearest_oxy_donor = nearest_donor_candidate(
            protein_atoms,
            carbonyl_oxygen,
            donor_names=DEFAULT_DONOR_NAMES,
        )
        nearest_oxy_backbone_n_donor = nearest_donor_candidate(
            protein_atoms,
            carbonyl_oxygen,
            donor_names={"N"},
        )
        expected_oxy_backbone_hits = expected_backbone_n_donor_hits(
            resolved_expected_donors,
            carbonyl_oxygen,
            max_distance=float(thresholds["max_oxyanion_donor_distance_angstrom"]),
        )
        nearest_expected_oxy_backbone = nearest_expected_backbone_n_donor(
            resolved_expected_donors,
            carbonyl_oxygen,
        )
        unit_contacts = pair_count_metrics(atom_records_for_resid(pet_atoms, carbonyl.resid), protein_atoms)
        his_ser = min_distance_to_atoms(ser_atom, his_atoms)
        his_carbonyl_o = min_distance_to_atoms(carbonyl_oxygen, his_atoms)
        acid_his = None
        if his_atoms and acid_atoms:
            nearest_pair = min(
                ((distance(his.xyz, acid.xyz), his, acid) for his in his_atoms for acid in acid_atoms),
                key=lambda item: item[0],
            )
            acid_his = {
                "distance_angstrom": round(nearest_pair[0], 3),
                "his_atom": nearest_pair[1].label,
                "acid_atom": nearest_pair[2].label,
            }

        distance_pass = ser_carbonyl_c <= float(thresholds["max_ser_carbonyl_c_distance_angstrom"])
        angle_pass = (
            float(thresholds["min_ser_attack_angle_deg"])
            <= attack_angle
            <= float(thresholds["max_ser_attack_angle_deg"])
        )
        oxyanion_pass = len(oxy_donors) >= int(thresholds["min_oxyanion_donors"])
        oxyanion_backbone_pass = len(oxy_backbone_n_donors) >= int(
            thresholds["min_oxyanion_backbone_n_donors"]
        )
        min_expected_backbone_n = int(
            thresholds.get("min_expected_oxyanion_backbone_n_donors", 1)
        )
        expected_donor_pass = (
            not expected_donors
            or min_expected_backbone_n <= 0
            or (
                not missing_expected_donors
                and len(expected_oxy_backbone_hits) >= min_expected_backbone_n
            )
        )
        clash_pass = int(whole_pet_contacts["heavy_pairs_lt_1a"]) <= int(thresholds["max_severe_clashes_lt_1a"])
        his_ser_pass = (
            his_ser is None
            or float(his_ser["distance_angstrom"]) <= float(thresholds["max_ser_his_distance_angstrom"])
        )
        pose_pass = (
            distance_pass
            and angle_pass
            and oxyanion_pass
            and oxyanion_backbone_pass
            and expected_donor_pass
            and clash_pass
            and his_ser_pass
        )
        flags = {
            "ser_carbonyl_distance_pass": distance_pass,
            "ser_attack_angle_pass": angle_pass,
            "oxyanion_donor_pass": oxyanion_pass,
            "oxyanion_backbone_n_donor_pass": oxyanion_backbone_pass,
            "expected_oxyanion_backbone_n_donor_pass": expected_donor_pass,
            "clash_pass": clash_pass,
            "ser_his_pass": his_ser_pass,
            "static_pose_qc_v2_pass": pose_pass,
        }

        candidates.append(
            {
                "carbonyl_unit": carbonyl.resid,
                "carbonyl_atom": carbonyl.name,
                "carbonyl_oxygen_atom": carbonyl_oxygen.name,
                "ester_oxygen_atom": ester_oxygen.name if ester_oxygen else None,
                "ser_og_to_carbonyl_c_angstrom": round(ser_carbonyl_c, 3),
                "ser_og_to_carbonyl_o_angstrom": round(ser_carbonyl_o, 3),
                "ser_og_carbonyl_c_o_angle_deg": round(attack_angle, 1),
                "ser_og_carbonyl_c_ester_o_angle_deg": round(leaving_angle, 1)
                if leaving_angle is not None
                else None,
                "oxyanion_donors": oxy_donors[:10],
                "oxyanion_donor_count": len(oxy_donors),
                "oxyanion_backbone_n_donors": oxy_backbone_n_donors[:10],
                "oxyanion_backbone_n_donor_count": len(oxy_backbone_n_donors),
                "nearest_oxyanion_donor": nearest_oxy_donor,
                "nearest_oxyanion_backbone_n_donor": nearest_oxy_backbone_n_donor,
                "expected_oxyanion_backbone_n_donors": expected_donors,
                "missing_expected_oxyanion_backbone_n_donors": missing_expected_donors,
                "expected_oxyanion_backbone_n_donor_hits": expected_oxy_backbone_hits,
                "expected_oxyanion_backbone_n_donor_hit_count": len(expected_oxy_backbone_hits),
                "nearest_expected_oxyanion_backbone_n_donor": nearest_expected_oxy_backbone,
                "his_to_ser": his_ser,
                "his_to_carbonyl_o": his_carbonyl_o,
                "his_to_acid": acid_his,
                "reactive_unit_contacts": unit_contacts,
                "gate_fail_reasons": gate_fail_reasons(flags),
                "flags": flags,
            }
        )

    if not candidates:
        raise ValueError("No PET carbonyl C1/C8 atoms with carbonyl oxygen partners found")
    candidates = sorted(candidates, key=candidate_sort_key)
    passing = [item for item in candidates if item["flags"]["static_pose_qc_v2_pass"]]
    distance_only = [
        item
        for item in candidates
        if item["flags"]["ser_carbonyl_distance_pass"] and item["flags"]["clash_pass"]
    ]
    verdict = "PASS" if passing else "WARN" if distance_only else "FAIL"
    expected_note = (
        " Expected oxyanion-hole donor annotations are enforced when present."
        if expected_donors
        else ""
    )
    interpretation = (
        "At least one PET carbonyl satisfies distance, attack-angle, oxyanion, clash, and Ser-His gates."
        + expected_note
        if passing
        else "At least one PET carbonyl is close to Ser without severe clashes, but v2 catalytic geometry is incomplete."
        + expected_note
        if distance_only
        else "No PET carbonyl satisfies the minimum Ser-carbonyl distance/clash screen."
        + expected_note
    )
    return {
        "schema_version": 1,
        "verdict": verdict,
        "interpretation": interpretation,
        "active_site": {
            "ser_resid": ser_resid,
            "ser_atom_name": ser_atom_name,
            "his_resid": his_resid,
            "acid_resid": acid_resid,
            "expected_oxyanion_backbone_n_donors": expected_donors,
            "missing_expected_oxyanion_backbone_n_donors": missing_expected_donors,
        },
        "thresholds": thresholds,
        "whole_pet_contacts": whole_pet_contacts,
        "best_candidate": candidates[0],
        "passing_candidates": passing,
        "candidates": candidates,
    }


def candidate_sort_key(candidate: dict[str, object]) -> tuple:
    flags = candidate["flags"]
    nearest_backbone_n = candidate.get("nearest_oxyanion_backbone_n_donor") or {}
    nearest_backbone_n_distance = float(nearest_backbone_n.get("distance_angstrom", 999.0))
    nearest_expected_backbone_n = candidate.get("nearest_expected_oxyanion_backbone_n_donor") or {}
    nearest_expected_backbone_n_distance = float(
        nearest_expected_backbone_n.get("distance_angstrom", 999.0)
    )
    return (
        0 if flags["static_pose_qc_v2_pass"] else 1,
        0 if flags["ser_carbonyl_distance_pass"] else 1,
        0 if flags["ser_attack_angle_pass"] else 1,
        0 if flags["expected_oxyanion_backbone_n_donor_pass"] else 1,
        nearest_expected_backbone_n_distance,
        -int(candidate["expected_oxyanion_backbone_n_donor_hit_count"]),
        0 if flags["oxyanion_backbone_n_donor_pass"] else 1,
        0 if flags["oxyanion_donor_pass"] else 1,
        nearest_backbone_n_distance,
        float(candidate["ser_og_to_carbonyl_c_angstrom"]),
        abs(float(candidate["ser_og_carbonyl_c_o_angle_deg"]) - 105.0),
        -int(candidate["oxyanion_backbone_n_donor_count"]),
        -int(candidate["oxyanion_donor_count"]),
    )


def default_thresholds(args: argparse.Namespace) -> dict[str, float | int]:
    return {
        "max_ser_carbonyl_c_distance_angstrom": args.max_ser_carbonyl_c_distance,
        "min_ser_attack_angle_deg": args.min_ser_attack_angle,
        "max_ser_attack_angle_deg": args.max_ser_attack_angle,
        "max_oxyanion_donor_distance_angstrom": args.max_oxyanion_donor_distance,
        "min_oxyanion_donors": args.min_oxyanion_donors,
        "min_oxyanion_backbone_n_donors": args.min_oxyanion_backbone_n_donors,
        "min_expected_oxyanion_backbone_n_donors": args.min_expected_oxyanion_backbone_n_donors,
        "max_ser_his_distance_angstrom": args.max_ser_his_distance,
        "max_severe_clashes_lt_1a": args.max_severe_clashes,
    }


def write_csv(path: Path, result: dict[str, object]) -> None:
    fieldnames = [
        "rank",
        "static_pose_qc_v2_pass",
        "carbonyl_unit",
        "carbonyl_atom",
        "ser_og_to_carbonyl_c_angstrom",
        "ser_og_carbonyl_c_o_angle_deg",
        "oxyanion_donor_count",
        "oxyanion_backbone_n_donor_count",
        "expected_oxyanion_backbone_n_donor_hit_count",
        "nearest_oxyanion_donor",
        "nearest_oxyanion_backbone_n_donor",
        "nearest_expected_oxyanion_backbone_n_donor",
        "nearest_oxyanion_donor_distance_angstrom",
        "nearest_oxyanion_backbone_n_donor_distance_angstrom",
        "nearest_expected_oxyanion_backbone_n_donor_distance_angstrom",
        "gate_fail_reasons",
        "whole_pet_heavy_pairs_lt_1a",
        "whole_pet_heavy_pairs_lt_4a",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for rank, candidate in enumerate(result["candidates"], start=1):
            donors = candidate["oxyanion_donors"]
            backbone_donors = candidate["oxyanion_backbone_n_donors"]
            nearest_donor = donors[0]["label"] if donors else ""
            nearest_backbone = backbone_donors[0]["label"] if backbone_donors else ""
            expected_hits = candidate["expected_oxyanion_backbone_n_donor_hits"]
            nearest_expected = expected_hits[0]["label"] if expected_hits else ""
            nearest_any_donor = candidate.get("nearest_oxyanion_donor") or {}
            nearest_any_backbone = candidate.get("nearest_oxyanion_backbone_n_donor") or {}
            nearest_any_expected = candidate.get("nearest_expected_oxyanion_backbone_n_donor") or {}
            writer.writerow(
                {
                    "rank": rank,
                    "static_pose_qc_v2_pass": candidate["flags"]["static_pose_qc_v2_pass"],
                    "carbonyl_unit": candidate["carbonyl_unit"],
                    "carbonyl_atom": candidate["carbonyl_atom"],
                    "ser_og_to_carbonyl_c_angstrom": candidate["ser_og_to_carbonyl_c_angstrom"],
                    "ser_og_carbonyl_c_o_angle_deg": candidate["ser_og_carbonyl_c_o_angle_deg"],
                    "oxyanion_donor_count": candidate["oxyanion_donor_count"],
                    "oxyanion_backbone_n_donor_count": candidate[
                        "oxyanion_backbone_n_donor_count"
                    ],
                    "expected_oxyanion_backbone_n_donor_hit_count": candidate[
                        "expected_oxyanion_backbone_n_donor_hit_count"
                    ],
                    "nearest_oxyanion_donor": nearest_donor,
                    "nearest_oxyanion_backbone_n_donor": nearest_backbone,
                    "nearest_expected_oxyanion_backbone_n_donor": nearest_expected,
                    "nearest_oxyanion_donor_distance_angstrom": nearest_any_donor.get("distance_angstrom"),
                    "nearest_oxyanion_backbone_n_donor_distance_angstrom": nearest_any_backbone.get("distance_angstrom"),
                    "nearest_expected_oxyanion_backbone_n_donor_distance_angstrom": nearest_any_expected.get("distance_angstrom"),
                    "gate_fail_reasons": ";".join(candidate.get("gate_fail_reasons", [])),
                    "whole_pet_heavy_pairs_lt_1a": result["whole_pet_contacts"]["heavy_pairs_lt_1a"],
                    "whole_pet_heavy_pairs_lt_4a": result["whole_pet_contacts"]["heavy_pairs_lt_4a"],
                }
            )


def write_markdown(path: Path, result: dict[str, object]) -> None:
    best = result["best_candidate"]
    donors = best["oxyanion_donors"]
    backbone_donors = best["oxyanion_backbone_n_donors"]
    expected_hits = best["expected_oxyanion_backbone_n_donor_hits"]
    nearest_donor = donors[0]["label"] if donors else "none"
    nearest_backbone = backbone_donors[0]["label"] if backbone_donors else "none"
    nearest_expected = expected_hits[0]["label"] if expected_hits else "none"
    nearest_any_donor = best.get("nearest_oxyanion_donor") or {}
    nearest_any_backbone = best.get("nearest_oxyanion_backbone_n_donor") or {}
    nearest_any_donor_text = (
        f"{nearest_any_donor.get('label')} at {nearest_any_donor.get('distance_angstrom')} A"
        if nearest_any_donor
        else "none"
    )
    nearest_any_backbone_text = (
        f"{nearest_any_backbone.get('label')} at {nearest_any_backbone.get('distance_angstrom')} A"
        if nearest_any_backbone
        else "none"
    )
    nearest_any_expected = best.get("nearest_expected_oxyanion_backbone_n_donor") or {}
    nearest_any_expected_text = (
        f"{nearest_any_expected.get('label')} at {nearest_any_expected.get('distance_angstrom')} A"
        if nearest_any_expected
        else "none"
    )
    lines = [
        "# Catalytic Pose QC v2",
        "",
        f"- Verdict: `{result['verdict']}`",
        f"- Interpretation: {result['interpretation']}",
        f"- Best carbonyl: unit `{best['carbonyl_unit']}` atom `{best['carbonyl_atom']}`",
        f"- Ser-carbonyl C distance: `{best['ser_og_to_carbonyl_c_angstrom']} A`",
        f"- Ser-CarbonylC-O angle: `{best['ser_og_carbonyl_c_o_angle_deg']} deg`",
        f"- Oxyanion donors within threshold: `{best['oxyanion_donor_count']}`; nearest `{nearest_donor}`",
        f"- Backbone N oxyanion donors within threshold: `{best['oxyanion_backbone_n_donor_count']}`; nearest `{nearest_backbone}`",
        f"- Expected backbone N donor hits: `{best['expected_oxyanion_backbone_n_donor_hit_count']}`; nearest hit `{nearest_expected}`",
        f"- Nearest donor overall: `{nearest_any_donor_text}`",
        f"- Nearest backbone N overall: `{nearest_any_backbone_text}`",
        f"- Nearest expected backbone N overall: `{nearest_any_expected_text}`",
        f"- Gate fail reasons: `{', '.join(best.get('gate_fail_reasons', [])) or 'none'}`",
        f"- Whole PET severe `<1 A` pairs: `{result['whole_pet_contacts']['heavy_pairs_lt_1a']}`",
        "",
        "## Candidate Carbonyls",
        "",
        "| rank | pass | unit | atom | Ser-C A | attack angle | oxy N donors | expected N hits | nearest expected N A | nearest backbone N A | <1A clashes |",
        "|---:|:---:|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for rank, candidate in enumerate(result["candidates"], start=1):
        nearest_candidate_backbone = candidate.get("nearest_oxyanion_backbone_n_donor") or {}
        nearest_candidate_expected = candidate.get("nearest_expected_oxyanion_backbone_n_donor") or {}
        lines.append(
            f"| {rank} | {'yes' if candidate['flags']['static_pose_qc_v2_pass'] else 'no'} | "
            f"{candidate['carbonyl_unit']} | {candidate['carbonyl_atom']} | "
            f"{candidate['ser_og_to_carbonyl_c_angstrom']} | "
            f"{candidate['ser_og_carbonyl_c_o_angle_deg']} | "
            f"{candidate['oxyanion_backbone_n_donor_count']} | "
            f"{candidate['expected_oxyanion_backbone_n_donor_hit_count']} | "
            f"{nearest_candidate_expected.get('distance_angstrom', 'NA')} | "
            f"{nearest_candidate_backbone.get('distance_angstrom', 'NA')} | "
            f"{result['whole_pet_contacts']['heavy_pairs_lt_1a']} |"
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- This QC uses PET coordinates and protein geometry only.",
            "- Vina charges, atom types, and absolute scores are not MD force-field inputs.",
            "- A `PASS` here is still not production approval; it only means the pose is worth EM/NVT/NPT and release-gate testing.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Score static PET catalytic-pose geometry.")
    parser.add_argument("--complex-pdb", required=True, type=Path)
    parser.add_argument("--protein-id")
    parser.add_argument("--active-site-registry", default=DEFAULT_ACTIVE_SITE_REGISTRY, type=Path)
    parser.add_argument("--ser-resid", type=int)
    parser.add_argument("--ser-atom-name", default=None)
    parser.add_argument("--his-resid", type=int)
    parser.add_argument("--acid-resid", type=int)
    parser.add_argument("--pet-chain", default="B")
    parser.add_argument("--pet-resname", action="append")
    parser.add_argument("--max-ser-carbonyl-c-distance", type=float, default=5.0)
    parser.add_argument("--min-ser-attack-angle", type=float, default=80.0)
    parser.add_argument("--max-ser-attack-angle", type=float, default=130.0)
    parser.add_argument("--max-oxyanion-donor-distance", type=float, default=3.5)
    parser.add_argument("--min-oxyanion-donors", type=int, default=1)
    parser.add_argument("--min-oxyanion-backbone-n-donors", type=int, default=1)
    parser.add_argument("--min-expected-oxyanion-backbone-n-donors", type=int, default=1)
    parser.add_argument("--max-ser-his-distance", type=float, default=4.0)
    parser.add_argument("--max-severe-clashes", type=int, default=0)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    root = repo_root()
    complex_pdb = resolve_repo_path(args.complex_pdb, root)
    out_dir = resolve_repo_path(args.out_dir, root)
    if out_dir.exists() and not args.force:
        raise FileExistsError(f"Output directory exists: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    site = {}
    if args.protein_id:
        site = load_active_site_from_registry(
            resolve_repo_path(args.active_site_registry, root),
            args.protein_id,
        )
    ser_resid = args.ser_resid if args.ser_resid is not None else site.get("ser_resid")
    if ser_resid is None:
        raise ValueError("--ser-resid or --protein-id is required")
    ser_atom_name = args.ser_atom_name or site.get("ser_atom_name") or "OG"
    his_resid = args.his_resid if args.his_resid is not None else site.get("his_resid")
    acid_resid = args.acid_resid if args.acid_resid is not None else site.get("acid_resid")
    pet_resnames = set(DEFAULT_PET_RESNAMES)
    if args.pet_resname:
        pet_resnames.update(item.upper() for item in args.pet_resname)

    result = score_static_pose(
        atoms=parse_pdb_atoms(complex_pdb),
        ser_resid=int(ser_resid),
        ser_atom_name=str(ser_atom_name),
        his_resid=int(his_resid) if his_resid is not None else None,
        acid_resid=int(acid_resid) if acid_resid is not None else None,
        pet_chain=args.pet_chain,
        pet_resnames=pet_resnames,
        thresholds=default_thresholds(args),
        expected_oxyanion_backbone_n_donors=site.get(
            "expected_oxyanion_backbone_n_donors",
            [],
        ),
    )
    result["inputs"] = {
        "complex_pdb": str(complex_pdb),
        "protein_id": args.protein_id,
        "active_site_source": site.get("source") if site else "cli",
    }
    (out_dir / "catalytic_pose_qc_v2.json").write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    (out_dir / "catalytic_pose_qc_v2.yaml").write_text(
        yaml.safe_dump(result, sort_keys=False),
        encoding="utf-8",
    )
    write_csv(out_dir / "catalytic_pose_qc_v2.csv", result)
    write_markdown(out_dir / "catalytic_pose_qc_v2.md", result)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
