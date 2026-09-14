#!/usr/bin/env python3
"""
Build the GROMACS system_build directory from protein/PET inputs and mdp/parameterization.

What this script does (pipeline):
  1) pdb2gmx: generate protein topology (topol.top), coords (protein.gro), posre (posre.itp)
  2) editconf: create simulation box (box.gro)
  3) insert-molecules: insert PET into box (conf_pet.gro)
  4) solvate: add water (box_sol.gro) and update topol.top
  5) grompp + genion: add ions (box_sol_ion.gro, ions.tpr) and update topol.top
  6) make_ndx: create index.ndx and ensure PETL / Protein_LIG / Water_and_ions groups

Default inputs (relative to generalGromacs_PL):
  - Protein_data/00057.pdb
  - inputs/forcefields/charmmgui_gmx/step3_input.gro
  - inputs/forcefields/charmmgui_gmx/toppar/P1.itp
  - inputs/forcefields/charmm36-jul2021.ff
  - inputs/forcefields/charmmgui_gmx/toppar/forcefield.itp (PET add-on)
  - configurations/templates/md/ions.mdp

Outputs in code/system_build/:
  protein.gro, topol.top, posre.itp
  box.gro, conf_pet.gro, box_sol.gro, box_sol_ion.gro
  ions.tpr, index.ndx
  PET.gro, P1.itp, pet_ff_add.itp

Notes:
  - By default, the forcefield directory is NOT copied into system_build.
    GMXLIB is set to include parameterization so GROMACS can resolve includes.
    Use --copy-ff if you want a local forcefield copy.
  - auto_box supports long ligands: off/suggest/apply.
  - Optional post-check runs check_system_build.py automatically.
  - This script does NOT run EM/NVT/NPT/MD; use run_pipeline.py for that.
  - Script config (optional): configurations/templates/md/build_system.params.mdp
    CLI args override config values.
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path


def find_gmx() -> str:
    for name in ("gmx_mpi", "gmx"):
        path = shutil.which(name)
        if path:
            return path
    raise SystemExit("Cannot find gmx_mpi/gmx in PATH. Please source GMXRC.")


def supports_color() -> bool:
    return sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


USE_COLOR = supports_color()


def highlight(tag: str, msg: str, color: str) -> None:
    prefix = f"[{tag}]"
    if USE_COLOR:
        prefix = f"\033[{color}m{prefix}\033[0m"
    print(f"{prefix} {msg}")


def run(cmd: list[str], cwd: Path | None = None, stdin: str | None = None, env: dict | None = None) -> None:
    cmd_str = " ".join(cmd)
    if cwd:
        highlight("CMD", f"{cmd_str}  (cwd={cwd})", "1;36")
    else:
        highlight("CMD", cmd_str, "1;36")
    subprocess.run(
        cmd,
        check=True,
        cwd=cwd,
        input=stdin,
        text=True,
        env=env,
    )


def parse_bool(val: str) -> bool:
    return val.strip().lower() in ("1", "true", "yes", "y", "on")


def parse_config(path: Path) -> dict[str, str]:
    cfg: dict[str, str] = {}
    if not path.exists():
        return cfg
    for raw in path.read_text().splitlines():
        line = raw.split(";", 1)[0].split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, val = line.split("=", 1)
        cfg[key.strip().lower()] = val.strip()
    return cfg


def resolve_path(val: str, root: Path) -> str:
    p = Path(val)
    if p.is_absolute():
        return str(p)
    return str(root / p)


def parse_moleculetype(itp_path: Path) -> str:
    in_section = False
    for line in itp_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith("["):
            in_section = stripped.lower().startswith("[ moleculetype")
            continue
        if in_section:
            return stripped.split()[0]
    raise SystemExit(f"Could not find [ moleculetype ] in {itp_path}")


def collect_dihedraltypes(ffbonded: Path) -> set[tuple[str, str, str, str, str]]:
    dihs: set[tuple[str, str, str, str, str]] = set()
    in_dih = False
    for line in ffbonded.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith("["):
            in_dih = stripped.lower().startswith("[ dihedraltypes")
            continue
        if not in_dih:
            continue
        parts = stripped.split()
        if len(parts) < 5:
            continue
        key = tuple(parts[:5])
        dihs.add(key)
    return dihs


def canonical_pair(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a <= b else (b, a)


def param_key(section: str, parts: list[str]) -> tuple[str, ...] | None:
    if section == "atomtypes":
        return (parts[0],) if len(parts) >= 1 else None
    if section == "bondtypes":
        return (*canonical_pair(parts[0], parts[1]), parts[2]) if len(parts) >= 3 else None
    if section == "pairtypes":
        return (*canonical_pair(parts[0], parts[1]), parts[2]) if len(parts) >= 3 else None
    if section == "angletypes":
        if len(parts) < 4:
            return None
        i, j, k = parts[0], parts[1], parts[2]
        a, c = sorted((i, k))
        return (a, j, c, parts[3])
    if section == "dihedraltypes":
        return tuple(parts[:5]) if len(parts) >= 5 else None
    return None


def collect_param_keys_from_file(path: Path) -> dict[str, set[tuple[str, ...]]]:
    keys: dict[str, set[tuple[str, ...]]] = {
        "atomtypes": set(),
        "bondtypes": set(),
        "pairtypes": set(),
        "angletypes": set(),
        "dihedraltypes": set(),
    }
    if not path.exists():
        return keys

    section = ""
    for raw in path.read_text().splitlines():
        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            continue
        if stripped.startswith("["):
            section = stripped.strip("[]").strip().lower()
            continue
        if section not in keys:
            continue
        parts = stripped.split()
        key = param_key(section, parts)
        if key is not None:
            keys[section].add(key)
    return keys


def merge_param_keys(*key_maps: dict[str, set[tuple[str, ...]]]) -> dict[str, set[tuple[str, ...]]]:
    merged: dict[str, set[tuple[str, ...]]] = {
        "atomtypes": set(),
        "bondtypes": set(),
        "pairtypes": set(),
        "angletypes": set(),
        "dihedraltypes": set(),
    }
    for km in key_maps:
        for section in merged:
            merged[section].update(km.get(section, set()))
    return merged


def write_pet_ff_add(src: Path, dst: Path, ffbonded: Path, ffnonbonded: Path) -> None:
    existing = merge_param_keys(
        collect_param_keys_from_file(ffbonded),
        collect_param_keys_from_file(ffnonbonded),
    )
    lines = src.read_text().splitlines()
    out = []
    skip_defaults = False
    section = ""
    seen_local: dict[str, set[tuple[str, ...]]] = {
        "atomtypes": set(),
        "bondtypes": set(),
        "pairtypes": set(),
        "angletypes": set(),
        "dihedraltypes": set(),
    }
    for line in lines:
        stripped = line.strip()
        lower = stripped.lower()
        if lower.startswith("[") and "defaults" in lower:
            skip_defaults = True
            continue
        if skip_defaults and lower.startswith("[") and "defaults" not in lower:
            skip_defaults = False
        if skip_defaults:
            continue

        if lower.startswith("["):
            section = stripped.strip("[]").strip().lower()
            out.append(line)
            continue

        if section in seen_local and stripped and not stripped.startswith(";"):
            parts = stripped.split()
            key = param_key(section, parts)
            if key is not None:
                if key in existing[section] or key in seen_local[section]:
                    continue
                seen_local[section].add(key)
        out.append(line)
    dst.write_text("\n".join(out) + "\n")


def read_gro_coords(path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    with path.open() as f:
        next(f)  # title
        natoms = int(next(f).strip())
        for _ in range(natoms):
            line = next(f)
            coords.append((float(line[20:28]), float(line[28:36]), float(line[36:44])))
    return coords


def coord_extent(coords: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    if not coords:
        return (0.0, 0.0, 0.0)
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def dist3(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def suggest_box_distance(
    protein_gro: Path,
    pet_gro: Path,
    current_box_distance: float,
    pbc_buffer: float,
) -> dict[str, float]:
    prot = read_gro_coords(protein_gro)
    pet = read_gro_coords(pet_gro)

    prot_ex = coord_extent(prot)
    pet_ex = coord_extent(pet)
    pet_span_max = max(pet_ex)
    pet_e2e = dist3(pet[0], pet[-1]) if len(pet) >= 2 else pet_span_max
    target_min_box = max(pet_span_max, pet_e2e) + 2.0 * pbc_buffer
    required_d = max(0.0, 0.5 * (target_min_box - min(prot_ex)))
    suggested_d = max(current_box_distance, required_d)

    return {
        "protein_x": prot_ex[0],
        "protein_y": prot_ex[1],
        "protein_z": prot_ex[2],
        "pet_x": pet_ex[0],
        "pet_y": pet_ex[1],
        "pet_z": pet_ex[2],
        "pet_e2e": pet_e2e,
        "target_min_box": target_min_box,
        "required_d": required_d,
        "suggested_d": suggested_d,
    }


def patch_topol(
    topol_path: Path,
    itp_name: str,
    molname: str,
    ff_add_name: str | None = None,
) -> None:
    lines = topol_path.read_text().splitlines()
    posres_define = [
        "#ifndef POSRES_FC_BB",
        "#define POSRES_FC_BB 1000",
        "#endif",
    ]
    if ff_add_name:
        for i, line in enumerate(lines):
            if "forcefield.itp" in line and "include" in line:
                add_line = f'#include "{ff_add_name}"'
                if add_line not in lines:
                    lines.insert(i + 1, add_line)
                break

    include_line = f'#include "{itp_name}"'
    if any("POSRES_FC_BB" in line for line in lines) is False:
        insert_idx = 0
        for i, line in enumerate(lines):
            if line.strip().startswith("[ system ]"):
                insert_idx = i
                break
        lines[insert_idx:insert_idx] = posres_define
    if include_line not in lines:
        insert_idx = 0
        for i, line in enumerate(lines):
            if line.strip().startswith("[ system ]"):
                insert_idx = i
                break
        lines.insert(insert_idx, include_line)

    # Ensure molecule exists in [ molecules ]
    mol_section = None
    for i, line in enumerate(lines):
        if line.strip().lower().startswith("[ molecules ]"):
            mol_section = i
            break
    if mol_section is None:
        lines.append("")
        lines.append("[ molecules ]")
        mol_section = len(lines) - 1

    found = False
    for i in range(mol_section + 1, len(lines)):
        if lines[i].strip().startswith("["):
            break
        if not lines[i].strip() or lines[i].strip().startswith(";"):
            continue
        parts = lines[i].split()
        if parts and parts[0] == molname:
            lines[i] = f"{molname:<8} 1"
            found = True
            break
    if not found:
        insert_idx = len(lines)
        for i in range(mol_section + 1, len(lines)):
            if lines[i].strip().startswith("["):
                insert_idx = i
                break
        lines.insert(insert_idx, f"{molname:<8} 1")

    topol_path.write_text("\n".join(lines) + "\n")


def read_gro_resnames(gro_path: Path) -> dict[str, list[int]]:
    res_atoms: dict[str, list[int]] = {}
    with gro_path.open() as f:
        next(f)  # title
        natoms = int(next(f).strip())
        for _ in range(natoms):
            line = next(f)
            resname = line[5:10].strip()
            atom_idx = int(line[15:20])
            res_atoms.setdefault(resname, []).append(atom_idx)
    return res_atoms


def read_index(index_path: Path) -> OrderedDict[str, list[int]]:
    groups: OrderedDict[str, list[int]] = OrderedDict()
    current = None
    for line in index_path.read_text().splitlines():
        if line.strip().startswith("["):
            name = line.strip().strip("[]").strip()
            current = name
            groups[current] = []
            continue
        if current is None:
            continue
        if not line.strip():
            continue
        groups[current].extend(int(x) for x in line.split())
    return groups


def write_index(index_path: Path, groups: OrderedDict[str, list[int]]) -> None:
    lines = []
    for name, atoms in groups.items():
        lines.append(f"[ {name} ]")
        for i in range(0, len(atoms), 15):
            chunk = atoms[i : i + 15]
            lines.append(" ".join(f"{a:5d}" for a in chunk))
        lines.append("")
    index_path.write_text("\n".join(lines))


def ensure_groups(index_path: Path, gro_path: Path) -> None:
    groups = read_index(index_path)
    res_atoms = read_gro_resnames(gro_path)

    def add_group(name: str, atoms: list[int]) -> None:
        if name in groups:
            return
        groups[name] = sorted(set(atoms))

    # PETL group
    if "PETL" not in groups and "PETL" in res_atoms:
        add_group("PETL", res_atoms["PETL"])
    if "r PETL" in groups and "PETL" not in groups:
        add_group("PETL", groups["r PETL"])

    # Protein_LIG / Protein_PETL
    if "Protein" in groups and "PETL" in groups:
        union = sorted(set(groups["Protein"]) | set(groups["PETL"]))
        add_group("Protein_LIG", union)
        add_group("Protein_PETL", union)

    # Water_and_ions
    water_atoms = []
    for water_name in ("Water", "SOL", "WAT", "HOH"):
        if water_name in groups:
            water_atoms = groups[water_name]
            break
    ion_atoms = []
    for ion_name in ("SOD", "CLA", "NA", "CL"):
        if ion_name in groups:
            ion_atoms.extend(groups[ion_name])
    if water_atoms:
        add_group("Water_and_ions", sorted(set(water_atoms + ion_atoms)))

    write_index(index_path, groups)


def find_default_pet_gro(root: Path) -> Path:
    candidates = [
        root / "inputs/forcefields/charmmgui_gmx/step3_input.gro",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise SystemExit("Cannot find PET .gro; provide --pet-gro explicitly.")


def find_default_pet_itp(root: Path) -> Path:
    candidates = [
        root / "inputs/forcefields/charmmgui_gmx/toppar/P1.itp",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise SystemExit("Cannot find PET .itp; provide --pet-itp explicitly.")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--config", help="Config file (mdp-style key=value)")
    p.add_argument("--gmx-bin", help="GROMACS binary (gmx_mpi, gmx, or absolute path)")
    p.add_argument("--protein", help="Protein PDB/GRO path")
    p.add_argument("--pet-gro", help="PET .gro path")
    p.add_argument("--pet-itp", help="PET .itp path")
    p.add_argument("--ff-dir", help="Forcefield directory (charmm36-jul2021.ff)")
    p.add_argument("--ff-itp", help="Forcefield .itp (CHARMM-GUI forcefield.itp)")
    p.add_argument(
        "--copy-ff",
        action="store_true",
        default=None,
        help="Copy forcefield into output dir (default: use parameterization path)",
    )
    p.add_argument("--mdp-dir", help="MDP directory (for ions.mdp)")
    p.add_argument("--outdir", help="Output directory (system_build)")
    p.add_argument("--salt", type=float, default=None, help="Salt concentration (M)")
    p.add_argument("--water", default=None, help="Water model for pdb2gmx (default tip3p)")
    p.add_argument("--water-coord", default=None, help="Water coord file for solvate")
    p.add_argument("--box-type", default=None, help="Box type for editconf")
    p.add_argument("--box-distance", type=float, default=None, help="Box distance in nm")
    p.add_argument(
        "--auto-box",
        choices=("off", "suggest", "apply"),
        default=None,
        help="Auto-tune box_distance for long ligands: off/suggest/apply",
    )
    p.add_argument(
        "--pbc-buffer",
        type=float,
        default=None,
        help="Target ligand-image buffer (nm) used by auto-box logic",
    )
    p.add_argument(
        "--insert-seed",
        type=int,
        default=None,
        help="Random seed for insert-molecules (omit for random)",
    )
    p.add_argument(
        "--genion-seed",
        type=int,
        default=None,
        help="Random seed for genion (omit for random)",
    )
    p.add_argument("--genion-group", default=None, help="Group name to replace in genion")
    chk = p.add_mutually_exclusive_group()
    chk.add_argument(
        "--post-check",
        dest="post_check",
        action="store_true",
        default=None,
        help="Run check_system_build.py after building",
    )
    chk.add_argument(
        "--no-post-check",
        dest="post_check",
        action="store_false",
        default=None,
        help="Skip post-build check_system_build.py",
    )
    strict = p.add_mutually_exclusive_group()
    strict.add_argument(
        "--check-strict",
        dest="check_strict",
        action="store_true",
        default=None,
        help="Fail post-check on WARN/FAIL (passes --strict)",
    )
    strict.add_argument(
        "--no-check-strict",
        dest="check_strict",
        action="store_false",
        default=None,
        help="Post-check fails only on FAIL",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve()
    root = script_dir.parents[2]  # packaged simulation root

    default_cfg = root / "configurations/templates/md/build_system.params.mdp"
    cfg_path = Path(args.config) if args.config else default_cfg
    cfg = parse_config(cfg_path)

    def pick(attr: str, key: str, cast=None, is_path: bool = False) -> None:
        if getattr(args, attr) is not None:
            return
        if key not in cfg:
            return
        val = cfg[key]
        if is_path:
            val = resolve_path(val, root)
        if cast:
            setattr(args, attr, cast(val))
        else:
            setattr(args, attr, val)

    pick("gmx_bin", "gmx_bin")
    pet_choice = cfg.get("pet_choice", "").strip()
    pick("protein", "protein", is_path=True)
    pick("pet_gro", "pet_gro", is_path=True)
    pick("pet_itp", "pet_itp", is_path=True)
    pick("ff_dir", "ff_dir", is_path=True)
    pick("ff_itp", "ff_itp", is_path=True)
    pick("mdp_dir", "mdp_dir", is_path=True)
    pick("outdir", "outdir", is_path=True)
    pick("salt", "salt", float)
    pick("water", "water")
    pick("water_coord", "water_coord")
    pick("box_type", "box_type")
    pick("box_distance", "box_distance", float)
    pick("pbc_buffer", "pbc_buffer", float)
    pick("insert_seed", "insert_seed", int)
    pick("genion_seed", "genion_seed", int)
    pick("genion_group", "genion_group")
    if args.copy_ff is None and "copy_ff" in cfg:
        args.copy_ff = parse_bool(cfg["copy_ff"])
    if args.post_check is None and "post_check" in cfg:
        args.post_check = parse_bool(cfg["post_check"])
    if args.check_strict is None and "check_strict" in cfg:
        args.check_strict = parse_bool(cfg["check_strict"])
    if args.auto_box is None and "auto_box" in cfg:
        args.auto_box = cfg["auto_box"].strip().lower()

    if args.salt is None:
        args.salt = 0.1
    if args.water is None:
        args.water = "tip3p"
    if args.water_coord is None:
        args.water_coord = "spc216.gro"
    if args.box_type is None:
        args.box_type = "triclinic"
    if args.box_distance is None:
        args.box_distance = 1.0
    if args.pbc_buffer is None:
        args.pbc_buffer = 1.2
    if args.genion_group is None:
        args.genion_group = "SOL"
    if args.copy_ff is None:
        args.copy_ff = False
    if args.post_check is None:
        args.post_check = True
    if args.check_strict is None:
        args.check_strict = False
    if args.auto_box is None:
        args.auto_box = "suggest"

    if args.auto_box in ("yes", "true", "on", "1"):
        args.auto_box = "apply"
    elif args.auto_box in ("no", "false", "off", "0"):
        args.auto_box = "off"
    elif args.auto_box not in ("off", "suggest", "apply"):
        raise SystemExit("auto_box must be one of: off, suggest, apply")

    protein = Path(args.protein) if args.protein else root / "Protein_data/00057.pdb"
    if pet_choice:
        if args.pet_gro is None:
            args.pet_gro = root / f"PET_data/{pet_choice}/gromacs/step3_input.gro"
        if args.pet_itp is None:
            args.pet_itp = root / f"PET_data/{pet_choice}/gromacs/toppar/P1.itp"
        if args.ff_itp is None:
            args.ff_itp = root / f"PET_data/{pet_choice}/gromacs/toppar/forcefield.itp"
    pet_gro = Path(args.pet_gro) if args.pet_gro else find_default_pet_gro(root)
    pet_itp = Path(args.pet_itp) if args.pet_itp else find_default_pet_itp(root)
    ff_dir = Path(args.ff_dir) if args.ff_dir else root / "inputs/forcefields/charmm36-jul2021.ff"
    ff_itp = Path(args.ff_itp) if args.ff_itp else root / "inputs/forcefields/charmmgui_gmx/toppar/forcefield.itp"
    mdp_dir = Path(args.mdp_dir) if args.mdp_dir else root / "configurations/templates/md"
    outdir = Path(args.outdir) if args.outdir else root / "code/system_build"

    ions_mdp = mdp_dir / "ions.mdp"

    if not protein.exists():
        raise SystemExit(f"Protein file not found: {protein}")
    if not pet_gro.exists():
        raise SystemExit(f"PET gro not found: {pet_gro}")
    if not pet_itp.exists():
        raise SystemExit(f"PET itp not found: {pet_itp}")
    if not ff_dir.exists():
        raise SystemExit(f"Forcefield dir not found: {ff_dir}")
    if not ff_itp.exists():
        raise SystemExit(f"Forcefield itp not found: {ff_itp}")
    if not ions_mdp.exists():
        raise SystemExit(f"ions.mdp not found: {ions_mdp}")

    outdir.mkdir(parents=True, exist_ok=True)
    if args.gmx_bin:
        if Path(args.gmx_bin).is_absolute():
            gmx_path = Path(args.gmx_bin)
            if not gmx_path.exists():
                raise SystemExit(f"gmx_bin not found: {gmx_path}")
            gmx = str(gmx_path)
        else:
            gmx_path = shutil.which(args.gmx_bin)
            if not gmx_path:
                raise SystemExit(f"gmx_bin not found in PATH: {args.gmx_bin}")
            gmx = gmx_path
    else:
        gmx = find_gmx()

    highlight("INFO", "Starting system build...", "1;34")
    highlight("INFO", f"Output dir: {outdir}", "1;34")
    highlight("INFO", f"auto_box: {args.auto_box} (pbc_buffer={args.pbc_buffer} nm)", "1;34")
    if args.insert_seed is not None:
        highlight("INFO", f"insert_seed: {args.insert_seed}", "1;34")
    if args.genion_seed is not None:
        highlight("INFO", f"genion_seed: {args.genion_seed}", "1;34")

    # Forcefield location (default: parameterization path)
    ff_src = ff_dir
    if args.copy_ff:
        highlight("INFO", f"Copying forcefield into {outdir}", "1;34")
        ff_dst = outdir / ff_dir.name
        if not ff_dst.exists():
            shutil.copytree(ff_dir, ff_dst)
        ff_src = ff_dst

    # PET itp + add-on forcefield for PET
    ff_add_dst = outdir / "pet_ff_add.itp"
    ffbonded = ff_src / "ffbonded.itp"
    ffnonbonded = ff_src / "ffnonbonded.itp"
    write_pet_ff_add(ff_itp, ff_add_dst, ffbonded, ffnonbonded)
    itp_dst = outdir / pet_itp.name
    if not itp_dst.exists():
        shutil.copy2(pet_itp, itp_dst)
    shutil.copy2(pet_gro, outdir / "PET.gro")

    env = os.environ.copy()
    env["GMXLIB"] = os.pathsep.join([str(outdir), str(ff_dir.parent)])

    # 1) Protein topology
    highlight("STEP 1/6", "pdb2gmx: generate protein topology and coordinates", "1;32")
    run(
        [
            gmx,
            "pdb2gmx",
            "-f",
            str(protein),
            "-o",
            "protein.gro",
            "-p",
            "topol.top",
            "-i",
            "posre.itp",
            "-ff",
            ff_dir.stem,
            "-water",
            args.water,
            "-ignh",
        ],
        cwd=outdir,
        env=env,
    )

    molname = parse_moleculetype(itp_dst)
    patch_topol(outdir / "topol.top", itp_dst.name, molname, ff_add_dst.name)

    # Auto-box suggestion/apply based on protein + PET geometry.
    box_hint = suggest_box_distance(outdir / "protein.gro", outdir / "PET.gro", args.box_distance, args.pbc_buffer)
    highlight(
        "INFO",
        (
            "protein span=({:.3f},{:.3f},{:.3f}) nm; "
            "PET span=({:.3f},{:.3f},{:.3f}) nm; PET end-to-end={:.3f} nm"
        ).format(
            box_hint["protein_x"],
            box_hint["protein_y"],
            box_hint["protein_z"],
            box_hint["pet_x"],
            box_hint["pet_y"],
            box_hint["pet_z"],
            box_hint["pet_e2e"],
        ),
        "1;34",
    )
    highlight(
        "INFO",
        (
            "current box_distance={:.3f} nm, required_d>={:.3f} nm "
            "(target min box {:.3f} nm)"
        ).format(args.box_distance, box_hint["required_d"], box_hint["target_min_box"]),
        "1;34",
    )
    if args.auto_box == "apply" and box_hint["suggested_d"] > args.box_distance:
        args.box_distance = box_hint["suggested_d"]
        highlight("INFO", f"auto_box=apply -> using box_distance={args.box_distance:.3f} nm", "1;33")
    elif args.auto_box == "suggest" and box_hint["suggested_d"] > args.box_distance:
        highlight(
            "WARN",
            f"Suggested box_distance >= {box_hint['suggested_d']:.3f} nm for this protein/PET pair",
            "1;33",
        )

    # 2) Define box
    highlight("STEP 2/6", "editconf: define box", "1;32")
    run(
        [
            gmx,
            "editconf",
            "-f",
            "protein.gro",
            "-o",
            "box.gro",
            "-c",
            "-d",
            str(args.box_distance),
            "-bt",
            args.box_type,
        ],
        cwd=outdir,
        env=env,
    )

    # 3) Insert PET into box
    highlight("STEP 3/6", "insert-molecules: insert PET into box", "1;32")
    insert_cmd = [
        gmx,
        "insert-molecules",
        "-f",
        "box.gro",
        "-ci",
        "PET.gro",
        "-nmol",
        "1",
        "-o",
        "conf_pet.gro",
        "-try",
        "1000",
    ]
    if args.insert_seed is not None:
        insert_cmd += ["-seed", str(args.insert_seed)]
    run(
        insert_cmd,
        cwd=outdir,
        env=env,
    )

    # 4) Solvate
    highlight("STEP 4/6", "solvate: add water", "1;32")
    run(
        [
            gmx,
            "solvate",
            "-cp",
            "conf_pet.gro",
            "-cs",
            args.water_coord,
            "-p",
            "topol.top",
            "-o",
            "box_sol.gro",
        ],
        cwd=outdir,
        env=env,
    )

    # 5) Ions
    highlight("STEP 5/6", "grompp + genion: add ions", "1;32")
    run(
        [
            gmx,
            "grompp",
            "-f",
            str(ions_mdp),
            "-c",
            "box_sol.gro",
            "-p",
            "topol.top",
            "-o",
            "ions.tpr",
            "-maxwarn",
            "20",
        ],
        cwd=outdir,
        env=env,
    )
    genion_cmd = [
        gmx,
        "genion",
        "-s",
        "ions.tpr",
        "-p",
        "topol.top",
        "-neutral",
        "-conc",
        str(args.salt),
        "-pname",
        "SOD",
        "-nname",
        "CLA",
        "-o",
        "box_sol_ion.gro",
    ]
    if args.genion_seed is not None:
        genion_cmd += ["-seed", str(args.genion_seed)]
    run(
        genion_cmd,
        cwd=outdir,
        env=env,
        stdin=f"{args.genion_group}\n",
    )

    # 6) Index file with PETL + union groups
    highlight("STEP 6/6", "make_ndx: build index groups", "1;32")
    run(
        [gmx, "make_ndx", "-f", "box_sol_ion.gro", "-o", "index.ndx"],
        cwd=outdir,
        env=env,
        stdin="q\n",
    )
    ensure_groups(outdir / "index.ndx", outdir / "box_sol_ion.gro")

    if args.post_check:
        check_script = script_dir.parent / "check_system_build.py"
        if check_script.exists():
            highlight("STEP 7/7", "check_system_build: sanity report", "1;32")
            check_cmd = [
                sys.executable,
                str(check_script),
                "--gro",
                "box_sol_ion.gro",
                "--itp",
                itp_dst.name,
                "--ligand-resname",
                "PETL",
            ]
            if args.check_strict:
                check_cmd.append("--strict")
            run(check_cmd, cwd=outdir, env=env)
        else:
            highlight("WARN", f"check_system_build.py not found: {check_script}", "1;33")

    highlight("DONE", f"system_build outputs in: {outdir}", "1;32")


if __name__ == "__main__":
    main()
