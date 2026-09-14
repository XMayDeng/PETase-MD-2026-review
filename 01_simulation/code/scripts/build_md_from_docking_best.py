#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import shutil
import subprocess
from pathlib import Path


def _load_build_system_module():
    path = Path(__file__).resolve().with_name("build_system.py")
    spec = importlib.util.spec_from_file_location("build_system", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


BS = _load_build_system_module()


def _load_box_protocol_module():
    path = Path(__file__).resolve().parents[1] / "pipelines" / "box_protocol.py"
    spec = importlib.util.spec_from_file_location("box_protocol", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


BP = _load_box_protocol_module()


def resolve_launch_path(value: str | Path, launch_cwd: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (launch_cwd / path).resolve()


def resolve_repo_path(value: str | Path, repo_root: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return repo_root / path


def run(cmd: list[str], cwd: Path | None = None, stdin: str | None = None, env: dict | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, input=stdin, text=True, check=True, env=env)


def split_best_model(src: Path, protein_out: Path, pet_out: Path) -> None:
    protein_lines = []
    pet_lines = []
    for line in Path(src).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21:22].strip()
        if chain == "A":
            protein_lines.append(line + "\n")
        elif chain == "B":
            pet_lines.append(line + "\n")
    if not protein_lines:
        raise ValueError(f"No chain A atoms found in {src}")
    if not pet_lines:
        raise ValueError(f"No chain B atoms found in {src}")
    protein_out.parent.mkdir(parents=True, exist_ok=True)
    pet_out.parent.mkdir(parents=True, exist_ok=True)
    protein_out.write_text("".join(protein_lines) + "TER\nEND\n", encoding="utf-8")
    pet_out.write_text("".join(pet_lines) + "TER\nEND\n", encoding="utf-8")


def convert_haddock_pet_to_petl(src: Path, dst: Path) -> None:
    out = []
    for raw in Path(src).read_text().splitlines():
        line = raw.rstrip("\n")
        if line.startswith(("ATOM", "HETATM")):
            if len(line) < 80:
                line = line.ljust(80)
            line = f"{line[:17]}PETL{line[21:]}"
            out.append(line)
        elif line.startswith(("TER", "END")):
            out.append(line)
    if not out:
        raise ValueError(f"No PET lines found in {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(out) + "\n", encoding="utf-8")


def merge_gro_files(protein_gro: Path, pet_gro: Path, merged_gro: Path, title: str = "complex_docked") -> None:
    prot_lines = Path(protein_gro).read_text().splitlines()
    pet_lines = Path(pet_gro).read_text().splitlines()
    prot_n = int(prot_lines[1].strip())
    pet_n = int(pet_lines[1].strip())
    atom_lines = prot_lines[2 : 2 + prot_n] + pet_lines[2 : 2 + pet_n]
    box_line = prot_lines[2 + prot_n]
    merged = [title, f"{prot_n + pet_n:5d}", *atom_lines, box_line]
    merged_gro.parent.mkdir(parents=True, exist_ok=True)
    merged_gro.write_text("\n".join(merged) + "\n", encoding="utf-8")


def build_stage_paths(build_root: Path) -> dict[str, Path]:
    return {
        "root": build_root,
        "protein": build_root / "protein",
        "PET": build_root / "PET",
        "assemble": build_root / "assemble",
        "solvate": build_root / "solvate",
        "ionize": build_root / "ionize",
        "checks": build_root / "checks",
        "outputs": build_root / "outputs",
    }


def pet_source_kind_from_kind(pet_kind: str) -> str:
    upper = pet_kind.upper()
    if upper.startswith("PET_L"):
        return upper
    match = re.fullmatch(r"PET(\d+)", upper)
    if match:
        return f"PET_L{match.group(1)}"
    return upper


def box_dim_nm_for(pet_kind: str, override: float | None = None) -> float:
    if override is not None:
        return float(override)
    return float(BP.box_dim_for(pet_source_kind_from_kind(pet_kind)))


def editconf_box_command(gmx: str, complex_gro: Path, box_gro: Path,
                         box_dim_nm: float) -> list[str]:
    dim = str(box_dim_nm)
    return [
        gmx,
        "editconf",
        "-f", str(complex_gro),
        "-o", str(box_gro),
        "-c",
        "-box", dim, dim, dim,
        "-bt", "cubic",
    ]


def write_build_provenance(path: Path, *, args: argparse.Namespace,
                           best_model: Path, pet_source_kind: str,
                           box_dim_nm: float, gmx: str) -> None:
    data = {
        "schema_version": "v1.1",
        "best_model": str(best_model),
        "protein_id": args.protein_id or None,
        "pet_kind": args.pet_kind,
        "pet_source_kind": pet_source_kind,
        "box_type": "cubic",
        "box_dim_nm": float(box_dim_nm),
        "salt_concentration_M": float(args.salt),
        "genion_seed": int(args.genion_seed),
        "gmx_bin": str(gmx),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--best-model", required=True)
    p.add_argument("--build-root")
    p.add_argument("--build-dir")
    p.add_argument("--gmx-bin", default="gmx")
    p.add_argument("--pet-kind", default="PET_L20")
    p.add_argument("--box-dim-nm", type=float, default=None,
                   help="Optional cubic edge override in nm. Defaults to the "
                        "v1.1 chain-length-aware box protocol.")
    p.add_argument("--protein-id", default="")
    p.add_argument("--ff-dir", default="./inputs/forcefields/charmm36-jul2021.ff")
    p.add_argument("--pet-itp")
    p.add_argument("--ff-itp")
    p.add_argument("--ions-mdp", default="./configurations/templates/md/ions.mdp")
    p.add_argument("--em-mdp", default="./configurations/templates/md/EM.mdp")
    p.add_argument("--nvt-mdp", default="./configurations/templates/md/NVT.mdp")
    p.add_argument("--salt", type=float, default=0.1)
    p.add_argument("--genion-seed", type=int, default=20260308)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    launch_cwd = Path.cwd()
    build_root_value = args.build_root or args.build_dir
    if not build_root_value:
        raise SystemExit("one of --build-root or --build-dir is required")
    build_root = resolve_launch_path(build_root_value, launch_cwd)
    layout = build_stage_paths(build_root)
    for stage_dir in layout.values():
        stage_dir.mkdir(parents=True, exist_ok=True)

    best_model = resolve_launch_path(args.best_model, launch_cwd)
    pet_source_kind = pet_source_kind_from_kind(args.pet_kind)
    box_dim_nm = box_dim_nm_for(args.pet_kind, args.box_dim_nm)
    ff_dir = resolve_repo_path(args.ff_dir, repo_root)
    pet_toppar_rel = (
        f"./inputs/pets/{pet_source_kind}"
        "/raw/charmm-gui-001/gromacs/toppar"
    )
    pet_itp_rel = args.pet_itp or f"{pet_toppar_rel}/P1.itp"
    ff_itp_rel = args.ff_itp or f"{pet_toppar_rel}/forcefield.itp"
    pet_itp = resolve_repo_path(pet_itp_rel, repo_root)
    ff_itp = resolve_repo_path(ff_itp_rel, repo_root)
    ions_mdp = resolve_repo_path(args.ions_mdp, repo_root)
    em_mdp = resolve_repo_path(args.em_mdp, repo_root)
    nvt_mdp = resolve_repo_path(args.nvt_mdp, repo_root)

    protein_pdb = layout["protein"] / "protein_docked_chainA.pdb"
    pet_haddock_pdb = layout["PET"] / "pet_docked_chainB_haddock.pdb"
    pet_raw_pdb = layout["PET"] / "pet_docked_chainB_raw.pdb"
    protein_gro = layout["protein"] / "protein.gro"
    pet_gro = layout["PET"] / "pet_docked_raw.gro"
    complex_gro = layout["assemble"] / "complex_docked.gro"
    box_gro = layout["solvate"] / "box.gro"
    box_sol = layout["solvate"] / "box_sol.gro"
    box_sol_ion = layout["ionize"] / "box_sol_ion.gro"
    ions_tpr = layout["ionize"] / "ions.tpr"
    index_ndx = layout["outputs"] / "index.ndx"
    topol = layout["outputs"] / "topol.top"
    posre = layout["outputs"] / "posre.itp"
    pet_itp_out = layout["outputs"] / pet_itp.name
    pet_ff_add = layout["outputs"] / "pet_ff_add.itp"
    check_em_tpr = layout["checks"] / "check_EM.tpr"
    check_nvt_tpr = layout["checks"] / "check_NVT.tpr"

    split_best_model(best_model, protein_pdb, pet_haddock_pdb)
    convert_haddock_pet_to_petl(pet_haddock_pdb, pet_raw_pdb)

    if not pet_itp.exists():
        raise FileNotFoundError(pet_itp)
    if not ff_itp.exists():
        raise FileNotFoundError(ff_itp)

    gmx = shutil.which(args.gmx_bin) or BS.find_gmx()
    env = os.environ.copy()
    env["GMXLIB"] = os.pathsep.join([str(layout["outputs"]), str(ff_dir.parent)])

    run(
        [
            gmx,
            "pdb2gmx",
            "-f",
            str(protein_pdb),
            "-o",
            str(protein_gro),
            "-p",
            str(topol),
            "-i",
            str(posre),
            "-ff",
            ff_dir.stem,
            "-water",
            "tip3p",
            "-ignh",
        ],
        cwd=layout["protein"],
        env=env,
    )
    run([gmx, "editconf", "-f", str(pet_raw_pdb), "-o", str(pet_gro)], cwd=layout["PET"], env=env)
    merge_gro_files(protein_gro, pet_gro, complex_gro)

    shutil.copy2(pet_itp, pet_itp_out)
    BS.write_pet_ff_add(ff_itp, pet_ff_add, ff_dir / "ffbonded.itp", ff_dir / "ffnonbonded.itp")
    molname = BS.parse_moleculetype(pet_itp_out)
    BS.patch_topol(topol, pet_itp.name, molname, "pet_ff_add.itp")

    run(editconf_box_command(gmx, complex_gro, box_gro, box_dim_nm), cwd=layout["assemble"], env=env)
    run([gmx, "solvate", "-cp", str(box_gro), "-cs", "spc216.gro", "-p", str(topol), "-o", str(box_sol)], cwd=layout["solvate"], env=env)
    run([gmx, "grompp", "-f", str(ions_mdp), "-c", str(box_sol), "-p", str(topol), "-o", str(ions_tpr), "-maxwarn", "20"], cwd=layout["ionize"], env=env)
    run([gmx, "genion", "-s", str(ions_tpr), "-o", str(box_sol_ion), "-p", str(topol), "-neutral", "-conc", str(args.salt), "-pname", "SOD", "-nname", "CLA", "-seed", str(args.genion_seed)], cwd=layout["ionize"], stdin="SOL\n", env=env)
    shutil.copy2(box_sol_ion, layout["outputs"] / "box_sol_ion.gro")
    run([gmx, "make_ndx", "-f", str(box_sol_ion), "-o", str(index_ndx)], cwd=layout["ionize"], stdin="q\n", env=env)
    BS.ensure_groups(index_ndx, box_sol_ion)
    run([gmx, "grompp", "-f", str(em_mdp), "-c", str(box_sol_ion), "-r", str(box_sol_ion), "-p", str(topol), "-n", str(index_ndx), "-o", str(check_em_tpr), "-maxwarn", "20"], cwd=layout["checks"], env=env)
    run([gmx, "grompp", "-f", str(nvt_mdp), "-c", str(box_sol_ion), "-r", str(box_sol_ion), "-p", str(topol), "-n", str(index_ndx), "-o", str(check_nvt_tpr), "-maxwarn", "20"], cwd=layout["checks"], env=env)
    write_build_provenance(
        layout["outputs"] / "build_provenance.json",
        args=args,
        best_model=best_model,
        pet_source_kind=pet_source_kind,
        box_dim_nm=box_dim_nm,
        gmx=gmx,
    )


if __name__ == "__main__":
    main()
