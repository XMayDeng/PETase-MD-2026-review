#!/usr/bin/env python3
"""R1 Phase B multi-conformer Kabsch screen orchestrator (spec §4.1).

Wraps the single-conformer ``screen_vina_pet_longchain_candidates.py`` to
sweep K=5 pre-equilibrated PET conformers from the R1 Phase A asset:

    for conformer in K=5 conformers of the long-chain PET asset:
        for window in candidate windows:
            for fwd_rev in {forward, reverse}:
                Kabsch align (conformer, window, direction) to docked L4
                score = motif_rmsd + clash + affinity + post_kabsch_envelope

Results are aggregated into ``screening_summary_multi_conformer.json`` and
tagged with:

* ``conformer_idx``                 — 1…K from R1 Phase A asset
* ``chain_extension_direction``     — head-side / tail-side / bidirectional
* ``post_kabsch_envelope_nm3``      — bounding-box volume of transformed PET
* ``ranking_score``                 — composite score (lower is better)

Downstream consumer: ``select_v1_binding_candidate.py --per-direction``.

Usage::

    python screen_vina_pet_multi_conformer.py \\
        --conformer-dir inputs/pets/PET_L10/preequilibrated/ \\
        --target-pet-kind PET_L10 \\
        --target-unit-map inputs/pets/PET_L10/alignment/pet_l10_unit_map.yaml \\
        --docking-dir <runtime>/02_docking/vina_pet_l4/outputs \\
        --protein-pdb <runtime>/02_docking/inputs/receptor.pdb \\
        --active-site-registry inputs/active_site_registry.yaml \\
        --ser-resid 131 --ser-atom-name OG --his-resid 209 --acid-resid 177 \\
        --out-dir <runtime>/03_pet_alignment_screen \\
        [--ranking-weights motif_rmsd=0.4,clash=0.3,affinity=0.2,envelope=0.1]
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import datetime
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = REPO_ROOT / "."
SCRIPT_DIR = Path(__file__).resolve().parent

DEFAULT_RANKING_WEIGHTS = {
    "motif_rmsd": 0.4,
    "clash": 0.3,
    "affinity": 0.2,
    "envelope": 0.1,
}

DEFAULT_PET_RESNAME_PREFIXES = ("PET", "PEA", "PEM", "PEB", "PETL")


def default_workers() -> int:
    """Return the default total worker budget for this screen.

    In Slurm, use the allocated CPU count automatically so existing
    ``run_case_v1_1.py`` calls become parallel without changing the command
    line. Outside Slurm, keep the old conservative serial default.
    """
    for key in ("AUTOPROPET_SCREEN_WORKERS", "SLURM_CPUS_PER_TASK"):
        value = os.environ.get(key)
        if not value:
            continue
        try:
            workers = int(value)
        except ValueError:
            continue
        if workers > 0:
            return workers
    return 1


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


R1 = _load_module(SIM_ROOT / "code/pipelines/r1_direction.py", "r1_direction")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_ranking_weights(spec: str | None) -> dict[str, float]:
    """Parse ``motif_rmsd=0.4,clash=0.3,affinity=0.2,envelope=0.1`` into a dict."""
    if not spec:
        return DEFAULT_RANKING_WEIGHTS.copy()
    weights = DEFAULT_RANKING_WEIGHTS.copy()
    for pair in spec.split(","):
        if not pair.strip():
            continue
        key, _, value = pair.partition("=")
        key = key.strip()
        if key not in weights:
            raise ValueError(f"Unknown ranking weight key: {key!r}")
        weights[key] = float(value.strip())
    return weights


def extract_pet_heavy_coords(
    pdb_path: Path,
    resname_prefixes: tuple[str, ...] = DEFAULT_PET_RESNAME_PREFIXES,
) -> list[tuple[float, float, float]]:
    """Return PET heavy-atom coordinates (nm) from a PDB file.

    Filters atoms whose residue name begins with any of ``resname_prefixes``
    (matched case-insensitively) and whose atom name does not start with H.
    PDB coordinates are in Å; converted to nm for envelope consistency
    with ``box_protocol``.
    """
    coords: list[tuple[float, float, float]] = []
    if not pdb_path.exists():
        return coords
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if len(line) < 54:
            continue
        atom_name = line[12:16].strip()
        if atom_name.upper().startswith("H"):
            continue
        resname = line[17:21].strip().upper()
        if not any(resname.startswith(p) for p in resname_prefixes):
            continue
        try:
            x = float(line[30:38]) * 0.1   # Å → nm
            y = float(line[38:46]) * 0.1
            z = float(line[46:54]) * 0.1
        except ValueError:
            continue
        coords.append((x, y, z))
    return coords


def run_subprocess(cmd: list[str], *, cwd: Path | None = None,
                   check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True,
                          text=True, check=check)


def convert_gro_to_pdb(gro: Path, pdb: Path, gmx: Path) -> None:
    """Convert a single-frame .gro to .pdb via gmx editconf."""
    run_subprocess([str(gmx), "editconf",
                    "-f", str(gro), "-o", str(pdb)])


def run_single_conformer_screen(
    *,
    script: Path,
    conformer_pdb: Path,
    target_unit_map: Path,
    target_pet_kind: str,
    source_pet_kind: str,
    docking_dir: Path,
    poses_pdbqt: Path | None,
    ligand_atom_map: Path | None,
    protein_pdb: Path,
    active_site_registry: Path,
    ser_resid: int,
    ser_atom_name: str,
    his_resid: int,
    acid_resid: int,
    out_dir: Path,
    extra_args: list[str],
) -> Path:
    """Invoke the single-conformer screen against one conformer PDB.

    Returns the path of the produced screening_summary.json.
    """
    cmd = [
        sys.executable, str(script),
        "--target-pet-pdb", str(conformer_pdb),
        "--target-unit-map", str(target_unit_map),
        "--target-pet-kind", target_pet_kind,
        "--source-pet-kind", source_pet_kind,
        "--protein-pdb", str(protein_pdb),
        "--active-site-registry", str(active_site_registry),
        "--ser-resid", str(ser_resid),
        "--ser-atom-name", str(ser_atom_name),
        "--his-resid", str(his_resid),
        "--acid-resid", str(acid_resid),
        "--docking-dir", str(docking_dir),
        "--out-dir", str(out_dir),
        "--force",
    ]
    if poses_pdbqt is not None:
        cmd.extend(["--poses-pdbqt", str(poses_pdbqt)])
    if ligand_atom_map is not None:
        cmd.extend(["--ligand-atom-map", str(ligand_atom_map)])
    cmd.extend(extra_args)

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr or "")
        raise RuntimeError(
            f"screen_vina_pet_longchain_candidates returned {proc.returncode} "
            f"for conformer {conformer_pdb.name}"
        )
    summary = out_dir / "screening_summary.json"
    if not summary.exists():
        raise FileNotFoundError(
            f"Expected screening_summary.json under {out_dir} but it is missing"
        )
    return summary


def screen_conformer_worker(payload: dict) -> tuple[int, str]:
    """Convert and screen one PET conformer.

    Each worker writes to a conformer-specific output directory. The
    single-conformer screen may itself use pose-level workers, but all writes
    stay below that conformer directory.
    """
    conformer_path = Path(payload["conformer_path"])
    conformer_idx = int(payload["conformer_idx"])
    out_dir = Path(payload["out_dir"])
    work_dir = Path(payload["work_dir"])
    conformer_pdb = work_dir / f"conformer_{conformer_idx:03d}.pdb"
    convert_gro_to_pdb(conformer_path, conformer_pdb, Path(payload["gmx"]))
    extra_args = list(payload.get("extra_args") or [])
    pose_workers = int(payload.get("pose_workers") or 1)
    extra_args.extend(["--workers", str(max(1, pose_workers))])
    summary = run_single_conformer_screen(
        script=Path(payload["script"]),
        conformer_pdb=conformer_pdb,
        target_unit_map=Path(payload["target_unit_map"]),
        target_pet_kind=str(payload["target_pet_kind"]),
        source_pet_kind=str(payload["source_pet_kind"]),
        docking_dir=Path(payload["docking_dir"]),
        poses_pdbqt=Path(payload["poses_pdbqt"]) if payload.get("poses_pdbqt") else None,
        ligand_atom_map=Path(payload["ligand_atom_map"]) if payload.get("ligand_atom_map") else None,
        protein_pdb=Path(payload["protein_pdb"]),
        active_site_registry=Path(payload["active_site_registry"]),
        ser_resid=int(payload["ser_resid"]),
        ser_atom_name=str(payload["ser_atom_name"]),
        his_resid=int(payload["his_resid"]),
        acid_resid=int(payload["acid_resid"]),
        out_dir=out_dir,
        extra_args=extra_args,
    )
    return conformer_idx, str(summary)


def aggregate_candidates(
    *,
    per_conformer_summaries: list[tuple[int, Path]],
    target_pet_kind: str,
    weights: dict[str, float],
) -> list[dict]:
    """Tag and combine candidate records from all conformer sub-runs."""
    aggregated: list[dict] = []
    for conformer_idx, summary_path in per_conformer_summaries:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        candidates = summary.get("top_candidates") or summary.get("ranked_candidates") or []
        for cand in candidates:
            window_start = cand.get("window_start")
            window_end = cand.get("window_end")
            if window_start is None or window_end is None:
                continue
            try:
                ext_dir = R1.infer_chain_extension_direction(
                    (int(window_start), int(window_end)), target_pet_kind
                )
            except ValueError:
                ext_dir = "unknown"

            envelope_nm3 = None
            complex_pdb = cand.get("complex_pdb")
            if complex_pdb:
                coords = extract_pet_heavy_coords(Path(complex_pdb))
                if coords:
                    envelope_nm3 = R1.compute_post_kabsch_envelope_nm3(coords)

            motif_rmsd = cand.get("motif_rmsd_angstrom")
            clash = cand.get("contacts_lt_1a")
            affinity = cand.get("affinity_kcal_mol")
            score = None
            if (motif_rmsd is not None and clash is not None
                    and affinity is not None and envelope_nm3 is not None):
                try:
                    score = R1.ranking_score(
                        motif_rmsd_angstrom=float(motif_rmsd),
                        contacts_lt_1a=int(clash),
                        affinity_kcal_mol=float(affinity),
                        envelope_nm3=float(envelope_nm3),
                        weights=weights,
                    )
                except (TypeError, ValueError):
                    score = None

            enriched = dict(cand)
            enriched["conformer_idx"] = conformer_idx
            enriched["chain_extension_direction"] = ext_dir
            enriched["post_kabsch_envelope_nm3"] = envelope_nm3
            enriched["ranking_score"] = score
            aggregated.append(enriched)
    aggregated.sort(
        key=lambda r: (r["ranking_score"] is None, r.get("ranking_score") or float("inf"))
    )
    return aggregated


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--conformer-dir", required=True, type=Path,
                   help="R1 Phase A asset dir containing conformer_001..N.gro.")
    p.add_argument("--target-pet-kind", required=True,
                   choices=sorted(R1.PET_RESIDUE_COUNT),
                   help="Long-chain PET target (PET_L10 or PET_L20).")
    p.add_argument("--target-unit-map", required=True, type=Path,
                   help="Per-target atom map yaml shared across all K conformers.")
    p.add_argument("--source-pet-kind", default="PET_L4",
                   help="Source PET kind docked by Vina (default PET_L4).")
    p.add_argument("--docking-dir", required=True, type=Path,
                   help="Directory containing Vina docking outputs.")
    p.add_argument("--poses-pdbqt", type=Path, default=None)
    p.add_argument("--ligand-atom-map", type=Path, default=None)
    p.add_argument("--protein-pdb", required=True, type=Path)
    p.add_argument("--active-site-registry", required=True, type=Path)
    p.add_argument("--ser-resid", required=True, type=int)
    p.add_argument("--ser-atom-name", default="OG",
                   help="Catalytic serine atom name forwarded to the single-conformer screen.")
    p.add_argument("--his-resid", required=True, type=int)
    p.add_argument("--acid-resid", required=True, type=int)
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Multi-conformer screen output dir.")
    p.add_argument("--ranking-weights", type=str, default=None,
                   help="Comma-separated weight overrides "
                        "(e.g. motif_rmsd=0.4,clash=0.3,affinity=0.2,envelope=0.1).")
    p.add_argument("--gmx", type=Path,
                   default=Path(os.environ.get("GMX_BIN", "/usr/local/gromacs/bin/gmx")))
    p.add_argument("--screen-script", type=Path,
                   default=SCRIPT_DIR / "screen_vina_pet_longchain_candidates.py")
    p.add_argument("--extra-screen-arg", action="append", default=None,
                   help="Verbatim extra arg forwarded to the single-conformer screen "
                        "(may be repeated).")
    p.add_argument("--workers", type=int, default=default_workers(),
                   help="Total worker budget. Defaults to AUTOPROPET_SCREEN_WORKERS, "
                        "then SLURM_CPUS_PER_TASK, then 1.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    conformer_paths = sorted(args.conformer_dir.glob("conformer_*.gro"))
    if not conformer_paths:
        print(f"ERROR: no conformer_*.gro files found under {args.conformer_dir}",
              file=sys.stderr)
        return 2
    if not args.target_unit_map.exists():
        print(f"ERROR: target unit map missing: {args.target_unit_map}",
              file=sys.stderr)
        return 2
    if not args.screen_script.exists():
        print(f"ERROR: single-conformer screen script missing: {args.screen_script}",
              file=sys.stderr)
        return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = args.out_dir / "_work"
    work_dir.mkdir(exist_ok=True)

    weights = parse_ranking_weights(args.ranking_weights)
    extra_args = args.extra_screen_arg or []

    total_workers = max(1, int(args.workers))
    conformer_workers = min(len(conformer_paths), total_workers)
    pose_workers = max(1, total_workers // max(1, conformer_workers))
    per_conformer_summaries: list[tuple[int, Path]] = []
    conformer_jobs: list[dict] = []
    for conformer_path in conformer_paths:
        m = re.search(r"conformer_(\d+)", conformer_path.stem)
        if not m:
            continue
        conformer_idx = int(m.group(1))
        sub_out = args.out_dir / f"conformer_{conformer_idx:03d}"
        conformer_jobs.append(
            {
                "conformer_path": str(conformer_path),
                "conformer_idx": conformer_idx,
                "out_dir": str(sub_out),
                "work_dir": str(work_dir),
                "gmx": str(args.gmx),
                "script": str(args.screen_script),
                "target_unit_map": str(args.target_unit_map),
                "target_pet_kind": args.target_pet_kind,
                "source_pet_kind": args.source_pet_kind,
                "docking_dir": str(args.docking_dir),
                "poses_pdbqt": str(args.poses_pdbqt) if args.poses_pdbqt else None,
                "ligand_atom_map": str(args.ligand_atom_map) if args.ligand_atom_map else None,
                "protein_pdb": str(args.protein_pdb),
                "active_site_registry": str(args.active_site_registry),
                "ser_resid": args.ser_resid,
                "ser_atom_name": args.ser_atom_name,
                "his_resid": args.his_resid,
                "acid_resid": args.acid_resid,
                "extra_args": list(extra_args),
                "pose_workers": pose_workers,
            }
        )

    print(
        f"[parallel] total_workers={total_workers} "
        f"conformer_workers={conformer_workers} pose_workers={pose_workers}",
        file=sys.stderr,
    )
    if conformer_workers > 1:
        with ProcessPoolExecutor(max_workers=conformer_workers) as executor:
            futures = [
                executor.submit(screen_conformer_worker, payload)
                for payload in conformer_jobs
            ]
            for future in as_completed(futures):
                conformer_idx, summary = future.result()
                per_conformer_summaries.append((conformer_idx, Path(summary)))
    else:
        for payload in conformer_jobs:
            conformer_idx, summary = screen_conformer_worker(payload)
            per_conformer_summaries.append((conformer_idx, Path(summary)))
    per_conformer_summaries.sort(key=lambda item: item[0])

    candidates = aggregate_candidates(
        per_conformer_summaries=per_conformer_summaries,
        target_pet_kind=args.target_pet_kind,
        weights=weights,
    )

    out = {
        "schema_version": "v1.1",
        "stage": "R1_phase_B_multi_conformer_screen",
        "target_pet_kind": args.target_pet_kind,
        "source_pet_kind": args.source_pet_kind,
        "n_conformers": len(per_conformer_summaries),
        "parallel": {
            "total_workers": total_workers,
            "conformer_workers": conformer_workers,
            "pose_workers_per_conformer": pose_workers,
        },
        "ranking_weights": weights,
        "top_candidates": candidates,
        "computed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    out_path = args.out_dir / "screening_summary_multi_conformer.json"
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")

    print(f"[done] aggregated {len(candidates)} candidates from "
          f"{len(per_conformer_summaries)} conformers → {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
