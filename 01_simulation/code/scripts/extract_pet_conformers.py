#!/usr/bin/env python3
"""Pool N PET-only trajectories, cluster, and emit the K=5 conformer asset.

PET-only preequilibration provides conformers shared by the docking screen.

Inputs are typically the 3 independent 50 ns trajectories produced by
``run_pet_solo_preeq.py`` for one chain length (L10 or L20). Trajectories
are PBC-unwrapped, concatenated, then clustered with gmx cluster
(linkage method, RMSD cutoff 0.2 nm). The K=5 most populated cluster
centroids become the frozen conformer asset feeding R1 Phase B (Kabsch
alignment screen).

Output asset layout::

    inputs/pets/<kind>/preequilibrated/
    ├── conformer_001.gro ... conformer_005.gro
    ├── cluster_summary.yaml          # Rg / E2E / population per conformer
    ├── cluster_table.txt             # raw gmx cluster population table
    └── MANIFEST.json                 # SHA-256 of every output file
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = REPO_ROOT / "."

DEFAULT_K = 5
DEFAULT_RMSD_CUTOFF_NM = 0.2
SUPPORTED_KINDS = {"PET_L10", "PET_L20"}


# ---------------------------------------------------------------------------
# Pure helpers (unit-test target)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GroAtom:
    resid: int
    resname: str
    atom_name: str
    atom_serial: int
    x_nm: float
    y_nm: float
    z_nm: float


# Fixed-width GROMACS .gro atom line layout (positions per gro file spec).
#   cols  1- 5 : resid
#   cols  6-10 : resname
#   cols 11-15 : atom_name
#   cols 16-20 : atom_serial
#   cols 21-28 : x
#   cols 29-36 : y
#   cols 37-44 : z

def _parse_gro_atom_line(line: str) -> GroAtom | None:
    if len(line) < 44:
        return None
    try:
        return GroAtom(
            resid=int(line[0:5]),
            resname=line[5:10].strip(),
            atom_name=line[10:15].strip(),
            atom_serial=int(line[15:20]),
            x_nm=float(line[20:28]),
            y_nm=float(line[28:36]),
            z_nm=float(line[36:44]),
        )
    except ValueError:
        return None


def parse_gro(path: Path) -> list[GroAtom]:
    """Parse a single-frame .gro file. Returns atom list (excludes box line)."""
    text = path.read_text(encoding="utf-8").splitlines()
    if len(text) < 3:
        raise ValueError(f"Truncated .gro: {path}")
    n_atoms = int(text[1].strip())
    atom_lines = text[2:2 + n_atoms]
    atoms: list[GroAtom] = []
    for line in atom_lines:
        parsed = _parse_gro_atom_line(line)
        if parsed is not None:
            atoms.append(parsed)
    if len(atoms) != n_atoms:
        raise ValueError(f"Parsed {len(atoms)} atoms but header says {n_atoms} ({path})")
    return atoms


def compute_end_to_end_distance(atoms: list[GroAtom]) -> float:
    """End-to-end distance (nm) of the polymer.

    Defined as the Euclidean distance between the heaviest non-hydrogen atom
    of the first PETL residue and that of the last PETL residue. Falls back
    to the first/last atoms if no PETL residue is found.
    """
    pet_atoms = [a for a in atoms if a.resname.upper().startswith("PET")]
    if not pet_atoms:
        pet_atoms = atoms
    first_resid = min(a.resid for a in pet_atoms)
    last_resid = max(a.resid for a in pet_atoms)

    def _pick(resid: int) -> GroAtom:
        candidates = [a for a in pet_atoms if a.resid == resid
                      and not a.atom_name.upper().startswith("H")]
        if not candidates:
            candidates = [a for a in pet_atoms if a.resid == resid]
        # Prefer terminal heavy-oxygen if present, else first atom.
        for hint in ("OZ", "O1", "O11", "OE1"):
            for a in candidates:
                if a.atom_name == hint:
                    return a
        return candidates[0]

    a = _pick(first_resid)
    b = _pick(last_resid)
    dx, dy, dz = (a.x_nm - b.x_nm, a.y_nm - b.y_nm, a.z_nm - b.z_nm)
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_cluster_log(log_path: Path) -> dict[int, int]:
    """Parse `gmx cluster` log lines like ``cl. |  #st  | ...`` to extract
    per-cluster member counts. Cluster IDs are 1-indexed.

    Returns ``{cluster_id: member_count}``.
    """
    populations: dict[int, int] = {}
    text = log_path.read_text(encoding="utf-8", errors="replace")
    pattern = re.compile(
        r"^\s*(\d+)\s*\|\s*(\d+)(?:\s+[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)?\s*\|"
    )
    in_table = False
    for line in text.splitlines():
        if "cl." in line and "#st" in line:
            in_table = True
            continue
        if in_table:
            m = pattern.match(line)
            if m:
                cid = int(m.group(1))
                count = int(m.group(2))
                populations[cid] = count
            elif line.strip() and not pattern.match(line):
                if line.strip().startswith("---") or "rmsd" in line.lower():
                    continue
                # End of table on first non-matching non-blank line.
                if not line.lstrip().startswith(("|", " ")):
                    break
    return populations


def select_top_cluster_ids(populations: dict[int, int], k: int) -> list[int]:
    """Return cluster IDs sorted by decreasing population, then cluster ID."""
    ranked = sorted(populations.items(), key=lambda item: (-item[1], item[0]))
    return [cid for cid, _count in ranked[:k]]


def build_cluster_summary(*, pet_kind: str, k: int, rmsd_cutoff_nm: float,
                          conformer_records: list[dict],
                          source_trajectories: list[str],
                          pool_frames: int | None,
                          gromacs_version: str,
                          begin_ps: float | None = None,
                          end_ps: float | None = None,
                          ) -> dict:
    return {
        "schema_version": "v1.1",
        "pet_kind": pet_kind,
        "asset_purpose": "R1 Phase A conformer asset (consumed by Phase B Kabsch)",
        "cluster_method": "linkage",
        "rmsd_cutoff_nm": rmsd_cutoff_nm,
        "k": k,
        "pool_strategy": "cross_trajectory",
        "frame_filter": {
            "begin_ps": begin_ps,
            "end_ps": end_ps,
        },
        "source_trajectories": list(source_trajectories),
        "total_frames_pooled": pool_frames,
        "gromacs_version": gromacs_version,
        "conformers": conformer_records,
        "computed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }


# ---------------------------------------------------------------------------
# Subprocess + cluster pipeline
# ---------------------------------------------------------------------------

class _Runner:
    def __init__(self, dry_run: bool):
        self.dry_run = dry_run

    def run(self, cmd: list[str], *, cwd: Path | None = None,
            stdin: str | None = None) -> None:
        printable = " ".join(str(c) for c in cmd)
        prefix = "[dry-run] " if self.dry_run else "[run] "
        cwd_label = f"  (cwd={cwd})" if cwd else ""
        print(f"{prefix}{printable}{cwd_label}", flush=True)
        if self.dry_run:
            return
        subprocess.run(cmd, cwd=cwd, input=stdin, text=True, check=True)


def gmx_version(gmx_bin: Path) -> str:
    try:
        out = subprocess.run([str(gmx_bin), "--version"],
                             capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            if line.strip().startswith("GROMACS"):
                return line.strip().rstrip(":")
        return out.splitlines()[0].strip() if out else "unknown"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def trjconv_pbc_whole(xtc: Path, tpr: Path, out_xtc: Path, gmx: Path,
                      runner: _Runner, *, begin_ps: float | None = None,
                      end_ps: float | None = None) -> None:
    # Restrict to PETL atoms only: cluster runs on PET conformation alone,
    # not water + ions which would swamp the RMSD signal and balloon the
    # pairwise distance matrix.
    cmd = [str(gmx), "trjconv",
           "-f", str(xtc), "-s", str(tpr),
           "-o", str(out_xtc), "-pbc", "whole"]
    if begin_ps is not None:
        cmd.extend(["-b", f"{begin_ps:g}"])
    if end_ps is not None:
        cmd.extend(["-e", f"{end_ps:g}"])
    runner.run(cmd, stdin="PETL\n")


def trjcat(xtcs: list[Path], out_xtc: Path, gmx: Path, runner: _Runner) -> None:
    runner.run([str(gmx), "trjcat", "-f", *[str(x) for x in xtcs],
                "-o", str(out_xtc), "-cat"])


def cluster_pool(pool_xtc: Path, tpr: Path, work_dir: Path, rmsd_cutoff_nm: float,
                 gmx: Path, runner: _Runner) -> tuple[Path, Path]:
    """Run ``gmx cluster -method linkage``. Returns (clusters.pdb, log_path)."""
    clusters_pdb = work_dir / "clusters.pdb"
    log_path = work_dir / "cluster.log"
    cmd = [str(gmx), "cluster",
           "-f", str(pool_xtc), "-s", str(tpr),
           "-method", "linkage", "-cutoff", f"{rmsd_cutoff_nm}",
           "-cl", str(clusters_pdb),
           "-g", str(log_path),
           "-wcl", str(DEFAULT_K)]
    runner.run(cmd, cwd=work_dir,
               stdin="PETL\nPETL\n")  # fit + RMSD on PET heavy-atom group
    return clusters_pdb, log_path


def split_clusters_pdb_to_gros(clusters_pdb: Path, out_dir: Path,
                               cluster_ids: list[int], gmx: Path,
                               runner: _Runner) -> list[Path]:
    """``gmx cluster -cl`` emits a multi-model PDB. Split selected cluster
    models into individual conformer_001..K .gro files via gmx editconf.
    """
    text = clusters_pdb.read_text(encoding="utf-8")
    model_re = re.compile(r"^MODEL\s+(\d+)\s*$", re.MULTILINE)
    matches = list(model_re.finditer(text))
    end_positions = [m.start() for m in matches] + [len(text)]
    model_chunks: dict[int, str] = {}
    for i, match in enumerate(matches):
        model_no = int(match.group(1))
        start = match.start()
        stop = end_positions[i + 1]
        chunk = text[start:stop]
        chunk = re.sub(r"^MODEL.*$", "", chunk, count=1, flags=re.MULTILINE)
        chunk = re.sub(r"^ENDMDL\s*$", "END\n", chunk, flags=re.MULTILINE)
        model_chunks[model_no] = chunk
    out_paths: list[Path] = []
    for i, cluster_id in enumerate(cluster_ids, start=1):
        chunk = model_chunks.get(cluster_id)
        if chunk is None:
            continue
        tmp_pdb = out_dir / f"_cluster_{cluster_id:03d}_conformer_{i:03d}.pdb"
        tmp_pdb.write_text(chunk.lstrip("\n"), encoding="utf-8")
        gro_path = out_dir / f"conformer_{i:03d}.gro"
        runner.run([str(gmx), "editconf",
                    "-f", str(tmp_pdb), "-o", str(gro_path)],
                   cwd=out_dir)
        if not runner.dry_run:
            tmp_pdb.unlink(missing_ok=True)
        out_paths.append(gro_path)
    return out_paths


def require_centroid_count(centroids: list[Path], k: int, source: Path) -> None:
    if len(centroids) != k:
        raise ValueError(
            f"expected {k} cluster centroids from {source}, got {len(centroids)}"
        )


def gyrate_rg_nm(gro: Path, gmx: Path, runner: _Runner, work_dir: Path) -> float:
    """Return the Rg (nm) of a single-frame .gro via gmx gyrate."""
    xvg = work_dir / f"_rg_{gro.stem}.xvg"
    runner.run([str(gmx), "gyrate",
                "-f", str(gro), "-s", str(gro),
                "-o", str(xvg)],
               stdin="System\n")
    if runner.dry_run or not xvg.exists():
        return float("nan")
    for line in xvg.read_text().splitlines():
        if line.startswith(("#", "@")) or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                return float(parts[1])
            except ValueError:
                continue
    return float("nan")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pet-kind", required=True, choices=sorted(SUPPORTED_KINDS))
    p.add_argument("--trajectories", required=True, nargs="+", type=Path,
                   help="One or more PET-only MD.xtc files to pool.")
    p.add_argument("--tpr", required=True, type=Path,
                   help="A representative TPR matching the topology (any run's MD.tpr).")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Destination asset directory "
                        "(typically inputs/pets/<kind>/preequilibrated/).")
    p.add_argument("--k", type=int, default=DEFAULT_K)
    p.add_argument("--rmsd-cutoff", type=float, default=DEFAULT_RMSD_CUTOFF_NM,
                   help=f"RMSD linkage cutoff (nm; default {DEFAULT_RMSD_CUTOFF_NM}).")
    p.add_argument("--begin-ps", type=float,
                   help="Drop frames before this time (ps) before pooling, "
                        "e.g. 25000 to cluster only the final 25 ns of a 50 ns run.")
    p.add_argument("--end-ps", type=float,
                   help="Drop frames after this time (ps) before pooling.")
    p.add_argument("--gmx", type=Path,
                   default=Path(os.environ.get("GMX_BIN", "/usr/local/gromacs/bin/gmx")))
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--keep-pool", action="store_true",
                   help="Retain intermediate pool.xtc + clusters.pdb (for debugging).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    # Resolve every input path to absolute up front: cluster_pool /
    # trjconv invocations change cwd to work_dir, and relative paths
    # would dangle.
    args.trajectories = [p.resolve() for p in args.trajectories]
    args.tpr = args.tpr.resolve()
    args.gmx = args.gmx.resolve() if args.gmx.exists() else args.gmx
    if not args.dry_run:
        for xtc in args.trajectories:
            if not xtc.exists():
                print(f"ERROR: trajectory missing: {xtc}", file=sys.stderr)
                return 2
        if not args.tpr.exists():
            print(f"ERROR: tpr missing: {args.tpr}", file=sys.stderr)
            return 2

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = out_dir / "_work"
    work_dir.mkdir(exist_ok=True)
    runner = _Runner(dry_run=args.dry_run)

    # 1) PBC-unwrap each input.
    unwrapped: list[Path] = []
    for i, xtc in enumerate(args.trajectories, start=1):
        out_xtc = work_dir / f"unwrapped_{i:03d}.xtc"
        trjconv_pbc_whole(
            xtc, args.tpr, out_xtc, args.gmx, runner,
            begin_ps=args.begin_ps, end_ps=args.end_ps,
        )
        unwrapped.append(out_xtc)

    # 2) Pool concatenate.
    pool_xtc = work_dir / "pool.xtc"
    trjcat(unwrapped, pool_xtc, args.gmx, runner)

    # 3) Cluster.
    clusters_pdb, cluster_log = cluster_pool(
        pool_xtc, args.tpr, work_dir, args.rmsd_cutoff, args.gmx, runner)

    if args.dry_run:
        print("[dry-run] downstream centroid split + metrics require real cluster output; "
              "skipping. Re-run without --dry-run after cluster pipeline produces "
              "clusters.pdb + cluster log.")
        return 0

    # 4) Parse populations, then split top-K populated cluster centroids
    # → conformer_001..K.gro. Cluster IDs in GROMACS' log are not guaranteed
    # to be population sorted.
    populations = parse_cluster_log(cluster_log) if not args.dry_run else {}
    if not populations:
        print(f"ERROR: no cluster populations parsed from {cluster_log}",
              file=sys.stderr)
        return 4
    selected_cluster_ids = select_top_cluster_ids(populations, args.k)
    centroids = split_clusters_pdb_to_gros(
        clusters_pdb, out_dir, selected_cluster_ids, args.gmx, runner)
    try:
        require_centroid_count(centroids, args.k, clusters_pdb)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3

    # 5) Per-conformer metrics (Rg + E2E + population).
    pool_frames = sum(populations.values())
    asset_cluster_log = out_dir / "cluster_table.txt"
    if cluster_log.exists():
        shutil.copy2(cluster_log, asset_cluster_log)

    conformer_records: list[dict] = []
    for i, gro in enumerate(centroids, start=1):
        cluster_id = selected_cluster_ids[i - 1]
        rg = gyrate_rg_nm(gro, args.gmx, runner, work_dir)
        if args.dry_run or not gro.exists():
            e2e = float("nan")
        else:
            atoms = parse_gro(gro)
            e2e = compute_end_to_end_distance(atoms)
        pop = populations.get(cluster_id)
        conformer_records.append({
            "id": f"{i:03d}",
            "gro_path": str(gro.relative_to(out_dir)),
            "cluster_id": cluster_id,
            "population_frames": pop,
            "population_fraction": (pop / pool_frames) if (pop is not None
                                                            and pool_frames) else None,
            "rg_nm": None if rg != rg else round(rg, 4),  # NaN-safe
            "end_to_end_nm": None if e2e != e2e else round(e2e, 4),
            "gro_sha256": sha256_file(gro) if (not args.dry_run and gro.exists())
                          else "<dry-run>",
        })

    # 6) cluster_summary.yaml
    try:
        import yaml  # type: ignore
        yaml_loader = yaml
    except ImportError:
        yaml_loader = None

    summary = build_cluster_summary(
        pet_kind=args.pet_kind,
        k=args.k,
        rmsd_cutoff_nm=args.rmsd_cutoff,
        conformer_records=conformer_records,
        source_trajectories=[str(p) for p in args.trajectories],
        pool_frames=pool_frames,
        gromacs_version=gmx_version(args.gmx) if not args.dry_run else "<dry-run>",
        begin_ps=args.begin_ps,
        end_ps=args.end_ps,
    )
    summary_path = out_dir / "cluster_summary.yaml"
    if yaml_loader is not None:
        summary_path.write_text(yaml_loader.safe_dump(summary, sort_keys=False),
                                encoding="utf-8")
    else:
        # Fallback to JSON-flavoured YAML if PyYAML missing.
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # 7) MANIFEST.json — asset-layer SHA-256 of every output file.
    manifest_files = sorted(p for p in out_dir.iterdir()
                            if p.is_file() and p.name not in {"MANIFEST.json"})
    manifest = {
        "schema_version": "v1.1",
        "pet_kind": args.pet_kind,
        "files": {
            p.name: {
                "sha256": sha256_file(p) if not args.dry_run else "<dry-run>",
                "size_bytes": p.stat().st_size if p.exists() else None,
            }
            for p in manifest_files
        },
        "computed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    (out_dir / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    # 8) Clean intermediate unless --keep-pool.
    if not args.keep_pool and not args.dry_run and work_dir.exists():
        shutil.rmtree(work_dir, ignore_errors=True)

    print(f"[done] {args.pet_kind}: wrote K={args.k} conformers + summary → {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
