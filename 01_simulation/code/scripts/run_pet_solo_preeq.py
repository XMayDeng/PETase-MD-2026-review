#!/usr/bin/env python3
"""Run a single 50 ns PET-only pre-equilibration trajectory (R1 Phase A).

Spec: docs/plans/2026-05-03_md_workflow_v1_binding_attribution_spec.md §4.1
      "长链 dedicated preeq 协议 (R1)".

This script runs one independent trajectory for one chain length. The full
R1 Phase A asset is built by invoking this N=3 times per chain length, then
running extract_pet_conformers.py to pool and cluster.

Inputs are drawn from the frozen CHARMM-GUI raw asset
(inputs/pets/<kind>/raw/charmm-gui-001/gromacs/). Output is written to a
runtime directory in 01_preeq/ stage layout.

The NVT gen-seed and genion ion seed are derived from a deterministic
SHA-256 hash of (kind, run_idx); identical inputs always reproduce the
same trajectory.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = REPO_ROOT / "."
PREEQ_MDP_DIR = SIM_ROOT / "configurations/templates/preequilibration"
FF_DIR = SIM_ROOT / "inputs/forcefields/charmm36-jul2021.ff"


def _load_build_system_module():
    """Reuse build_system.write_pet_ff_add for PET FF deduplication.

    Loaded dynamically to avoid making this script a package import target
    (mirrors the pattern in build_md_from_docking_best.py).
    """
    path = Path(__file__).resolve().with_name("build_system.py")
    spec = importlib.util.spec_from_file_location("build_system", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module

SUPPORTED_KINDS = {"PET_L10", "PET_L20"}

# Chain-length-aware preeq box (triclinic, per spec §4.1 Phase A table).
# Unit: nm. Distinct from the production-MD cubic box (spec §4.5).
PREEQ_BOX_DIM_NM = {
    "PET_L10": (20.0, 10.0, 10.0),
    "PET_L20": (26.0, 10.0, 10.0),
}

# Default 50 ns MD at 2 fs time step = 25_000_000 steps. Overrides mdp's
# default 10M (20 ns) value via `mdrun -nsteps`.
DEFAULT_MD_NSTEPS = 25_000_000

# Genion ion concentration (target). v1 fingerprint pins to 0.10 M.
DEFAULT_ION_CONC_M = 0.10

# PET-only Phase A production timestep. Kept explicit because provenance
# stores duration in ns.
DEFAULT_MD_DT_FS = 2.0


# ---------------------------------------------------------------------------
# Pure helpers (unit-test target)
# ---------------------------------------------------------------------------

def replica_seed(pet_kind: str, run_idx: int) -> int:
    """Deterministic positive 31-bit gen-seed for NVT velocity assignment.

    GROMACS parses ``gen-seed`` and ``genion -seed`` as a signed 32-bit
    integer, so values must fit in ``[0, 2**31 - 1]``. Masking the SHA-256
    prefix to 31 bits guarantees this while preserving the deterministic
    namespace and 1 / 2^31 ~ 4.7e-10 collision probability across the
    ~hundred (case, replica) pairs in v1.

    Namespace ``pet_solo_preeq:`` keeps these distinct from production-MD
    replica seeds (see spec §4.6) and from genion ion seeds (below).
    """
    payload = f"pet_solo_preeq:{pet_kind}:r{run_idx}"
    return int(hashlib.sha256(payload.encode()).hexdigest()[:8], 16) & 0x7FFFFFFF


def ion_seed(pet_kind: str, run_idx: int) -> int:
    """Deterministic positive 31-bit seed for genion ion placement.

    Distinct namespace from NVT velocity seed so that velocity randomisation
    and ion-placement randomisation are independently auditable.
    """
    payload = f"pet_solo_preeq_ion:{pet_kind}:r{run_idx}"
    return int(hashlib.sha256(payload.encode()).hexdigest()[:8], 16) & 0x7FFFFFFF


def box_dim_for(pet_kind: str) -> tuple[float, float, float]:
    if pet_kind not in PREEQ_BOX_DIM_NM:
        raise ValueError(
            f"Unsupported pet_kind={pet_kind!r}. R1 Phase A applies only to "
            f"L10/L20 (short chains skip Phase A per spec §4.1)."
        )
    return PREEQ_BOX_DIM_NM[pet_kind]


def md_duration_ns(nsteps: int, dt_fs: float = DEFAULT_MD_DT_FS) -> float:
    """Convert MD steps × femtosecond timestep to nanoseconds."""
    return nsteps * dt_fs / 1_000_000.0


def rewrite_nvt_mdp_seed(template_text: str, gen_seed: int) -> str:
    """Return ``template_text`` with the gen-seed line replaced.

    Recognises ``gen-seed = <int>`` and ``gen_seed = <int>``; preserves
    surrounding whitespace and a trailing comment if present.
    """
    out_lines = []
    found = False
    for raw in template_text.splitlines(keepends=False):
        stripped = raw.lstrip()
        if stripped.startswith(("gen-seed", "gen_seed")):
            head, _, comment = raw.partition(";")
            indent_len = len(raw) - len(stripped)
            indent = raw[:indent_len]
            comment_suffix = f"  ;{comment}" if comment else ""
            out_lines.append(f"{indent}gen-seed                = {gen_seed}{comment_suffix}")
            found = True
        else:
            out_lines.append(raw)
    if not found:
        raise ValueError("Template NVT.mdp has no gen-seed line to override")
    return "\n".join(out_lines) + "\n"


# ---------------------------------------------------------------------------
# Subprocess helpers
# ---------------------------------------------------------------------------

class _Runner:
    def __init__(self, dry_run: bool, env: dict[str, str] | None = None):
        self.dry_run = dry_run
        self.env = env

    def run(self, cmd: list[str], *, cwd: Path | None = None,
            stdin: str | None = None) -> None:
        printable = " ".join(str(c) for c in cmd)
        prefix = "[dry-run] " if self.dry_run else "[run] "
        cwd_label = f"  (cwd={cwd})" if cwd else ""
        print(f"{prefix}{printable}{cwd_label}", flush=True)
        if self.dry_run:
            return
        subprocess.run(cmd, cwd=cwd, input=stdin, text=True, check=True,
                       env=self.env)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def setup_runtime(runtime_dir: Path) -> dict[str, Path]:
    preeq = runtime_dir / "01_preeq"
    dirs = {
        "preeq": preeq,
        "em": preeq / "EM",
        "nvt": preeq / "NVT",
        "npt": preeq / "NPT",
        "md": preeq / "MD",
        "outputs": preeq / "outputs",
        "qc": preeq / "qc",
        "inputs": preeq / "inputs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def sync_topology(src_dir: Path, dst_dir: Path) -> None:
    """Copy topol.top + toppar/ + pet_ff_add.itp (if present) from src to dst.

    Mirrors the bash sync_topology() helper used by the legacy preeq template.
    """
    shutil.copy2(src_dir / "topol.top", dst_dir / "topol.top")
    pet_ff_add = src_dir / "pet_ff_add.itp"
    if pet_ff_add.exists():
        shutil.copy2(pet_ff_add, dst_dir / "pet_ff_add.itp")
    dst_toppar = dst_dir / "toppar"
    if dst_toppar.exists():
        shutil.rmtree(dst_toppar)
    shutil.copytree(src_dir / "toppar", dst_toppar)


def stage_v1_topology(*, pet_kind: str, source_gromacs: Path, dst_dir: Path) -> None:
    """Generate a v1.1 PET-only topology in ``dst_dir``.

    The CHARMM-GUI raw asset's slim ``toppar/forcefield.itp`` only defines
    PET-specific parameters; it does not include water or ion moleculetypes,
    so a stock ``gmx solvate`` + ``gmx grompp`` chain fails with
    ``No such moleculetype SOL``. Following the v1.0 archive's working pattern,
    we replace the topology with one that:

      - pulls forcefield / TIP3P-CHARMM water / ions from
        ``inputs/forcefields/charmm36-jul2021.ff`` (the locked v1.1 FF tree),
      - keeps the CHARMM-GUI PET molecule topology (``toppar/P1.itp``),
      - adds a ``pet_ff_add.itp`` containing any PET-specific FF additions
        that are not already present in the locked FF tree (built via
        ``build_system.write_pet_ff_add``).

    Absolute paths are used for the FF includes; the topology is only consumed
    locally during R1 Phase A (asset build, one-time per chain length).
    """
    dst_toppar = dst_dir / "toppar"
    if dst_toppar.exists():
        shutil.rmtree(dst_toppar)
    dst_toppar.mkdir(parents=True)
    shutil.copy2(source_gromacs / "toppar" / "P1.itp", dst_toppar / "P1.itp")

    bs = _load_build_system_module()
    bs.write_pet_ff_add(
        src=source_gromacs / "toppar" / "forcefield.itp",
        dst=dst_dir / "pet_ff_add.itp",
        ffbonded=FF_DIR / "ffbonded.itp",
        ffnonbonded=FF_DIR / "ffnonbonded.itp",
    )

    ff_abs = FF_DIR.resolve()
    topol = (
        f'; {pet_kind} pre-equilibration topology (v1.1, R1 Phase A)\n'
        f'; generated by run_pet_solo_preeq.py\n'
        f'#include "{ff_abs}/forcefield.itp"\n'
        f'#include "pet_ff_add.itp"\n'
        f'#include "toppar/P1.itp"\n'
        f'#include "{ff_abs}/tip3p.itp"\n'
        f'#include "{ff_abs}/ions.itp"\n'
        f'\n'
        f'[ system ]\n'
        f'{pet_kind} in water (R1 Phase A)\n'
        f'\n'
        f'[ molecules ]\n'
        f'; Compound  #mols\n'
        f'P1          1\n'
    )
    (dst_dir / "topol.top").write_text(topol, encoding="utf-8")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def gmx_version(gmx_bin: Path) -> str:
    try:
        out = subprocess.run(
            [str(gmx_bin), "--version"], capture_output=True, text=True, check=True
        ).stdout
        for line in out.splitlines():
            if line.strip().startswith("GROMACS"):
                return line.strip().rstrip(":")
        return out.splitlines()[0].strip() if out else "unknown"
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def rg_convergence_check(qc_dir: Path, md_dir: Path, gmx_bin: Path,
                         runner: _Runner) -> dict[str, str | float]:
    """Run PET-only gmx gyrate on MD trajectory; emit a drift verdict."""
    md_traj = md_dir / "MD.xtc"
    if not md_traj.exists():
        if (md_dir / "MD.trr").exists():
            md_traj = md_dir / "MD.trr"
        else:
            return {"verdict": "NO_TRAJECTORY", "message": "neither MD.xtc nor MD.trr"}

    rg_xvg = qc_dir / "rg.xvg"
    runner.run([str(gmx_bin), "gyrate",
                "-f", str(md_traj), "-s", str(md_dir / "MD.tpr"),
                "-o", str(rg_xvg)],
               stdin="PETL\n")
    if runner.dry_run or not rg_xvg.exists():
        return {"verdict": "SKIPPED", "message": "dry-run or rg.xvg missing"}

    rgs: list[float] = []
    for line in rg_xvg.read_text().splitlines():
        if line.startswith(("#", "@")) or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                rgs.append(float(parts[1]))
            except ValueError:
                continue
    if len(rgs) < 20:
        return {"verdict": "INSUFFICIENT_DATA",
                "message": f"only {len(rgs)} frames",
                "rg_n_frames": len(rgs)}
    mid = len(rgs) // 2
    r1 = sum(rgs[:mid]) / mid
    r2 = sum(rgs[mid:]) / (len(rgs) - mid)
    drift_pct = abs(r2 - r1) / r1 * 100.0
    verdict = "FAIL" if drift_pct > 10 else ("WARN" if drift_pct > 5 else "PASS")
    (qc_dir / "pet_rg_summary.txt").write_text(
        f"verdict={verdict}\nRg_first={r1:.4f} Rg_last={r2:.4f} drift={drift_pct:.2f}%\n",
        encoding="utf-8")
    return {"verdict": verdict, "rg_first_half_nm": r1,
            "rg_second_half_nm": r2, "drift_pct": drift_pct}


def write_provenance(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pet-kind", required=True, choices=sorted(SUPPORTED_KINDS),
                   help="Chain length (PET_L10 or PET_L20; L4/L2 skip Phase A per spec §4.1).")
    p.add_argument("--run-idx", required=True, type=int,
                   help="Trajectory replica index 1..N (default N=3 in wrapper).")
    p.add_argument("--runtime-dir", required=True, type=Path,
                   help="Output directory; 01_preeq/ stage tree is created beneath it.")
    p.add_argument("--md-nsteps", type=int, default=DEFAULT_MD_NSTEPS,
                   help=f"Production MD nsteps (default {DEFAULT_MD_NSTEPS} = 50 ns @ 2 fs).")
    p.add_argument("--ion-conc", type=float, default=DEFAULT_ION_CONC_M,
                   help=f"NaCl target concentration in M (default {DEFAULT_ION_CONC_M}).")
    p.add_argument("--gmx", type=Path, default=Path(os.environ.get("GMX_BIN", "/usr/local/gromacs/bin/gmx")),
                   help="Path to gmx binary.")
    p.add_argument("--gpu-id", type=str, default=os.environ.get("GPU_ID", ""),
                   help="GPU id for mdrun -gpu_id; empty string disables.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    pet_kind = args.pet_kind
    gen_seed = replica_seed(pet_kind, args.run_idx)
    ion_seed_val = ion_seed(pet_kind, args.run_idx)
    box_dim = box_dim_for(pet_kind)

    source_gromacs = SIM_ROOT / f"inputs/pets/{pet_kind}/raw/charmm-gui-001/gromacs"
    source_gro = source_gromacs / "step3_input.gro"
    if not source_gro.exists():
        print(f"ERROR: source asset missing: {source_gro}", file=sys.stderr)
        return 2
    if not PREEQ_MDP_DIR.is_dir():
        print(f"ERROR: preeq mdp dir missing: {PREEQ_MDP_DIR}", file=sys.stderr)
        return 2

    runtime_dir = args.runtime_dir.resolve()
    dirs = setup_runtime(runtime_dir)
    runner = _Runner(dry_run=args.dry_run)

    # Stage inputs as a v1.1 topology snapshot (raw CHARMM-GUI topology is
    # not solvate-ready — see stage_v1_topology docstring).
    stage_v1_topology(pet_kind=pet_kind,
                      source_gromacs=source_gromacs,
                      dst_dir=dirs["inputs"])
    shutil.copy2(source_gro, dirs["inputs"] / "step3_input.gro")

    # ----- EM -----
    sync_topology(dirs["inputs"], dirs["em"])
    shutil.copy2(source_gro, dirs["em"] / "pet0.gro")
    gmx = args.gmx

    runner.run([str(gmx), "editconf",
                "-f", "pet0.gro", "-o", "pet_box.gro",
                "-c", "-box", *(f"{d:.2f}" for d in box_dim),
                "-bt", "triclinic"],
               cwd=dirs["em"])
    runner.run([str(gmx), "solvate",
                "-cp", "pet_box.gro", "-cs", "spc216.gro",
                "-p", "topol.top", "-o", "pet_solv.gro"],
               cwd=dirs["em"])
    runner.run([str(gmx), "grompp",
                "-f", str(PREEQ_MDP_DIR / "em_preeq.mdp"),
                "-c", "pet_solv.gro", "-p", "topol.top",
                "-o", "ions.tpr", "-maxwarn", "20"],
               cwd=dirs["em"])
    runner.run([str(gmx), "genion",
                "-s", "ions.tpr", "-o", "pet_solv_ion.gro",
                "-p", "topol.top", "-neutral",
                "-conc", f"{args.ion_conc}",
                "-pname", "SOD", "-nname", "CLA",
                "-seed", str(ion_seed_val)],
               cwd=dirs["em"], stdin="SOL\n")
    runner.run([str(gmx), "grompp",
                "-f", str(PREEQ_MDP_DIR / "em_preeq.mdp"),
                "-c", "pet_solv_ion.gro", "-p", "topol.top",
                "-o", "EM.tpr", "-maxwarn", "20"],
               cwd=dirs["em"])
    em_cmd = [str(gmx), "mdrun", "-deffnm", "EM"]
    if args.gpu_id:
        em_cmd.extend(["-gpu_id", args.gpu_id])
    runner.run(em_cmd, cwd=dirs["em"])

    # ----- NVT (gen-seed override via temp mdp) -----
    sync_topology(dirs["em"], dirs["nvt"])
    nvt_template = (PREEQ_MDP_DIR / "nvt_preeq.mdp").read_text(encoding="utf-8")
    nvt_mdp_text = rewrite_nvt_mdp_seed(nvt_template, gen_seed)
    nvt_mdp_path = dirs["nvt"] / "nvt_preeq.mdp"
    nvt_mdp_path.write_text(nvt_mdp_text, encoding="utf-8")

    runner.run([str(gmx), "grompp",
                "-f", "nvt_preeq.mdp",
                "-c", str(dirs["em"] / "EM.gro"),
                "-r", str(dirs["em"] / "EM.gro"),
                "-p", "topol.top",
                "-o", "NVT.tpr", "-maxwarn", "20"],
               cwd=dirs["nvt"])
    nvt_cmd = [str(gmx), "mdrun", "-deffnm", "NVT"]
    if args.gpu_id:
        nvt_cmd.extend(["-gpu_id", args.gpu_id])
    runner.run(nvt_cmd, cwd=dirs["nvt"])

    # ----- NPT -----
    sync_topology(dirs["em"], dirs["npt"])
    runner.run([str(gmx), "grompp",
                "-f", str(PREEQ_MDP_DIR / "npt_preeq.mdp"),
                "-c", str(dirs["nvt"] / "NVT.gro"),
                "-r", str(dirs["nvt"] / "NVT.gro"),
                "-t", str(dirs["nvt"] / "NVT.cpt"),
                "-p", "topol.top",
                "-o", "NPT.tpr", "-maxwarn", "20"],
               cwd=dirs["npt"])
    npt_cmd = [str(gmx), "mdrun", "-deffnm", "NPT"]
    if args.gpu_id:
        npt_cmd.extend(["-gpu_id", args.gpu_id])
    runner.run(npt_cmd, cwd=dirs["npt"])

    # ----- Production MD (50 ns by default) -----
    sync_topology(dirs["em"], dirs["md"])
    runner.run([str(gmx), "grompp",
                "-f", str(PREEQ_MDP_DIR / "md_preeq.mdp"),
                "-c", str(dirs["npt"] / "NPT.gro"),
                "-t", str(dirs["npt"] / "NPT.cpt"),
                "-p", "topol.top",
                "-o", "MD.tpr", "-maxwarn", "20"],
               cwd=dirs["md"])
    md_cmd = [str(gmx), "mdrun", "-deffnm", "MD", "-nsteps", str(args.md_nsteps)]
    if args.gpu_id:
        md_cmd.extend(["-gpu_id", args.gpu_id])
    runner.run(md_cmd, cwd=dirs["md"])

    # ----- QC -----
    rg_summary = rg_convergence_check(dirs["qc"], dirs["md"], gmx, runner)

    # ----- Provenance -----
    provenance = {
        "schema_version": "v1.1",
        "stage": "R1_phase_A_pet_solo_preeq",
        "pet_kind": pet_kind,
        "run_idx": args.run_idx,
        "gen_seed_nvt": gen_seed,
        "ion_seed_genion": ion_seed_val,
        "box_dim_nm": list(box_dim),
        "box_type": "triclinic",
        "md_nsteps": args.md_nsteps,
        "md_dt_fs": DEFAULT_MD_DT_FS,
        "md_duration_ns": md_duration_ns(args.md_nsteps, DEFAULT_MD_DT_FS),
        "ion_target_concentration_M": args.ion_conc,
        "source_gro_sha256": sha256_file(source_gro) if not args.dry_run else "<dry-run>",
        "source_topol_sha256": (
            sha256_file(source_gromacs / "topol.top") if not args.dry_run else "<dry-run>"
        ),
        "gromacs_version": gmx_version(gmx) if not args.dry_run else "<dry-run>",
        "rg_qc": rg_summary,
        "computed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    write_provenance(dirs["preeq"] / "provenance.json", provenance)

    # Stage outputs (canonical asset hand-off).
    if not args.dry_run and (dirs["md"] / "MD.gro").exists():
        shutil.copy2(dirs["md"] / "MD.gro",
                     dirs["outputs"] / "pet_md_final.gro")

    print(f"[done] {pet_kind} r{args.run_idx} → {dirs['preeq']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
