#!/usr/bin/env python3
"""v1.1 case orchestrator for binding-attribution MD cases.

Drives a single case.yaml v1.1 through the full pipeline:

    materialize → docking → align → select → build → validate_box
        → for each replica: NVT → NPT → production MD
            → postprocess → per-replica gate → fingerprint
        → aggregate gate

Long-chain (L10/L20) cases add R1 Phase B multi-conformer screening,
per-direction build, and Phase C chain_relax before replica production.

Usage::

    python run_case_v1_1.py <case.yaml> \\
        [--runtime-root <path>] \\
        [--from-stage <name>] \\
        [--to-stage <name>] \\
        [--dry-run] \\
        [--gmx /usr/local/gromacs/bin/gmx] \\
        [--gpu-id 0]
"""
from __future__ import annotations

import argparse
import dataclasses
import datetime
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Callable

REPO_ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = REPO_ROOT / "."
SCRIPTS_DIR = Path(__file__).resolve().parent
MDP_PREEQ_DIR = SIM_ROOT / "configurations/templates/preequilibration"
MDP_BUILD_DIR = SIM_ROOT / "configurations/templates/md"
MDP_CHAIN_RELAX_DIR = SIM_ROOT / "configurations/templates/chain_relaxation"
DEFAULT_PRODUCTION_BURN_IN_NS = 20.0


# ---------------------------------------------------------------------------
# Dynamic-load helpers (mirrors the run_pet_solo_preeq pattern)
# ---------------------------------------------------------------------------

def _load_module(rel_path: str, name: str):
    path = SIM_ROOT / rel_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_mat():
    return _load_module("code/pipelines/materialize_v1_1.py", "materialize_v1_1")


def _load_box():
    return _load_module("code/pipelines/box_protocol.py", "box_protocol")


def deterministic_seed(*parts: object) -> int:
    """Stable positive 31-bit seed for tools requiring signed int seeds."""
    payload = ":".join(str(p) for p in parts)
    return int(hashlib.sha256(payload.encode()).hexdigest()[:8], 16) & 0x7FFFFFFF


def _default_vina_env() -> Path:
    return Path(os.environ.get("AUTOPROPET_VINA_ENV") or
                REPO_ROOT / "extra_env/autodock_vina_2026.04/vina")


def _vina_env_tool(name: str) -> Path | None:
    candidate = _default_vina_env() / "bin" / name
    return candidate if candidate.exists() else None


# ---------------------------------------------------------------------------
# Stage / context dataclasses
# ---------------------------------------------------------------------------

PIPELINE_STAGES_SHORT_CHAIN: tuple[str, ...] = (
    "materialize",
    "docking",
    "align",
    "select",
    "build",
    "validate_box",
    "em",
    "replica_md",
    "postprocess_completed",
    "postprocess",
    "gate",
    "aggregate_gate",
    "fingerprint",
)

PIPELINE_STAGES_LONG_CHAIN: tuple[str, ...] = (
    "materialize",
    "docking",
    "screen_r1_phase_b",
    "select_per_direction",
    "build_per_direction",
    "validate_box",
    "chain_relax",
    "replica_md",
    "postprocess_completed",
    "postprocess",
    "gate",
    "aggregate_gate",
    "fingerprint",
)


@dataclasses.dataclass
class CaseContext:
    case_yaml: Path
    case_data: dict
    case_id: str
    pet_kind: str
    is_long_chain: bool
    case_dir: Path
    replica_count: int
    directions: list[str]
    box_dim_nm: float
    chain_relax_enabled: bool
    gmx_bin: Path
    gpu_id: str
    dry_run: bool


@dataclasses.dataclass
class StageResult:
    name: str
    status: str            # "done" / "skipped" / "failed"
    detail: str = ""
    artifacts: list[str] = dataclasses.field(default_factory=list)


# ---------------------------------------------------------------------------
# Subprocess wrapper + summary I/O
# ---------------------------------------------------------------------------

def _run(cmd: list[str], *, cwd: Path | None = None, stdin: str | None = None,
         dry_run: bool, log_path: Path | None = None) -> int:
    printable = " ".join(str(c) for c in cmd)
    prefix = "[dry-run] " if dry_run else "[run] "
    cwd_label = f"  (cwd={cwd})" if cwd else ""
    line = f"{prefix}{printable}{cwd_label}"
    print(line, flush=True)
    if dry_run:
        return 0
    log_handle = log_path.open("a", encoding="utf-8") if log_path else None
    try:
        if log_handle:
            log_handle.write(f"\n=== {datetime.datetime.now(datetime.timezone.utc).isoformat()} ===\n")
            log_handle.write(line + "\n")
            log_handle.flush()
        proc = subprocess.run(cmd, cwd=cwd, input=stdin, text=True,
                              check=False,
                              stdout=log_handle if log_handle else None,
                              stderr=log_handle if log_handle else None)
        return proc.returncode
    finally:
        if log_handle:
            log_handle.close()


def _summary_path(case_dir: Path) -> Path:
    return case_dir / "qc" / "orchestration_summary.json"


def _load_summary(case_dir: Path) -> dict:
    p = _summary_path(case_dir)
    if not p.exists():
        return {"stages": {}}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {"stages": {}}


def _save_summary(case_dir: Path, summary: dict) -> None:
    p = _summary_path(case_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    summary["updated_at_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    p.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")


# ---------------------------------------------------------------------------
# Stage executors
# ---------------------------------------------------------------------------

def stage_materialize(ctx: CaseContext) -> StageResult:
    mat = _load_mat()
    try:
        mat.assert_conformer_asset_present(ctx.pet_kind)
    except FileNotFoundError as exc:
        return StageResult(
            "materialize", "failed",
            detail=f"long-chain conformer asset missing: {exc}",
        )
    if not ctx.dry_run:
        manifest = mat.materialize_v1_1_subruntime_tree(
            case_dir=ctx.case_dir,
            pet_kind_value=ctx.pet_kind,
            directions_list=ctx.directions,
            replica_count_value=ctx.replica_count,
            chain_relax=ctx.chain_relax_enabled,
        )
        nvt_template = MDP_BUILD_DIR / "NVT.mdp"
        if nvt_template.exists():
            mat.materialize_v1_1_nvt_mdps(
                case_dir=ctx.case_dir,
                case_id_value=ctx.case_id,
                pet_kind_value=ctx.pet_kind,
                directions_list=ctx.directions,
                replica_count_value=ctx.replica_count,
                nvt_template_path=nvt_template,
            )
        return StageResult("materialize", "done",
                           detail=f"{len(manifest['sub_runtimes'])} sub-runtime",
                           artifacts=[str(ctx.case_dir)])
    return StageResult("materialize", "done", detail="dry-run")


def _docking_fragment_kind(ctx: CaseContext) -> str:
    docking_cfg = ctx.case_data.get("docking") or {}
    explicit = docking_cfg.get("fragment")
    if explicit:
        return str(explicit)
    return "PET_L4" if ctx.is_long_chain else ctx.pet_kind


def _docking_seed(ctx: CaseContext, fragment_kind: str) -> int:
    docking_cfg = ctx.case_data.get("docking") or {}
    if docking_cfg.get("seed") is not None:
        return int(docking_cfg["seed"])
    return deterministic_seed("v1.1", "vina", ctx.case_id, fragment_kind)


def _build_genion_seed(ctx: CaseContext, direction: str | None = None) -> int:
    build_cfg = ctx.case_data.get("build") or {}
    raw = build_cfg.get("genion_seed")
    if isinstance(raw, dict):
        if direction and direction in raw:
            return int(raw[direction])
        if "default" in raw:
            return int(raw["default"])
    elif raw is not None:
        return int(raw)
    return deterministic_seed("v1.1", "genion", ctx.case_id, direction or "shared")


def pet_unit_map_path(pet_kind: str) -> Path:
    suffix = pet_kind.lower().replace("pet_", "")
    return SIM_ROOT / f"inputs/pets/{pet_kind}/alignment/pet_{suffix}_unit_map.yaml"


def ensure_pet_unit_map(ctx: CaseContext, *, log_path: Path) -> tuple[Path, str | None]:
    """Return the PET alignment map path, generating it if a long-chain map
    is missing.

    L10/L20 screen and chain_relax both require the same repeat-unit map.
    Keeping the check here prevents a case from failing halfway through the
    expensive long-chain workflow because the asset was never materialized.
    """
    unit_map = pet_unit_map_path(ctx.pet_kind)
    if unit_map.exists() or not ctx.is_long_chain:
        return unit_map, None

    script = SCRIPTS_DIR / "build_pet_alignment_maps.py"
    asset_root = SIM_ROOT / f"inputs/pets/{ctx.pet_kind}"
    source_pdb = asset_root / "raw/charmm-gui-001/gromacs/step3_input.pdb"
    if not ctx.dry_run and not source_pdb.exists():
        return unit_map, f"missing PET source PDB for unit map: {source_pdb}"

    log_path.parent.mkdir(parents=True, exist_ok=True)
    rc = _run([sys.executable, str(script),
               "--pet-kind", ctx.pet_kind,
               "--asset-root", str(asset_root)],
              dry_run=ctx.dry_run, log_path=log_path)
    if rc != 0 and not ctx.dry_run:
        return unit_map, f"unit map generation rc={rc}; see {log_path}"
    if not ctx.dry_run and not unit_map.exists():
        return unit_map, f"unit map generation did not create {unit_map}"
    return unit_map, None


def stage_docking(ctx: CaseContext) -> StageResult:
    """Vina docking of the PET fragment into the active-site grid.

    Delegates to ``run_vina_pet_fragment_docking.py``. Required inputs are
    resolved from case.yaml: protein PDB from
    ``inputs/proteins/<id>/source/<id>.pdb`` and ligand PDB from the
    fragment asset (L2/L4 dock directly; L10/L20 default to the PET_L4
    anchor unless case.yaml overrides ``docking.fragment``).
    Active-site grid centre comes from the registry via ``--site-source
    registry``.
    """
    script = SCRIPTS_DIR / "run_vina_pet_fragment_docking.py"
    if not script.exists():
        return StageResult("docking", "failed",
                           detail=f"missing docking script: {script}")
    out_dir = ctx.case_dir / "02_docking"
    log = ctx.case_dir / "logs" / "docking.log"

    protein_cfg = ctx.case_data.get("protein") or {}
    protein_id = str(protein_cfg.get("id") or "00000")
    protein_pdb = SIM_ROOT / (protein_cfg.get("asset")
                              or f"inputs/proteins/{protein_id}") \
                          / "source" / f"{protein_id}.pdb"
    fragment_kind = _docking_fragment_kind(ctx)
    ligand_pdb = (SIM_ROOT /
                  f"inputs/pets/{fragment_kind}/raw/charmm-gui-001/"
                  f"gromacs/step3_input.pdb")
    docking_cfg = ctx.case_data.get("docking") or {}
    max_clash = int(docking_cfg.get("max_severe_clashes", 50))
    seed = _docking_seed(ctx, fragment_kind)

    cmd = [sys.executable, str(script),
           "--protein-id", protein_id,
           "--protein-pdb", str(protein_pdb),
           "--ligand-pdb", str(ligand_pdb),
           "--out-dir", str(out_dir),
           "--site-source", "registry",
           "--max-severe-clashes", str(max_clash),
           "--seed", str(seed),
           "--force"]
    if docking_cfg.get("require_trusted_active_site", True):
        cmd.append("--require-trusted-active-site")
    # AUTOPROPET_VINA_ENV points to the conda env directory containing
    # vina + obabel + mk_prepare_{ligand,receptor}.py under bin/. Default
    # to the in-tree extra_env install when the env var is unset.
    for flag, name in (("--vina-bin", "vina"),
                       ("--obabel-bin", "obabel"),
                       ("--prepare-ligand-bin", "mk_prepare_ligand.py"),
                       ("--prepare-receptor-bin", "mk_prepare_receptor.py")):
        candidate = _vina_env_tool(name)
        if candidate is not None:
            cmd.extend([flag, str(candidate)])

    rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("docking", "failed", detail=f"rc={rc}; see {log}")
    return StageResult("docking", "done", artifacts=[str(out_dir)])


def stage_align(ctx: CaseContext) -> StageResult:
    """Single-conformer Kabsch alignment screen for short-chain L4/L2.

    Long-chain L10/L20 uses ``stage_screen_r1_phase_b`` instead.
    """
    if ctx.is_long_chain:
        return StageResult("align", "skipped",
                           detail="long-chain uses screen_r1_phase_b")
    script = SCRIPTS_DIR / "screen_vina_pet_longchain_candidates.py"
    out_dir = ctx.case_dir / "03_pet_alignment_screen"
    log = ctx.case_dir / "logs" / "align.log"
    cmd = [sys.executable, str(script),
           "--source-pet-kind", ctx.pet_kind,
           "--target-pet-kind", ctx.pet_kind,
           "--docking-dir", str(ctx.case_dir / "02_docking"),
           "--out-dir", str(out_dir),
           "--force"]
    rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("align", "failed", detail=f"rc={rc}; see {log}")
    return StageResult("align", "done", artifacts=[str(out_dir)])


def stage_select(ctx: CaseContext) -> StageResult:
    """Pick the v1 best clash-free, motif-aligned candidate."""
    script = SCRIPTS_DIR / "select_v1_binding_candidate.py"
    screen_dir = ctx.case_dir / "03_pet_alignment_screen"
    target = ctx.case_dir / "02_docking" / "best_rigidbody.pdb"
    meta = ctx.case_dir / "qc" / "v1_candidate_selection.json"
    log = ctx.case_dir / "logs" / "select.log"
    cmd = [sys.executable, str(script),
           "--screen-dir", str(screen_dir),
           "--target-pdb", str(target),
           "--meta", str(meta)]
    rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("select", "failed", detail=f"rc={rc}; see {log}")
    return StageResult("select", "done", artifacts=[str(target), str(meta)])


def stage_build(ctx: CaseContext) -> StageResult:
    """CHARMM topology + v1.1 chain-length-aware cubic box."""
    script = SCRIPTS_DIR / "build_md_from_docking_best.py"
    log = ctx.case_dir / "logs" / "build.log"
    best_model = ctx.case_dir / "02_docking" / "best_rigidbody.pdb"
    build_dir = ctx.case_dir / "03_md_build"
    protein_id = str((ctx.case_data.get("protein") or {}).get("id") or "00000")
    cmd = [sys.executable, str(script),
           "--best-model", str(best_model),
           "--build-dir", str(build_dir),
           "--pet-kind", ctx.pet_kind,
           "--box-dim-nm", str(ctx.box_dim_nm),
           "--protein-id", protein_id,
           "--gmx-bin", str(ctx.gmx_bin),
           "--genion-seed", str(_build_genion_seed(ctx))]
    rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("build", "failed", detail=f"rc={rc}; see {log}")
    return StageResult("build", "done", artifacts=[str(build_dir)])


def _preferred_build_gro(build_outputs: Path) -> Path:
    for name in ("box_sol_ion.gro", "box.gro"):
        candidate = build_outputs / name
        if candidate.exists():
            return candidate
    candidates = sorted(build_outputs.glob("*.gro")) if build_outputs.exists() else []
    return candidates[0] if candidates else build_outputs / "box_sol_ion.gro"


def stage_validate_box(ctx: CaseContext) -> StageResult:
    """Reject build outputs whose solute envelope is within 1.5 nm of any
    box face (the build stage handles fallback; this stage only verifies)."""
    script = SCRIPTS_DIR / "validate_complex_box.py"
    log = ctx.case_dir / "logs" / "validate_box.log"
    targets = (
        [(direction, _build_dir(ctx, direction)) for direction in ctx.directions]
        if ctx.is_long_chain
        else [(None, _build_dir(ctx, None))]
    )
    artifacts: list[str] = []
    for direction, build_dir in targets:
        build_outputs = build_dir / "outputs"
        primary = _preferred_build_gro(build_outputs)
        label = direction or "shared"
        if not primary.exists() and not ctx.dry_run:
            return StageResult("validate_box", "failed",
                               detail=f"{label}: no .gro under {build_outputs}")
        out_name = "validate_complex_box.json"
        if direction:
            out_name = f"validate_complex_box_{direction}.json"
        out = ctx.case_dir / "qc" / out_name
        cmd = [sys.executable, str(script),
               "--gro", str(primary),
               "--pet-kind", ctx.pet_kind,
               "--strict-cubic",
               "--out", str(out)]
        rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("validate_box", "failed",
                               detail=f"{label}: rc={rc}; see {log}")
        artifacts.append(str(out))
    return StageResult("validate_box", "done",
                       detail=f"{len(targets)} build(s)",
                       artifacts=artifacts)


def _direction_dirname(direction: str) -> str:
    return f"dir_{direction}"


def _replica_runtime_dir(ctx: CaseContext, direction: str | None, idx: int) -> Path:
    """Match materialize_v1_1.materialize_v1_1_subruntime_tree layout."""
    if direction:
        return ctx.case_dir / "05_md_prod" / _direction_dirname(direction) / f"r{idx}"
    return ctx.case_dir / "05_md_prod" / f"r{idx}"


def production_burn_in_ns(ctx: CaseContext) -> float:
    md_cfg = ctx.case_data.get("md") or {}
    return float(md_cfg.get("burn_in_ns", DEFAULT_PRODUCTION_BURN_IN_NS))


def _postprocess_runtime_dir(ctx: CaseContext, direction: str | None, idx: int) -> Path:
    if direction:
        return ctx.case_dir / "06_postprocess" / _direction_dirname(direction) / f"r{idx}"
    return ctx.case_dir / "06_postprocess" / f"r{idx}"


def _md_output_dir(rdir: Path) -> Path:
    md_dir = rdir / "MD"
    if (md_dir / "MD.xtc").exists() and (md_dir / "MD.tpr").exists():
        return md_dir
    return rdir


def _iter_sub_runtimes(ctx: CaseContext):
    """Yield (direction, replica_idx) tuples covering this case.

    Long-chain: 3 directions x N replicas.
    Short-chain: direction=None x N replicas.
    """
    if ctx.is_long_chain:
        for direction in ctx.directions:
            for rep in range(1, ctx.replica_count + 1):
                yield direction, rep
    else:
        for rep in range(1, ctx.replica_count + 1):
            yield None, rep


def _verdict_filename(direction: str | None, idx: int) -> str:
    if direction:
        return f"v1_production_gate_verdict_{direction}_r{idx}.json"
    return f"v1_production_gate_verdict_r{idx}.json"


def _fingerprint_filename(direction: str | None, idx: int) -> str:
    if direction:
        return f"protocol_fingerprint_{direction}_r{idx}.json"
    return f"protocol_fingerprint_r{idx}.json"


def _build_dir(ctx: CaseContext, direction: str | None) -> Path:
    """Per-direction build output (long-chain) or shared case-level build (short)."""
    if direction:
        return ctx.case_dir / "03_md_build" / _direction_dirname(direction)
    return ctx.case_dir / "03_md_build"


def _shared_em_gro(ctx: CaseContext, direction: str | None) -> Path:
    if direction:
        return ctx.case_dir / "05_md_prod" / _direction_dirname(direction) / "EM" / "EM.gro"
    return ctx.case_dir / "05_md_prod" / "EM" / "EM.gro"


def ensure_pet_anchor_posres_include(topol_path: Path,
                                     include_name: str = "posre_pet_anchor.itp") -> None:
    """Ensure PET anchor POSRES is included in the PET moleculetype."""
    include_line = f'#include "{include_name}"'
    lines = topol_path.read_text(encoding="utf-8").splitlines()
    if include_line in lines:
        return
    block = ["#ifdef POSRES", include_line, "#endif"]
    insert_idx = None
    for i, line in enumerate(lines):
        if line.strip() == '#include "P1.itp"':
            insert_idx = i + 1
            break
    if insert_idx is None:
        for i, line in enumerate(lines):
            if line.strip().startswith("[ system ]"):
                insert_idx = i
                break
    if insert_idx is None:
        insert_idx = len(lines)
    lines[insert_idx:insert_idx] = block
    topol_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_grompp_mdrun_chain(ctx: CaseContext, rdir: Path, *,
                            label: str, npt_mdp: Path, md_mdp: Path,
                            topol_path: Path, index_path: Path,
                            log_path: Path) -> tuple[int, str]:
    """Run NVT → NPT → production MD under ``rdir/{NVT,NPT,MD}``.

    Assumes:
      * ``rdir/NVT.mdp`` already materialised (replica-specific gen-seed)
      * ``rdir/../EM/EM.gro`` exists (shared EM upstream)
      * ``topol_path`` absolute (build outputs dir holds .itp includes)
      * ``index_path`` absolute (provides ``Protein_LIG`` / ``Water_and_ions``
        tc-grps referenced by NVT/NPT/MD.mdp)

    Returns ``(rc, stage)`` where ``stage`` names the failing step or
    ``"done"`` on full success.
    """
    em_gro = rdir.parent / "EM" / "EM.gro"
    nvt_mdp = rdir / "NVT.mdp"
    nvt_dir = rdir / "NVT"
    npt_dir = rdir / "NPT"
    md_dir = rdir / "MD"
    for subdir in (nvt_dir, npt_dir, md_dir):
        subdir.mkdir(parents=True, exist_ok=True)
    if not ctx.dry_run:
        if not em_gro.exists():
            return 2, f"{label}: EM not done ({em_gro})"
        if not nvt_mdp.exists():
            return 2, f"{label}: NVT.mdp missing in {rdir}"
        if not topol_path.exists():
            return 2, f"{label}: topol.top missing ({topol_path})"
        if not index_path.exists():
            return 2, f"{label}: index.ndx missing ({index_path})"

    def _complete(phase_dir: Path, phase: str) -> bool:
        if ctx.dry_run:
            return False
        required = [phase_dir / f"{phase}.gro"]
        if phase in ("NVT", "NPT"):
            required.append(phase_dir / f"{phase}.cpt")
        if phase == "MD":
            required.extend([
                phase_dir / "MD.xtc",
                phase_dir / "MD.edr",
            ])
        return all(p.exists() for p in required)

    def _mdrun_cmd(phase_dir: Path, phase: str) -> list[str]:
        cmd = [str(ctx.gmx_bin), "mdrun", "-deffnm", phase]
        if not ctx.dry_run:
            cpt = phase_dir / f"{phase}.cpt"
            if cpt.exists() and not _complete(phase_dir, phase):
                cmd.extend(["-cpi", str(cpt), "-append"])
        if ctx.gpu_id:
            cmd.extend(["-gpu_id", ctx.gpu_id])
        return cmd

    # NVT
    if not _complete(nvt_dir, "NVT"):
        if ctx.dry_run or not (nvt_dir / "NVT.tpr").exists():
            rc = _run([str(ctx.gmx_bin), "grompp",
                       "-f", str(nvt_mdp),
                       "-c", str(em_gro), "-r", str(em_gro),
                       "-p", str(topol_path),
                       "-n", str(index_path),
                       "-o", "NVT.tpr", "-maxwarn", "20"],
                      cwd=nvt_dir, dry_run=ctx.dry_run, log_path=log_path)
            if rc != 0 and not ctx.dry_run:
                return rc, f"{label}: NVT grompp"
        rc = _run(_mdrun_cmd(nvt_dir, "NVT"),
                  cwd=nvt_dir, dry_run=ctx.dry_run, log_path=log_path)
        if rc != 0 and not ctx.dry_run:
            return rc, f"{label}: NVT mdrun"

    # NPT
    if not _complete(npt_dir, "NPT"):
        if ctx.dry_run or not (npt_dir / "NPT.tpr").exists():
            rc = _run([str(ctx.gmx_bin), "grompp",
                       "-f", str(npt_mdp),
                       "-c", str(nvt_dir / "NVT.gro"), "-r", str(nvt_dir / "NVT.gro"),
                       "-t", str(nvt_dir / "NVT.cpt"),
                       "-p", str(topol_path),
                       "-n", str(index_path),
                       "-o", "NPT.tpr", "-maxwarn", "20"],
                      cwd=npt_dir, dry_run=ctx.dry_run, log_path=log_path)
            if rc != 0 and not ctx.dry_run:
                return rc, f"{label}: NPT grompp"
        rc = _run(_mdrun_cmd(npt_dir, "NPT"),
                  cwd=npt_dir, dry_run=ctx.dry_run, log_path=log_path)
        if rc != 0 and not ctx.dry_run:
            return rc, f"{label}: NPT mdrun"

    # Production MD
    if not _complete(md_dir, "MD"):
        if ctx.dry_run or not (md_dir / "MD.tpr").exists():
            rc = _run([str(ctx.gmx_bin), "grompp",
                       "-f", str(md_mdp),
                       "-c", str(npt_dir / "NPT.gro"),
                       "-t", str(npt_dir / "NPT.cpt"),
                       "-p", str(topol_path),
                       "-n", str(index_path),
                       "-o", "MD.tpr", "-maxwarn", "20"],
                      cwd=md_dir, dry_run=ctx.dry_run, log_path=log_path)
            if rc != 0 and not ctx.dry_run:
                return rc, f"{label}: MD grompp"
        rc = _run(_mdrun_cmd(md_dir, "MD"),
                  cwd=md_dir, dry_run=ctx.dry_run, log_path=log_path)
        if rc != 0 and not ctx.dry_run:
            return rc, f"{label}: MD mdrun"
    return 0, "done"


def stage_em(ctx: CaseContext) -> StageResult:
    """Short-chain shared energy minimisation (no POSRES).

    Produces ``05_md_prod/EM/EM.gro`` from the build's solvated/ionised
    box. Long-chain cases reach the same artifact via ``chain_relax`` and
    must not invoke this stage.
    """
    if ctx.is_long_chain:
        return StageResult("em", "skipped",
                           detail="long-chain uses chain_relax")
    em_mdp = MDP_BUILD_DIR / "EM.mdp"
    build_outputs = _build_dir(ctx, None) / "outputs"
    box_gro = build_outputs / "box_sol_ion.gro"
    topol = build_outputs / "topol.top"
    em_dir = _shared_em_gro(ctx, None).parent
    em_dir.mkdir(parents=True, exist_ok=True)
    log = ctx.case_dir / "logs" / "em.log"

    if not ctx.dry_run:
        if not box_gro.exists():
            return StageResult("em", "failed",
                               detail=f"missing build box: {box_gro}")
        if not topol.exists():
            return StageResult("em", "failed",
                               detail=f"missing topol: {topol}")

    rc = _run([str(ctx.gmx_bin), "grompp",
               "-f", str(em_mdp),
               "-c", str(box_gro),
               "-p", str(topol),
               "-o", str(em_dir / "EM.tpr"), "-maxwarn", "20"],
              cwd=em_dir, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("em", "failed", detail=f"EM grompp rc={rc}")
    cmd = [str(ctx.gmx_bin), "mdrun", "-deffnm", "EM"]
    if ctx.gpu_id:
        cmd.extend(["-gpu_id", ctx.gpu_id])
    rc = _run(cmd, cwd=em_dir, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("em", "failed", detail=f"EM mdrun rc={rc}")
    return StageResult("em", "done",
                       artifacts=[str(em_dir / "EM.gro")])


def stage_replica_md(ctx: CaseContext) -> StageResult:
    """Per (direction × replica) NVT → NPT → production MD.

    Handles both short-chain (single 3-replica chain) and long-chain
    (3 directions × 3 replicas = 9 sub-runtime) cases.
    """
    log = ctx.case_dir / "logs" / "replica_md.log"
    npt_mdp = MDP_BUILD_DIR / "NPT.mdp"
    md_mdp = MDP_BUILD_DIR / "MD.mdp"
    n_ok = 0
    for direction, rep in _iter_sub_runtimes(ctx):
        rdir = _replica_runtime_dir(ctx, direction, rep)
        topol_path = _build_dir(ctx, direction) / "outputs" / "topol.top"
        index_path = _build_dir(ctx, direction) / "outputs" / "index.ndx"
        label = f"r{rep}" if not direction else f"{direction}/r{rep}"
        rc, info = _run_grompp_mdrun_chain(
            ctx, rdir,
            label=label,
            npt_mdp=npt_mdp, md_mdp=md_mdp,
            topol_path=topol_path,
            index_path=index_path,
            log_path=log,
        )
        if rc != 0:
            return StageResult("replica_md", "failed",
                               detail=f"{info} rc={rc}")
        n_ok += 1
    return StageResult("replica_md", "done",
                       detail=f"{n_ok} sub-runtime(s)")


REQUIRED_POSTPROCESS_XVG: tuple[str, ...] = (
    "rmsd_backbone.xvg",
    "rmsf_residue.xvg",
    "hbonds.xvg",
    "mindist.xvg",
    "contacts.xvg",
    "sasa_total.xvg",
)


def _postprocess_complete(xvg_dir: Path) -> bool:
    return all((xvg_dir / name).exists() and (xvg_dir / name).stat().st_size > 0
               for name in REQUIRED_POSTPROCESS_XVG)


def _md_finished_for_postprocess(md_dir: Path) -> tuple[bool, str]:
    required = ("MD.tpr", "MD.xtc", "MD.log")
    missing = [name for name in required
               if not (md_dir / name).exists() or (md_dir / name).stat().st_size <= 0]
    if missing:
        return False, f"missing/empty {','.join(missing)}"
    try:
        with (md_dir / "MD.log").open("r", encoding="utf-8",
                                      errors="replace") as handle:
            if any("Finished mdrun" in line for line in handle):
                return True, "finished"
    except OSError as exc:
        return False, f"cannot read MD.log: {exc}"
    return False, "MD.log lacks Finished mdrun"


def _run_postprocess_subruntime(ctx: CaseContext, direction: str | None,
                                rep: int, log: Path) -> tuple[int, list[str]]:
    script = SCRIPTS_DIR / "run_post_pipeline.py"
    rdir = _replica_runtime_dir(ctx, direction, rep)
    md_dir = _md_output_dir(rdir)
    out_dir = _postprocess_runtime_dir(ctx, direction, rep)
    out_dir.mkdir(parents=True, exist_ok=True)
    rc = _run([sys.executable, str(script),
               "--tasks", "rmsd,rmsf,hbond,mindist,sasa",
               "--max-workers", "5",
               "--tpr", str(md_dir / "MD.tpr"),
               "--xtc", str(md_dir / "MD.xtc"),
               "--index", str(_build_dir(ctx, direction) / "outputs" / "index.ndx"),
               "--gmx", str(ctx.gmx_bin),
               "--outdir", str(out_dir)],
              dry_run=ctx.dry_run, log_path=log)
    artifacts = [str(out_dir / "xvg" / name)
                 for name in REQUIRED_POSTPROCESS_XVG]
    return rc, artifacts


def stage_postprocess_completed(ctx: CaseContext) -> StageResult:
    """Cache postprocess XVGs for replicas whose production MD is complete.

    This stage is safe to run while other replicas are still active. It only
    processes sub-runtimes with complete MD inputs and ``Finished mdrun`` in
    ``MD.log``. It does not run gate, aggregate, fingerprint, visualization,
    or final audits; those remain full-case stages.
    """
    log = ctx.case_dir / "logs" / "postprocess_completed.log"
    manifest = {
        "schema_version": "autopropet_v1_1_postprocess_completed_cache",
        "case_dir": str(ctx.case_dir),
        "processed": [],
        "skipped_already_done": [],
        "skipped_unfinished": [],
        "failed": [],
        "artifacts": [],
    }
    for direction, rep in _iter_sub_runtimes(ctx):
        label = f"{direction}/r{rep}" if direction else f"r{rep}"
        rdir = _replica_runtime_dir(ctx, direction, rep)
        md_dir = _md_output_dir(rdir)
        out_dir = _postprocess_runtime_dir(ctx, direction, rep)
        xvg_dir = out_dir / "xvg"
        if not ctx.dry_run and _postprocess_complete(xvg_dir):
            manifest["skipped_already_done"].append(label)
            manifest["artifacts"].extend(str(xvg_dir / name)
                                         for name in REQUIRED_POSTPROCESS_XVG)
            continue
        if not ctx.dry_run:
            ready, reason = _md_finished_for_postprocess(md_dir)
            if not ready:
                manifest["skipped_unfinished"].append({
                    "subruntime": label,
                    "reason": reason,
                    "md_dir": str(md_dir),
                })
                continue
        rc, artifacts = _run_postprocess_subruntime(ctx, direction, rep, log)
        if rc != 0 and not ctx.dry_run:
            manifest["failed"].append({
                "subruntime": label,
                "returncode": rc,
            })
            continue
        manifest["processed"].append(label)
        manifest["artifacts"].extend(artifacts)

    if not ctx.dry_run:
        qc_dir = ctx.case_dir / "qc"
        qc_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = qc_dir / "postprocess_cache_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2,
                                            sort_keys=True) + "\n",
                                 encoding="utf-8")
        md_path = qc_dir / "postprocess_cache_manifest.md"
        md_path.write_text(
            "\n".join([
                "# AutoProPET v1.1 postprocess completed cache",
                "",
                f"Processed: {len(manifest['processed'])}",
                f"Skipped already done: {len(manifest['skipped_already_done'])}",
                f"Skipped unfinished: {len(manifest['skipped_unfinished'])}",
                f"Failed: {len(manifest['failed'])}",
                "",
                f"JSON: `{manifest_path}`",
            ]) + "\n",
            encoding="utf-8",
        )

    if manifest["failed"]:
        return StageResult("postprocess_completed", "failed",
                           detail=f"{len(manifest['failed'])} failure(s)",
                           artifacts=manifest["artifacts"])
    return StageResult(
        "postprocess_completed",
        "done",
        detail=(
            f"processed={len(manifest['processed'])}; "
            f"skipped_unfinished={len(manifest['skipped_unfinished'])}; "
            f"skipped_done={len(manifest['skipped_already_done'])}"
        ),
        artifacts=manifest["artifacts"],
    )


def stage_postprocess(ctx: CaseContext) -> StageResult:
    """Trajectory post-processing per (direction × replica) sub-runtime.

    Produces the complete case-local xvg/png bundle used by the deeper
    figure scripts: RMSD, RMSF, H-bonds, min-distance/contact counts, and
    SASA.
    """
    log = ctx.case_dir / "logs" / "postprocess.log"
    n_ok = 0
    artifacts: list[str] = []
    for direction, rep in _iter_sub_runtimes(ctx):
        rdir = _replica_runtime_dir(ctx, direction, rep)
        md_dir = _md_output_dir(rdir)
        out_dir = _postprocess_runtime_dir(ctx, direction, rep)
        xvg_dir = out_dir / "xvg"
        md_xtc = md_dir / "MD.xtc"
        md_tpr = md_dir / "MD.tpr"
        if not md_xtc.exists() and not ctx.dry_run:
            return StageResult("postprocess", "failed",
                               detail=f"MD.xtc missing in {md_dir}")
        if not ctx.dry_run and _postprocess_complete(xvg_dir):
            n_ok += 1
            artifacts.extend(str(xvg_dir / name)
                             for name in REQUIRED_POSTPROCESS_XVG)
            continue
        rc, sub_artifacts = _run_postprocess_subruntime(ctx, direction, rep, log)
        if rc != 0 and not ctx.dry_run:
            label = f"{direction}/r{rep}" if direction else f"r{rep}"
            return StageResult("postprocess", "failed",
                               detail=f"postprocess {label} rc={rc}")
        n_ok += 1
        artifacts.extend(sub_artifacts)
    return StageResult("postprocess", "done",
                       detail=f"{n_ok} sub-runtime(s)",
                       artifacts=artifacts)


def stage_gate(ctx: CaseContext) -> StageResult:
    """Per (direction × replica) 4-layer production gate."""
    script = SCRIPTS_DIR / "check_v1_production_gate.py"
    log = ctx.case_dir / "logs" / "gate.log"
    paths: list[str] = []
    for direction, rep in _iter_sub_runtimes(ctx):
        rdir = _replica_runtime_dir(ctx, direction, rep)
        out_name = _verdict_filename(direction, rep)
        out_path = ctx.case_dir / "qc" / out_name
        cmd = [sys.executable, str(script),
               "--case-dir", str(rdir),
               "--replica-id", str(rep),
               "--pet-kind", ctx.pet_kind,
               "--burn-in-ns", str(production_burn_in_ns(ctx)),
               "--gmx-bin", str(ctx.gmx_bin),
               "--out", str(out_path)]
        if direction:
            cmd.extend(["--direction", direction])
        rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
        if rc not in (0, 1):
            label = f"{direction}/r{rep}" if direction else f"r{rep}"
            return StageResult("gate", "failed",
                               detail=f"gate {label} unexpected rc={rc}")
        paths.append(str(out_path))
    return StageResult("gate", "done", artifacts=paths)


def stage_aggregate_gate(ctx: CaseContext) -> StageResult:
    """Combine every per (direction × replica) verdict into one case-level
    verdict_aggregated.json."""
    script = SCRIPTS_DIR / "check_v1_production_gate.py"
    log = ctx.case_dir / "logs" / "aggregate_gate.log"
    verdict_files = [str(ctx.case_dir / "qc" / _verdict_filename(d, r))
                     for d, r in _iter_sub_runtimes(ctx)]
    expected = ctx.replica_count * (len(ctx.directions) if ctx.is_long_chain else 1)
    out = ctx.case_dir / "qc" / "v1_production_gate_verdict_aggregated.json"
    cmd = [sys.executable, str(script),
           "--aggregate",
           "--pet-kind", ctx.pet_kind,
           "--expected-replica-count", str(expected),
           "--verdict-files", *verdict_files,
           "--out", str(out)]
    rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
    if rc not in (0, 1):
        return StageResult("aggregate_gate", "failed",
                           detail=f"unexpected rc={rc}")
    return StageResult("aggregate_gate", "done", artifacts=[str(out)])


def stage_fingerprint(ctx: CaseContext) -> StageResult:
    """Per (direction × replica) v1.1 2-layer protocol fingerprint."""
    script = SCRIPTS_DIR / "compute_protocol_fingerprint.py"
    protein_id = (ctx.case_data.get("protein") or {}).get("id", "00000")
    protein_pdb = (ctx.case_data.get("protein") or {}).get("asset")
    if protein_pdb:
        protein_pdb = SIM_ROOT / protein_pdb / "source" / f"{protein_id}.pdb"
    log = ctx.case_dir / "logs" / "fingerprint.log"
    artifacts: list[str] = []
    for direction, rep in _iter_sub_runtimes(ctx):
        rdir = _replica_runtime_dir(ctx, direction, rep)
        nvt_mdp = rdir / "NVT.mdp"
        out = ctx.case_dir / "qc" / _fingerprint_filename(direction, rep)
        selection_meta = ctx.case_dir / "qc" / (
            "v1_candidate_selection_per_direction.json"
            if ctx.is_long_chain
            else "v1_candidate_selection.json"
        )
        cmd = [sys.executable, str(script),
               "--case-dir", str(ctx.case_dir),
               "--case-yaml", str(ctx.case_yaml),
               "--pet-kind", ctx.pet_kind,
               "--protein-pdb", str(protein_pdb) if protein_pdb else "/dev/null",
               "--replica-idx", str(rep),
               "--box-dim-nm", str(ctx.box_dim_nm),
               "--nvt-mdp", str(nvt_mdp),
               "--selection-meta", str(selection_meta),
               "--build-provenance", str(_build_dir(ctx, direction) / "outputs" / "build_provenance.json"),
               "--gmx-binary", str(ctx.gmx_bin),
               "--output", str(out)]
        vina_binary = _vina_env_tool("vina")
        if vina_binary is not None:
            cmd.extend(["--vina-binary", str(vina_binary)])
        if direction:
            cmd.extend(["--direction", direction])
        rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            label = f"{direction}/r{rep}" if direction else f"r{rep}"
            return StageResult("fingerprint", "failed",
                               detail=f"fingerprint {label} rc={rc}")
        artifacts.append(str(out))
    return StageResult("fingerprint", "done", artifacts=artifacts)


# ---------------------------------------------------------------------------
# Long-chain stages (R1 Phase B + chain_relax)
# ---------------------------------------------------------------------------

def stage_screen_r1_phase_b(ctx: CaseContext) -> StageResult:
    """Multi-conformer Kabsch screen for L10/L20.

    Invokes ``screen_vina_pet_multi_conformer.py`` which loops over the K=5
    pre-equilibrated PET conformers from
    ``inputs/pets/<kind>/preequilibrated/`` and feeds each into the
    existing single-conformer screen. Emits
    ``screening_summary_multi_conformer.json`` covering ~70 (L10) or ~170
    (L20) candidate poses tagged with conformer_idx + chain_extension_direction
    + post_kabsch_envelope + ranking_score.
    """
    if not ctx.is_long_chain:
        return StageResult("screen_r1_phase_b", "skipped",
                           detail="short-chain uses single-conformer align stage")
    script = SCRIPTS_DIR / "screen_vina_pet_multi_conformer.py"
    conformer_dir = SIM_ROOT / f"inputs/pets/{ctx.pet_kind}/preequilibrated"
    out_dir = ctx.case_dir / "03_pet_alignment_screen"
    log = ctx.case_dir / "logs" / "screen_r1_phase_b.log"
    target_unit_map, map_error = ensure_pet_unit_map(ctx, log_path=log)
    if map_error:
        return StageResult("screen_r1_phase_b", "failed", detail=map_error)

    docking_cfg = ctx.case_data.get("docking") or {}
    active_site_registry = SIM_ROOT / "inputs/active_site_registry.yaml"
    protein_cfg = ctx.case_data.get("protein") or {}
    protein_id = str(protein_cfg.get("id") or "00000").zfill(5)
    registry_entry: dict = {}
    if active_site_registry.exists():
        try:
            import yaml  # type: ignore
            registry = yaml.safe_load(active_site_registry.read_text(encoding="utf-8")) or {}
            registry_entry = (
                (registry.get("proteins") or {}).get(protein_id)
                or registry.get(protein_id)
                or {}
            )
        except Exception:
            registry_entry = {}
    case_triad = protein_cfg.get("active_site_triad") or {}
    registry_triad = registry_entry.get("catalytic_triad") or {}
    ser = case_triad.get("ser") or docking_cfg.get("ser_resid") or registry_triad.get("ser") or 131
    his = case_triad.get("his") or docking_cfg.get("his_resid") or registry_triad.get("his") or 209
    acid = case_triad.get("acid") or docking_cfg.get("acid_resid") or registry_triad.get("acid") or 177
    ser_atom = (
        protein_cfg.get("ser_atom")
        or docking_cfg.get("ser_atom_name")
        or registry_entry.get("ser_atom")
        or "OG"
    )
    protein_pdb = protein_cfg.get("asset")
    if protein_pdb:
        protein_pdb = SIM_ROOT / protein_pdb / "source" / f"{protein_id}.pdb"

    cmd = [sys.executable, str(script),
           "--conformer-dir", str(conformer_dir),
           "--target-pet-kind", ctx.pet_kind,
           "--target-unit-map", str(target_unit_map),
           "--docking-dir", str(ctx.case_dir / "02_docking"),
           "--protein-pdb", str(protein_pdb) if protein_pdb else "/dev/null",
           "--active-site-registry", str(active_site_registry),
           "--ser-resid", str(ser),
           "--ser-atom-name", str(ser_atom),
           "--his-resid", str(his),
           "--acid-resid", str(acid),
           "--out-dir", str(out_dir),
           "--gmx", str(ctx.gmx_bin)]
    rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("screen_r1_phase_b", "failed",
                           detail=f"screen rc={rc}; see {log}")
    return StageResult("screen_r1_phase_b", "done",
                       artifacts=[str(out_dir /
                                       "screening_summary_multi_conformer.json")])


def stage_select_per_direction(ctx: CaseContext) -> StageResult:
    """Pick the lowest-ranking-score clash-free pose per chain-extension
    direction (long-chain only)."""
    if not ctx.is_long_chain:
        return StageResult("select_per_direction", "skipped",
                           detail="short-chain uses single-pick select stage")
    script = SCRIPTS_DIR / "select_v1_binding_candidate.py"
    screen_dir = ctx.case_dir / "03_pet_alignment_screen"
    target = ctx.case_dir / "02_docking" / "best_rigidbody.pdb"
    meta = ctx.case_dir / "qc" / "v1_candidate_selection_per_direction.json"
    log = ctx.case_dir / "logs" / "select_per_direction.log"
    cmd = [sys.executable, str(script),
           "--screen-dir", str(screen_dir),
           "--target-pdb", str(target),
           "--meta", str(meta),
           "--per-direction"]
    rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
    if rc != 0 and not ctx.dry_run:
        return StageResult("select_per_direction", "failed",
                           detail=f"rc={rc}; see {log}")
    return StageResult("select_per_direction", "done",
                       artifacts=[str(target), str(meta)])


def stage_build_per_direction(ctx: CaseContext) -> StageResult:
    """Build CHARMM topology + cubic box for each chain-extension direction.

    Long-chain only — each direction gets its own
    ``03_md_build/dir_<direction>/`` with the appropriate selected pose
    handed in by the previous selector stage.
    """
    if not ctx.is_long_chain:
        return StageResult("build_per_direction", "skipped",
                           detail="short-chain uses single build stage")
    script = SCRIPTS_DIR / "build_md_from_docking_best.py"
    log = ctx.case_dir / "logs" / "build_per_direction.log"
    target_stem = ctx.case_dir / "02_docking" / "best_rigidbody.pdb"
    protein_id = str((ctx.case_data.get("protein") or {}).get("id") or "00000")
    artifacts: list[str] = []
    for direction in ctx.directions:
        per_dir_input = target_stem.with_name(
            f"best_rigidbody_{direction}.pdb")
        per_dir_build_dir = _build_dir(ctx, direction)
        per_dir_build_dir.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(script),
               "--best-model", str(per_dir_input),
               "--build-dir", str(per_dir_build_dir),
               "--pet-kind", ctx.pet_kind,
               "--box-dim-nm", str(ctx.box_dim_nm),
               "--protein-id", protein_id,
               "--gmx-bin", str(ctx.gmx_bin),
               "--genion-seed", str(_build_genion_seed(ctx, direction))]
        rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("build_per_direction", "failed",
                               detail=f"build {direction} rc={rc}")
        artifacts.append(str(per_dir_build_dir))
    return StageResult("build_per_direction", "done", artifacts=artifacts)


def stage_chain_relax(ctx: CaseContext) -> StageResult:
    """Phase C chain_relax stage (L10/L20 only).

    Per direction: generate ``posre_pet_anchor.itp`` for the selected
    window, then EM + NVT 200 ps + NPT 200 ps with POSRES on the protein
    backbone (K=1000) and the L4 anchor (K=500). Output feeds production
    MD. This stage uses the chain-relaxation MDP templates.
    """
    if not ctx.is_long_chain or not ctx.chain_relax_enabled:
        return StageResult("chain_relax", "skipped",
                           detail="not a long-chain case or chain_relax disabled")
    posres_script = SCRIPTS_DIR / "generate_pet_anchor_posres.py"
    em_mdp = MDP_CHAIN_RELAX_DIR / "em_chain_relax.mdp"
    nvt_mdp = MDP_CHAIN_RELAX_DIR / "nvt_chain_relax.mdp"
    npt_mdp = MDP_CHAIN_RELAX_DIR / "npt_chain_relax.mdp"
    log = ctx.case_dir / "logs" / "chain_relax.log"

    # Selection meta from the per-direction selector tells us which window
    # to anchor on per direction.
    sel_meta_path = ctx.case_dir / "qc" / "v1_candidate_selection_per_direction.json"
    sel_meta: dict = {}
    if sel_meta_path.exists() and not ctx.dry_run:
        try:
            sel_meta = json.loads(sel_meta_path.read_text(encoding="utf-8")) \
                       .get("per_direction") or {}
        except json.JSONDecodeError:
            sel_meta = {}

    unit_map, map_error = ensure_pet_unit_map(ctx, log_path=log)
    if map_error:
        return StageResult("chain_relax", "failed", detail=map_error)
    p1_itp = SIM_ROOT / (f"inputs/pets/{ctx.pet_kind}/raw/charmm-gui-001/"
                          f"gromacs/toppar/P1.itp")

    n_ok = 0
    for direction in ctx.directions:
        dirname = _direction_dirname(direction)
        crelax_dir = (ctx.case_dir / "04_chain_relax" / dirname).resolve()
        crelax_dir.mkdir(parents=True, exist_ok=True)
        for sub in ("EM", "NVT", "NPT"):
            (crelax_dir / sub).mkdir(parents=True, exist_ok=True)

        # 1) generate posre_pet_anchor.itp per direction
        sel = sel_meta.get(direction) or {}
        ws = sel.get("window_start"); we = sel.get("window_end")
        kabsch_dir = sel.get("kabsch_direction") or "forward"
        if ws is None or we is None:
            # Default to head-end anchor when no selection meta is available
            # (allows dry-run + sanity checks to proceed).
            ws, we = 1, 4
        posre_out = (_build_dir(ctx, direction) / "outputs"
                     / "posre_pet_anchor.itp").resolve()
        posre_out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(posres_script),
               "--unit-map", str(unit_map),
               "--p1-itp", str(p1_itp),
               "--window", str(ws), str(we),
               "--direction", kabsch_dir,
               "--output", str(posre_out),
               "--case-id", f"{ctx.case_id}__{direction}"]
        rc = _run(cmd, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("chain_relax", "failed",
                               detail=f"posres gen {direction} rc={rc}")

        # 2) chain_relax EM + NVT + NPT (uses mdp_chain_relax/* with POSRES)
        build_outputs = (_build_dir(ctx, direction) / "outputs").resolve()
        complex_gro = build_outputs / "box_sol_ion.gro"
        topol = build_outputs / "topol.top"
        index_path = build_outputs / "index.ndx"
        if not ctx.dry_run:
            for required in (complex_gro, topol, index_path):
                if not required.exists():
                    return StageResult("chain_relax", "failed",
                                       detail=f"missing {required}")
            ensure_pet_anchor_posres_include(topol)
        em_dir = crelax_dir / "EM"
        nvt_dir = crelax_dir / "NVT"
        npt_dir = crelax_dir / "NPT"

        # EM
        rc = _run([str(ctx.gmx_bin), "grompp",
                   "-f", str(em_mdp),
                   "-c", str(complex_gro),
                   "-p", str(topol),
                   "-o", str(em_dir / "EM.tpr"), "-maxwarn", "20"],
                  cwd=em_dir, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("chain_relax", "failed",
                               detail=f"EM grompp {direction} rc={rc}")
        cmd = [str(ctx.gmx_bin), "mdrun", "-deffnm", "EM"]
        if ctx.gpu_id:
            cmd.extend(["-gpu_id", ctx.gpu_id])
        rc = _run(cmd, cwd=em_dir, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("chain_relax", "failed",
                               detail=f"EM mdrun {direction} rc={rc}")

        # NVT (POSRES protein + L4 anchor)
        rc = _run([str(ctx.gmx_bin), "grompp",
                   "-f", str(nvt_mdp),
                   "-c", str(em_dir / "EM.gro"),
                   "-r", str(em_dir / "EM.gro"),
                   "-p", str(topol),
                   "-n", str(index_path),
                   "-o", str(nvt_dir / "NVT.tpr"), "-maxwarn", "20"],
                  cwd=nvt_dir, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("chain_relax", "failed",
                               detail=f"NVT grompp {direction} rc={rc}")
        cmd = [str(ctx.gmx_bin), "mdrun", "-deffnm", "NVT"]
        if ctx.gpu_id:
            cmd.extend(["-gpu_id", ctx.gpu_id])
        rc = _run(cmd, cwd=nvt_dir, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("chain_relax", "failed",
                               detail=f"NVT mdrun {direction} rc={rc}")

        # NPT (POSRES continues)
        rc = _run([str(ctx.gmx_bin), "grompp",
                   "-f", str(npt_mdp),
                   "-c", str(nvt_dir / "NVT.gro"),
                   "-r", str(nvt_dir / "NVT.gro"),
                   "-t", str(nvt_dir / "NVT.cpt"),
                   "-p", str(topol),
                   "-n", str(index_path),
                   "-o", str(npt_dir / "NPT.tpr"), "-maxwarn", "20"],
                  cwd=npt_dir, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("chain_relax", "failed",
                               detail=f"NPT grompp {direction} rc={rc}")
        cmd = [str(ctx.gmx_bin), "mdrun", "-deffnm", "NPT"]
        if ctx.gpu_id:
            cmd.extend(["-gpu_id", ctx.gpu_id])
        rc = _run(cmd, cwd=npt_dir, dry_run=ctx.dry_run, log_path=log)
        if rc != 0 and not ctx.dry_run:
            return StageResult("chain_relax", "failed",
                               detail=f"NPT mdrun {direction} rc={rc}")

        # Hand off relaxed coords + checkpoint to the production EM slot
        # so stage_replica_md sees the per-direction shared EM.gro it needs.
        shared_em = _shared_em_gro(ctx, direction)
        shared_em.parent.mkdir(parents=True, exist_ok=True)
        if not ctx.dry_run and (npt_dir / "NPT.gro").exists():
            shared_em.write_bytes((npt_dir / "NPT.gro").read_bytes())
        n_ok += 1
    return StageResult("chain_relax", "done",
                       detail=f"{n_ok} direction(s)")


STAGE_REGISTRY: dict[str, Callable[[CaseContext], StageResult]] = {
    "materialize":           stage_materialize,
    "docking":               stage_docking,
    "align":                 stage_align,
    "select":                stage_select,
    "build":                 stage_build,
    "validate_box":          stage_validate_box,
    "em":                    stage_em,
    "replica_md":            stage_replica_md,
    "postprocess_completed":  stage_postprocess_completed,
    "postprocess":           stage_postprocess,
    "gate":                  stage_gate,
    "aggregate_gate":        stage_aggregate_gate,
    "fingerprint":           stage_fingerprint,
    "screen_r1_phase_b":     stage_screen_r1_phase_b,
    "select_per_direction":  stage_select_per_direction,
    "build_per_direction":   stage_build_per_direction,
    "chain_relax":           stage_chain_relax,
}


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def build_context(*, case_yaml: Path, runtime_root: Path,
                  gmx_bin: Path, gpu_id: str, dry_run: bool) -> CaseContext:
    mat = _load_mat()
    import yaml  # type: ignore
    case_data = yaml.safe_load(case_yaml.read_text(encoding="utf-8"))
    pet_kind = mat.pet_kind(case_data)
    case_id = mat.case_id(case_data)
    return CaseContext(
        case_yaml=case_yaml,
        case_data=case_data,
        case_id=case_id,
        pet_kind=pet_kind,
        is_long_chain=pet_kind in mat.LONG_CHAINS,
        case_dir=runtime_root / case_id,
        replica_count=mat.replica_count(case_data),
        directions=mat.directions(case_data, pet_kind),
        box_dim_nm=mat.box_dim_nm(case_data, pet_kind),
        chain_relax_enabled=mat.chain_relax_enabled(case_data, pet_kind),
        gmx_bin=gmx_bin,
        gpu_id=gpu_id,
        dry_run=dry_run,
    )


def select_pipeline(ctx: CaseContext) -> tuple[str, ...]:
    return PIPELINE_STAGES_LONG_CHAIN if ctx.is_long_chain \
           else PIPELINE_STAGES_SHORT_CHAIN


def run_pipeline(ctx: CaseContext, *,
                 stages: tuple[str, ...] | None = None,
                 from_stage: str | None = None,
                 to_stage: str | None = None,
                 force: bool = False) -> dict:
    """Resume-aware pipeline driver.

    Skips stages already marked ``done`` in
    ``qc/orchestration_summary.json`` unless ``force=True`` or
    ``from_stage`` is provided (in which case all stages from that point
    are re-run, regardless of prior status).
    """
    stages = stages or select_pipeline(ctx)
    if from_stage and from_stage not in stages:
        raise ValueError(f"--from-stage {from_stage!r} not in pipeline {stages}")
    if to_stage and to_stage not in stages:
        raise ValueError(f"--to-stage {to_stage!r} not in pipeline {stages}")

    if not ctx.dry_run:
        (ctx.case_dir / "logs").mkdir(parents=True, exist_ok=True)
    summary = _load_summary(ctx.case_dir)
    summary.setdefault("stages", {})
    summary.setdefault("case_id", ctx.case_id)
    summary.setdefault("case_yaml", str(ctx.case_yaml))

    started_from_idx = stages.index(from_stage) if from_stage else 0
    stop_after_idx = stages.index(to_stage) if to_stage else len(stages) - 1
    if started_from_idx > stop_after_idx:
        raise ValueError(
            f"--from-stage {from_stage!r} comes after --to-stage {to_stage!r}"
        )

    for i, name in enumerate(stages):
        if i < started_from_idx:
            print(f"[skip-before-from] {name}")
            continue
        prior = summary["stages"].get(name, {})
        if (not force and from_stage is None
                and prior.get("status") == "done"):
            print(f"[resume-skip] {name}")
            if i == stop_after_idx:
                break
            continue
        runner = STAGE_REGISTRY.get(name)
        if runner is None:
            summary["stages"][name] = {"status": "skipped",
                                        "detail": "no runner"}
            if not ctx.dry_run:
                _save_summary(ctx.case_dir, summary)
            if i == stop_after_idx:
                break
            continue
        try:
            result = runner(ctx)
        except NotImplementedError as exc:
            result = StageResult(name, "skipped", detail=str(exc))
        summary["stages"][name] = dataclasses.asdict(result)
        if not ctx.dry_run:
            _save_summary(ctx.case_dir, summary)
        if result.status == "failed":
            print(f"[fail] stage={name} detail={result.detail}", file=sys.stderr)
            return summary
        print(f"[done] stage={name} detail={result.detail}")
        if i == stop_after_idx:
            break
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("case_yaml", type=Path)
    p.add_argument("--runtime-root", type=Path,
                   default=SIM_ROOT / "runtime/adhoc")
    p.add_argument("--from-stage", default=None,
                   help="Re-run starting at this stage (skips earlier).")
    p.add_argument("--to-stage", default=None,
                   help="Stop after this stage; useful for staged sanity runs.")
    p.add_argument("--force", action="store_true",
                   help="Re-run every stage regardless of prior 'done' state.")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--gmx", type=Path,
                   default=Path(os.environ.get("GMX_BIN",
                                               "/usr/local/gromacs/bin/gmx")))
    p.add_argument("--gpu-id", type=str,
                   default=os.environ.get("GPU_ID", ""))
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.case_yaml.exists():
        print(f"case.yaml missing: {args.case_yaml}", file=sys.stderr)
        return 2
    ctx = build_context(
        case_yaml=args.case_yaml,
        runtime_root=args.runtime_root,
        gmx_bin=args.gmx,
        gpu_id=args.gpu_id,
        dry_run=args.dry_run,
    )
    summary = run_pipeline(
        ctx,
        from_stage=args.from_stage,
        to_stage=args.to_stage,
        force=args.force,
    )
    failed = any(v.get("status") == "failed"
                 for v in summary.get("stages", {}).values())
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
