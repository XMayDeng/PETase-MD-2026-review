"""v1.1 case-materialisation helpers.

Helpers for the R1 replica fan-out, direction layout, and case.yaml v1.1 schema.

What this module covers
-----------------------

* **case.yaml v1.1 schema parsers** — read ``r1_phase_b``,
  ``box``, ``chain_relax``, ``classification``, ``md.replicas``,
  ``ligand`` (with ``pet:`` legacy alias) sections; resolve ``auto`` values
  per chain length.
* **Hash-derived NVT gen-seed** — deterministic seed algorithm,
  masked to a signed 32-bit positive integer (GROMACS requirement;
  same convention as ``run_pet_solo_preeq.replica_seed``).
* **NVT.mdp seed rewriter** — overwrite the ``gen-seed`` line in an
  existing mdp template.
* **Sub-runtime directory layout** — create the per-direction × per-replica
  nested tree (L10/L20 → 9 sub-runtimes,
  L4/L2 → 3 sub-runtimes).
* **R1 Phase A asset presence check** — fail early when a long-chain
  case is materialised but ``inputs/pets/<kind>/preequilibrated/`` has
  no conformer files.
"""
from __future__ import annotations

import hashlib
import importlib.util
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SIM_ROOT = REPO_ROOT / "."

DEFAULT_REPLICA_COUNT = 3
LONG_CHAINS: frozenset[str] = frozenset({"PET_L10", "PET_L20"})
LONG_CHAIN_DIRECTIONS: tuple[str, str, str] = (
    "head-side", "tail-side", "bidirectional",
)


def _load_box_protocol():
    path = SIM_ROOT / "code/pipelines/box_protocol.py"
    spec = importlib.util.spec_from_file_location("box_protocol", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# case.yaml v1.1 section accessors
# ---------------------------------------------------------------------------

def _section(case: dict, key: str) -> dict:
    v = case.get(key)
    return v if isinstance(v, dict) else {}


def case_id(case: dict) -> str:
    """Return the canonical case identifier; raises if missing."""
    cid = case.get("id")
    if not cid:
        raise ValueError("case yaml has no top-level 'id' field")
    return str(cid)


def ligand_section(case: dict) -> dict:
    """v1.1 prefers ``ligand:``; falls back to legacy ``pet:`` alias."""
    primary = case.get("ligand")
    if isinstance(primary, dict) and primary:
        return primary
    return _section(case, "pet")


def pet_kind(case: dict) -> str:
    """Resolve the PET chain length (``PET_L4`` etc.) from ligand.kind
    (preferred) or legacy fallbacks."""
    lig = ligand_section(case)
    kind = lig.get("kind") or case.get("ligand_kind") or case.get("pet_kind")
    if not kind:
        raise ValueError("case yaml does not specify a PET chain length")
    return _normalise_pet_kind(str(kind))


def _normalise_pet_kind(value: str) -> str:
    upper = value.upper()
    if upper.startswith("PET_L"):
        return upper
    if upper.startswith("PET") and upper[3:].isdigit():
        return f"PET_L{upper[3:]}"
    return upper


def case_schema_version(case: dict) -> str:
    return str(case.get("case_schema_version") or case.get("schema_version") or "v1.0")


def replica_count(case: dict) -> int:
    """``md.replicas`` (v1.1) or top-level ``md_repeats`` fall-back; default 3."""
    md = _section(case, "md")
    if "replicas" in md:
        return int(md["replicas"])
    if "md_repeats" in case:
        return int(case["md_repeats"])
    return DEFAULT_REPLICA_COUNT


def directions(case: dict, pet_kind_value: str | None = None) -> list[str]:
    """Resolve the chain-extension directions for this case.

    L10 / L20 default to the three canonical directions
    ({head-side, tail-side, bidirectional}); L4 / L2 default to ``[]``
    (single-pose path). ``r1_phase_b.directions`` may override.
    """
    k = pet_kind_value or pet_kind(case)
    r1 = _section(case, "r1_phase_b")
    explicit = r1.get("directions")
    if isinstance(explicit, list) and explicit:
        cleaned = [str(d) for d in explicit if d is not None]
        return cleaned
    if k in LONG_CHAINS:
        return list(LONG_CHAIN_DIRECTIONS)
    return []


def chain_relax_enabled(case: dict, pet_kind_value: str | None = None) -> bool:
    """``chain_relax.enabled`` honoured if set explicitly; ``auto``
    activates for L10/L20 only."""
    cr = _section(case, "chain_relax")
    flag = cr.get("enabled", "auto")
    if isinstance(flag, bool):
        return flag
    if isinstance(flag, str) and flag.lower() != "auto":
        return flag.lower() in ("true", "yes", "on", "1")
    k = pet_kind_value or pet_kind(case)
    return k in LONG_CHAINS


def box_dim_nm(case: dict, pet_kind_value: str | None = None) -> float:
    """Resolve the production cubic box edge (nm).

    ``box.box_dim_nm`` honoured if numeric; ``auto`` (default) defers to
    ``box_protocol.COMPLEX_BOX_DIM_NM`` for the case's chain length.
    """
    bx = _section(case, "box")
    raw = bx.get("box_dim_nm", "auto")
    if isinstance(raw, (int, float)):
        return float(raw)
    bp = _load_box_protocol()
    return bp.box_dim_for(pet_kind_value or pet_kind(case))


def is_negative_control(case: dict) -> bool:
    return bool(_section(case, "classification").get("is_negative_control", False))


def positive_protein_family(case: dict) -> str:
    return str(_section(case, "classification").get(
        "positive_protein_family") or "PETase")


# ---------------------------------------------------------------------------
# Hash-derived NVT gen-seed
# ---------------------------------------------------------------------------

def replica_seed(case_id_value: str, replica_idx: int) -> int:
    """Deterministic positive 31-bit gen-seed for production NVT.

    Spec section 4.6 algorithm: ``int(sha256(f"{case_id}:r{idx}")[:8], 16)``,
    masked to 31 bits so it fits in GROMACS's signed 32-bit ``gen-seed``
    field (genion-seed shares the same constraint; the run_pet_solo_preeq
    namespace differs only by prefix).
    """
    payload = f"{case_id_value}:r{int(replica_idx)}"
    return int(hashlib.sha256(payload.encode()).hexdigest()[:8], 16) & 0x7FFFFFFF


# ---------------------------------------------------------------------------
# NVT.mdp gen-seed rewriting
# ---------------------------------------------------------------------------

def rewrite_nvt_mdp_with_seed(mdp_text: str, gen_seed: int) -> str:
    """Return ``mdp_text`` with the ``gen-seed`` / ``gen_seed`` line
    replaced. Raises ``ValueError`` when no seed line is present."""
    out: list[str] = []
    replaced = False
    for raw in mdp_text.splitlines():
        stripped = raw.lstrip()
        if stripped.startswith(("gen-seed", "gen_seed")):
            head, _, comment = raw.partition(";")
            indent_len = len(raw) - len(stripped)
            indent = raw[:indent_len]
            comment_suffix = f"  ;{comment}" if comment else ""
            out.append(
                f"{indent}gen-seed                = {gen_seed}{comment_suffix}"
            )
            replaced = True
        else:
            out.append(raw)
    if not replaced:
        raise ValueError("NVT.mdp has no gen-seed line to rewrite")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# R1 Phase A asset presence check
# ---------------------------------------------------------------------------

def assert_conformer_asset_present(pet_kind_value: str, *,
                                   sim_root: Path | None = None) -> list[Path]:
    """For long-chain cases, ensure the R1 Phase A frozen conformer asset
    exists. Returns the sorted list of conformer .gro files."""
    if pet_kind_value not in LONG_CHAINS:
        return []
    root = sim_root or SIM_ROOT
    asset_dir = root / f"inputs/pets/{pet_kind_value}/preequilibrated"
    if not asset_dir.is_dir():
        raise FileNotFoundError(
            f"R1 Phase A conformer asset dir missing: {asset_dir}. "
            f"Run extract_pet_conformers.py before materialising "
            f"long-chain cases."
        )
    gros = sorted(asset_dir.glob("conformer_*.gro"))
    if not gros:
        raise FileNotFoundError(
            f"No conformer_*.gro files under {asset_dir}; run "
            f"extract_pet_conformers.py to build the asset."
        )
    return gros


# ---------------------------------------------------------------------------
# Sub-runtime directory layout
# ---------------------------------------------------------------------------

_DIR_SANITISER = re.compile(r"[^A-Za-z0-9_-]")


def _direction_dirname(direction: str) -> str:
    """Canonical sub-directory name for one chain-extension direction.

    The direction label (``head-side`` / ``tail-side`` / ``bidirectional``)
    is sanitised to keep the filesystem name shell-safe while preserving
    the dash separators used in the protocol fingerprint metadata.
    """
    return f"dir_{_DIR_SANITISER.sub('_', direction)}"


def _make_stage_tree(stage_dir: Path, leaves: list[str]) -> None:
    for leaf in leaves:
        (stage_dir / leaf).mkdir(parents=True, exist_ok=True)


def materialize_v1_1_subruntime_tree(
    *,
    case_dir: Path,
    pet_kind_value: str,
    directions_list: list[str],
    replica_count_value: int,
    chain_relax: bool,
) -> dict:
    """Create the v1.1 sub-runtime directory skeleton.

    Returns a manifest dict listing every created stage path, suitable for
    downstream stage runners + the fingerprint writer.

    L10 / L20 layout::

        case_dir/
          02_docking/                      (shared)
          03_pet_alignment_screen/         (shared)
          03_md_build/dir_<direction>/outputs/
          04_chain_relax/dir_<direction>/{EM,NVT,NPT}/
          05_md_prod/dir_<direction>/EM/
          05_md_prod/dir_<direction>/r<n>/{NVT,NPT,MD}/
          06_postprocess/dir_<direction>/r<n>/
          07_analysis/dir_<direction>/r<n>/
          08_visualization/case_summary/
          08_visualization/dir_<direction>/summary/
          qc/

    L4 / L2 layout — same minus the per-direction layer; chain_relax is
    skipped entirely; ``03_pet_alignment_screen`` is omitted.
    """
    if replica_count_value < 1:
        raise ValueError(f"replica_count must be ≥ 1, got {replica_count_value}")
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "02_docking").mkdir(parents=True, exist_ok=True)
    (case_dir / "qc").mkdir(parents=True, exist_ok=True)

    is_long = pet_kind_value in LONG_CHAINS
    manifest = {
        "case_dir": str(case_dir),
        "pet_kind": pet_kind_value,
        "is_long_chain": is_long,
        "replica_count": replica_count_value,
        "directions": list(directions_list) if is_long else [],
        "chain_relax_enabled": bool(chain_relax) and is_long,
        "sub_runtimes": [],
    }

    if is_long:
        (case_dir / "03_pet_alignment_screen").mkdir(parents=True, exist_ok=True)
        for direction in directions_list:
            d = _direction_dirname(direction)
            (case_dir / "03_md_build" / d / "outputs").mkdir(parents=True, exist_ok=True)
            if chain_relax:
                _make_stage_tree(case_dir / "04_chain_relax" / d,
                                 ["EM", "NVT", "NPT"])
            (case_dir / "05_md_prod" / d / "EM").mkdir(parents=True, exist_ok=True)
            (case_dir / "08_visualization" / d / "summary").mkdir(parents=True, exist_ok=True)
            for rep in range(1, replica_count_value + 1):
                rname = f"r{rep}"
                _make_stage_tree(case_dir / "05_md_prod" / d / rname,
                                 ["NVT", "NPT", "MD"])
                for stage in ("06_postprocess", "07_analysis"):
                    (case_dir / stage / d / rname).mkdir(parents=True, exist_ok=True)
                manifest["sub_runtimes"].append(
                    {"direction": direction, "replica_idx": rep,
                     "md_dir": str(case_dir / "05_md_prod" / d / rname)}
                )
        (case_dir / "08_visualization" / "case_summary").mkdir(parents=True, exist_ok=True)
    else:
        # Short-chain: single docked pose × N replicas, no Phase C, no screen.
        (case_dir / "03_md_build" / "outputs").mkdir(parents=True, exist_ok=True)
        (case_dir / "05_md_prod" / "EM").mkdir(parents=True, exist_ok=True)
        (case_dir / "08_visualization" / "case_summary").mkdir(parents=True, exist_ok=True)
        (case_dir / "08_visualization" / "summary").mkdir(parents=True, exist_ok=True)
        for rep in range(1, replica_count_value + 1):
            rname = f"r{rep}"
            _make_stage_tree(case_dir / "05_md_prod" / rname,
                             ["NVT", "NPT", "MD"])
            for stage in ("06_postprocess", "07_analysis"):
                (case_dir / stage / rname).mkdir(parents=True, exist_ok=True)
            manifest["sub_runtimes"].append(
                {"direction": None, "replica_idx": rep,
                 "md_dir": str(case_dir / "05_md_prod" / rname)}
            )
    return manifest


def materialize_v1_1_nvt_mdps(
    *,
    case_dir: Path,
    case_id_value: str,
    pet_kind_value: str,
    directions_list: list[str],
    replica_count_value: int,
    nvt_template_path: Path,
) -> list[Path]:
    """Write per-replica NVT.mdp into every sub-runtime, each with a
    hash-derived ``gen-seed``.

    Returns the list of output mdp paths in canonical (direction, replica)
    order.
    """
    template = nvt_template_path.read_text(encoding="utf-8")
    is_long = pet_kind_value in LONG_CHAINS
    out_paths: list[Path] = []
    if is_long:
        for direction in directions_list:
            d = _direction_dirname(direction)
            for rep in range(1, replica_count_value + 1):
                seed = replica_seed(f"{case_id_value}__{direction}", rep)
                target = case_dir / "05_md_prod" / d / f"r{rep}" / "NVT.mdp"
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(rewrite_nvt_mdp_with_seed(template, seed),
                                  encoding="utf-8")
                out_paths.append(target)
    else:
        for rep in range(1, replica_count_value + 1):
            seed = replica_seed(case_id_value, rep)
            target = case_dir / "05_md_prod" / f"r{rep}" / "NVT.mdp"
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(rewrite_nvt_mdp_with_seed(template, seed),
                              encoding="utf-8")
            out_paths.append(target)
    return out_paths
