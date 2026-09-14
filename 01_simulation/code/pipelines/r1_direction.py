"""R1 Phase B direction-mapping + envelope helpers (spec §4.1).

The long-chain Kabsch screen iterates over K=5 conformers × W windows ×
{forward, reverse}. The selector then picks the top-1 candidate per
**chain extension direction**: head-side / tail-side / bidirectional.

Chain-extension direction is *derived* from the L4 anchor window position
on the long-chain target: when the anchor occupies the head-most window,
all new monomers grow toward the tail (tail-side extension); the converse
defines head-side; a middle window yields bidirectional extension.

This module also computes the post-Kabsch bounding-box envelope used as
a compactness term in the screen ranking score. Smaller envelope ⇒ more
compact post-alignment placement ⇒ less likely to clash with the
periodic-image boundary in the chain-length-aware cubic box.

Both helpers are shared across:

* ``code/scripts/screen_vina_pet_multi_conformer.py`` — multi-conformer
  R1 Phase B orchestrator.
* ``code/scripts/select_v1_binding_candidate.py`` — per-direction selector
  (``--per-direction`` mode).
* ``code/pipelines/materialize_case.py`` (Task 5) — multi-direction
  build / chain_relax / production fan-out.
"""
from __future__ import annotations

from typing import Iterable, Sequence

# Chain-extension direction labels (3 fundamental sub-runtime branches
# per long-chain case, spec §4.1 "Chain extension 必跑 3 方向").
DIRECTIONS: tuple[str, str, str] = (
    "head-side",
    "tail-side",
    "bidirectional",
)

# Residue counts for supported long-chain PET targets.
PET_RESIDUE_COUNT: dict[str, int] = {
    "PET_L10": 10,
    "PET_L20": 20,
}

# Anchor window size = 4 PET residues (L4 motif).
ANCHOR_WINDOW_SIZE: int = 4


def infer_chain_extension_direction(
    window: tuple[int, int] | Sequence[int],
    target_pet_kind: str,
) -> str:
    """Return the chain-extension direction implied by an L4 anchor window.

    The window is given as a (start, end) residue pair, 1-indexed inclusive.

    Mapping rules (spec §4.1 "Chain extension 必跑 3 方向"):

        * window[0] == 1                 → tail-side  (anchor at head; tail grows)
        * window[1] == N_residues        → head-side  (anchor at tail; head grows)
        * otherwise (middle window)      → bidirectional

    Raises ``ValueError`` if the window is malformed or out-of-range.
    """
    if target_pet_kind not in PET_RESIDUE_COUNT:
        raise ValueError(
            f"Unknown target_pet_kind={target_pet_kind!r}; supported = "
            f"{sorted(PET_RESIDUE_COUNT)}"
        )
    start, end = int(window[0]), int(window[1])
    n_res = PET_RESIDUE_COUNT[target_pet_kind]
    if start < 1 or end > n_res or start > end:
        raise ValueError(
            f"Bad window {window!r} for {target_pet_kind} "
            f"(must be 1 <= start <= end <= {n_res})"
        )
    if end - start + 1 != ANCHOR_WINDOW_SIZE:
        raise ValueError(
            f"Window {window!r} size != {ANCHOR_WINDOW_SIZE} "
            f"(spec §4.1 anchors L4 motif)"
        )
    if start == 1:
        return "tail-side"
    if end == n_res:
        return "head-side"
    return "bidirectional"


def compute_post_kabsch_envelope_nm3(
    coords: Iterable[tuple[float, float, float]],
) -> float:
    """Bounding-box volume (nm³) of transformed PET coordinates.

    Lower is more compact. Used as the ``envelope`` term of the R1 Phase B
    ranking score (spec §4.1 "post_kabsch_box_envelope ranking term").
    """
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for x, y, z in coords:
        xs.append(x)
        ys.append(y)
        zs.append(z)
    if not xs:
        raise ValueError("compute_post_kabsch_envelope_nm3: empty coordinate set")
    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    dz = max(zs) - min(zs)
    return dx * dy * dz


def ranking_score(
    *,
    motif_rmsd_angstrom: float,
    contacts_lt_1a: int,
    affinity_kcal_mol: float,
    envelope_nm3: float,
    weights: dict[str, float] | None = None,
) -> float:
    """Composite ranking score; lower is better.

    Default weights (spec §5.1 ``kabsch_ranking_weights`` default block):

        motif_rmsd  : 0.4
        clash       : 0.3
        affinity    : 0.2
        envelope    : 0.1

    The Vina ``affinity_kcal_mol`` is most-negative-is-best, so it is
    used directly (more negative ⇒ lower contribution).
    """
    w = weights or {
        "motif_rmsd": 0.4,
        "clash": 0.3,
        "affinity": 0.2,
        "envelope": 0.1,
    }
    return (
        w["motif_rmsd"] * motif_rmsd_angstrom
        + w["clash"] * float(contacts_lt_1a)
        + w["affinity"] * float(affinity_kcal_mol)
        + w["envelope"] * envelope_nm3
    )
