"""Chain-length-aware cubic box protocol (spec §4.5).

Shared constants + lookup helpers, imported by:

* ``code/scripts/validate_complex_box.py`` — pre-MD reject when the solute
  envelope is too close to the periodic-image boundary.
* ``code/scripts/build_md_from_docking_best.py`` (Task 5) — sets the
  ``editconf -box <X X X> -bt cubic`` arguments per case.
* ``code/scripts/check_v1_production_gate.py`` (Task 7) — the PBC layer of
  the 4-layer production gate uses a chain-length-scaled jump threshold.

Spec authority: ``docs/plans/2026-05-03_md_workflow_v1_binding_attribution_spec.md``
§4.5 (box) and the production-gate PBC sub-layer.
"""
from __future__ import annotations

# Cubic-box edge length per chain length, in nanometres (spec §4.5 table).
# Values are based on PET-only 20 ns pre-equilibration end-to-end / Rg
# measurements and a +1.5 nm minimum margin (see ``DEFAULT_MARGIN_NM``).
COMPLEX_BOX_DIM_NM: dict[str, float] = {
    "PET_L2": 10.0,
    "PET_L4": 10.0,
    "PET_L10": 13.0,
    "PET_L20": 15.0,
}

# Production-gate PBC jump threshold per chain length (spec §4.5).
# Longer chains exhibit larger natural conformation jumps; a single
# fixed threshold would yield false positives on L20.
PBC_JUMP_THRESHOLD_NM: dict[str, float] = {
    "PET_L2": 1.0,
    "PET_L4": 1.0,
    "PET_L10": 1.5,
    "PET_L20": 2.0,
}

# Minimum required gap between the solute heavy-atom envelope and the
# nearest box face, in nanometres. Failing this gap triggers either a
# pre-MD reject (Task 2) or the R1 Phase B fallback retry loop (Task 5,
# PA-D).
DEFAULT_MARGIN_NM: float = 1.5

# Resnames excluded from the "solute envelope" calculation in
# validate_complex_box; everything else counts. Covers CHARMM (TIP3,
# SOD, CLA, POT, etc.) and the GROMACS canonical alternatives
# (SOL, NA, CL, K). Add new resnames here when a workflow stage
# introduces a new water/ion model.
SOLVENT_AND_ION_RESNAMES: frozenset[str] = frozenset({
    "SOL",   # GROMACS canonical water
    "TIP3",  # CHARMM TIP3P
    "TIPS3", # CHARMM-GUI variant
    "HOH",   # PDB water
    "WAT",   # alt water resname
    "NA",    # GROMACS sodium
    "CL",    # GROMACS chloride
    "SOD",   # CHARMM sodium
    "CLA",   # CHARMM chloride
    "K",     # potassium
    "POT",   # CHARMM potassium
    "MG",    # magnesium
    "CA",    # calcium (ambiguous with CA backbone atom name — context: this
             # list is matched against residue name, not atom name)
    "ZN",
})

SUPPORTED_PET_KINDS: frozenset[str] = frozenset(COMPLEX_BOX_DIM_NM.keys())


def box_dim_for(pet_kind: str) -> float:
    """Cubic edge length (nm) for the given PET chain length."""
    try:
        return COMPLEX_BOX_DIM_NM[pet_kind]
    except KeyError as exc:
        raise ValueError(
            f"Unknown pet_kind={pet_kind!r}; supported = "
            f"{sorted(COMPLEX_BOX_DIM_NM)}"
        ) from exc


def pbc_jump_threshold_for(pet_kind: str) -> float:
    """Production-gate PBC jump threshold (nm) for the given chain length."""
    try:
        return PBC_JUMP_THRESHOLD_NM[pet_kind]
    except KeyError as exc:
        raise ValueError(
            f"Unknown pet_kind={pet_kind!r}; supported = "
            f"{sorted(PBC_JUMP_THRESHOLD_NM)}"
        ) from exc


def is_solvent_or_ion(resname: str,
                      excluded: frozenset[str] | set[str] | None = None) -> bool:
    """Return ``True`` if ``resname`` should be excluded from the solute
    envelope (matches case-insensitively, trims whitespace)."""
    set_to_use = SOLVENT_AND_ION_RESNAMES if excluded is None else excluded
    return resname.strip().upper() in set_to_use
