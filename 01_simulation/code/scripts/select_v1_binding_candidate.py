#!/usr/bin/env python3
"""Select v1 binding-attribution candidate from 03_pet_alignment_screen output.

If the screening stage emitted ``queued_candidates`` for a short-MD retention
gate, use that queue first. Those candidates are the screener's explicit
"advance to MD gate" set. If no queue is available, fall back to the ranked
candidate list and choose a geometrically sane clash-free pose:

    1. prefer ``static_gate_pass == True``
    2. tie-break on lowest motif RMSD (best alignment quality)
    3. tie-break on best Vina affinity (most negative kcal/mol)

Output: copy the chosen `protein_pet_complex.pdb` to
``02_docking/best_rigidbody.pdb`` so ``run_build_from_docking.sh`` finds the
selected complex at the legacy HADDOCK path.

Usage::

    python select_v1_binding_candidate.py \\
        --screen-dir <runtime>/03_pet_alignment_screen \\
        --target-pdb  <runtime>/02_docking/best_rigidbody.pdb
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path


SUMMARY_FILE_CANDIDATES = (
    "screening_summary_multi_conformer.json",   # R1 Phase B aggregated (v1.1)
    "screening_summary.json",                    # single-conformer (legacy/back-compat)
)


def _load_summary(screen_dir: Path) -> dict:
    for name in SUMMARY_FILE_CANDIDATES:
        candidate = screen_dir / name
        if candidate.exists():
            return json.loads(candidate.read_text(encoding="utf-8"))
    raise SystemExit(
        f"select_v1_binding_candidate: no screening summary found under "
        f"{screen_dir} (looked for {SUMMARY_FILE_CANDIDATES})"
    )


def _selection_pool(summary: dict) -> tuple[list[dict], str, str]:
    """Return candidate records plus selector metadata.

    ``queued_candidates`` is the screener's explicit short-MD gate queue and
    must take precedence over the looser ranked list. Older summaries may not
    have a queue, so the ranked-list fallback is kept for compatibility.
    """
    queued = summary.get("queued_candidates")
    if queued:
        return (
            queued,
            "screening_short_md_gate_queue",
            "Screening emitted queued_candidates; select from the short-MD gate queue before production.",
        )
    for key in ("top_candidates", "ranked_candidates", "candidates", "all_candidates"):
        records = summary.get(key)
        if records:
            return (
                records,
                "static_screen_fallback",
                "No queued_candidates found; fallback selection uses static screen, motif RMSD, and affinity.",
            )
    return [], "none", "No candidate records found in screening summary."


def _candidate_records(summary: dict) -> list[dict]:
    return _selection_pool(summary)[0]


def _ranking_key(record: dict) -> tuple[bool, float]:
    """Sort key for v1.1 multi-conformer mode: prefer records with a valid
    ``ranking_score`` (set by the multi-conformer orchestrator); fall back to
    sentinel infinity for records without one."""
    score = record.get("ranking_score")
    if score is None:
        return (True, float("inf"))
    try:
        return (False, float(score))
    except (TypeError, ValueError):
        return (True, float("inf"))


def _group_by_direction(records: list[dict]) -> dict[str, list[dict]]:
    """Return {chain_extension_direction: [records]} for v1.1 multi-conformer
    selection. Records without the field are placed under ``"unknown"``."""
    groups: dict[str, list[dict]] = {}
    for r in records:
        key = r.get("chain_extension_direction") or "unknown"
        groups.setdefault(key, []).append(r)
    return groups


def select_per_direction(
    records: list[dict], *, max_severe_clashes: int,
) -> dict[str, dict | None]:
    """Pick the lowest-ranking-score clash-free candidate per chain-extension
    direction. Returns ``{"head-side": rec|None, "tail-side": ..., ...}``."""
    eligible = [r for r in records
                if _is_clash_free(r, max_severe_clashes) and r.get("complex_pdb")]
    groups = _group_by_direction(eligible)
    chosen: dict[str, dict | None] = {}
    for direction, recs in groups.items():
        if not recs:
            chosen[direction] = None
            continue
        recs_sorted = sorted(recs, key=_ranking_key)
        chosen[direction] = recs_sorted[0]
    return chosen


def _is_clash_free(record: dict, max_severe_clashes: int) -> bool:
    severe = record.get("contacts_lt_1a")
    if severe is None:
        return True  # screen did not compute clash count; trust upstream
    try:
        return int(severe) <= max_severe_clashes
    except (TypeError, ValueError):
        return True


def _sort_key(record: dict) -> tuple[bool, float, float]:
    """Lower is better. Sort priority:

    1. ``static_gate_pass == True`` first (screen itself flagged the candidate
       as geometrically sane: motif RMSD + clash + Ser-carbonyl thresholds).
    2. Lowest ``motif_rmsd_angstrom`` (best alignment quality between Vina
       PDBQT pose and CHARMM PET geometry).
    3. Lowest ``affinity_kcal_mol`` (most negative Vina score).

    Rationale (vs naive affinity-first ordering): Vina affinity differences
    < ~1 kcal/mol are within scoring-function noise; alignment-quality and
    clash-freeness are stronger signals for setup correctness in v1 binding-
    attribution. Earlier "affinity-first" sort picked a high-RMSD high-clash
    pose over a low-RMSD clash-free pose with marginally worse affinity.
    """
    static_pass = bool(record.get("static_gate_pass", False))
    motif_rmsd = record.get("motif_rmsd_angstrom")
    affinity = record.get("affinity_kcal_mol")
    return (
        not static_pass,  # False sorts before True; static_pass=True → not=False → first
        float(motif_rmsd) if motif_rmsd is not None else float("inf"),
        float(affinity) if affinity is not None else float("inf"),
    )


def _select(records: list[dict], max_severe_clashes: int) -> dict | None:
    eligible = [r for r in records if _is_clash_free(r, max_severe_clashes) and r.get("complex_pdb")]
    if not eligible:
        return None
    eligible.sort(key=_sort_key)
    return eligible[0]


def _candidate_meta(
    chosen: dict,
    src: Path,
    dst: Path,
    *,
    selector: str,
    selector_note: str | None = None,
) -> dict:
    return {
        "selected_pose": chosen.get("pose"),
        "selected_candidate": chosen.get("candidate"),
        "kabsch_direction": chosen.get("direction"),
        "chain_extension_direction": chosen.get("chain_extension_direction"),
        "conformer_idx": chosen.get("conformer_idx"),
        "window_start": chosen.get("window_start"),
        "window_end": chosen.get("window_end"),
        "ranking_score": chosen.get("ranking_score"),
        "post_kabsch_envelope_nm3": chosen.get("post_kabsch_envelope_nm3"),
        "affinity_kcal_mol": chosen.get("affinity_kcal_mol"),
        "motif_rmsd_angstrom": chosen.get("motif_rmsd_angstrom"),
        "contacts_lt_1a": chosen.get("contacts_lt_1a"),
        "contacts_lt_2a": chosen.get("contacts_lt_2a"),
        "contacts_lt_4a": chosen.get("contacts_lt_4a"),
        "ser_carbonyl_min_angstrom": chosen.get("ser_carbonyl_min_angstrom"),
        "source_complex_pdb": str(src),
        "target_pdb": str(dst),
        "selector": selector,
        "selector_note": selector_note or "",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--screen-dir", required=True, type=Path,
                        help="03_pet_alignment_screen directory containing screening_summary.json "
                             "(or screening_summary_multi_conformer.json in v1.1 mode)")
    parser.add_argument("--target-pdb", required=True, type=Path,
                        help="Destination path for the v1-selected protein_pet_complex.pdb. "
                             "In --per-direction mode the file stem is suffixed with "
                             "_<direction> (e.g. best_rigidbody_head-side.pdb) and 3 files are emitted.")
    parser.add_argument("--max-severe-clashes", type=int, default=50,
                        help="Reject candidates with strictly more than this many heavy-atom contacts "
                             "< 1.0 A. Note: Vina-Gasteiger vs CHARMM atom-center sub-Å differences "
                             "produce 4-10 'clashes' even for clean alignments; EM resolves these. "
                             "Default 50 is permissive — only catastrophic overlap (e.g. PET buried in "
                             "protein core) is rejected. Lower for stricter screening.")
    parser.add_argument("--meta", type=Path, default=None,
                        help="Optional path to write selection metadata as JSON")
    parser.add_argument("--per-direction", action="store_true",
                        help="v1.1 R1 Phase B mode: emit 3 best.pdb files (one per "
                             "chain_extension_direction) ranked by ranking_score.")
    args = parser.parse_args(argv)

    summary = _load_summary(args.screen_dir)
    records, selector_name, selector_note = _selection_pool(summary)
    if not records:
        print(
            "select_v1_binding_candidate: no candidate records in screening summary",
            file=sys.stderr,
        )
        return 11

    if args.per_direction:
        chosen_per_dir = select_per_direction(
            records, max_severe_clashes=args.max_severe_clashes
        )
        if not any(chosen_per_dir.values()):
            print(
                f"select_v1_binding_candidate: no clash-free candidate found in any direction "
                f"(severe_clashes ≤ {args.max_severe_clashes}) among {len(records)} records",
                file=sys.stderr,
            )
            return 12

        args.target_pdb.parent.mkdir(parents=True, exist_ok=True)
        stem = args.target_pdb.stem
        suffix = args.target_pdb.suffix
        per_direction_meta: dict[str, dict | None] = {}
        for direction, chosen in chosen_per_dir.items():
            if chosen is None:
                per_direction_meta[direction] = None
                continue
            src = Path(chosen["complex_pdb"])
            if not src.exists():
                print(
                    f"select_v1_binding_candidate: chosen complex_pdb missing for "
                    f"direction={direction}: {src}", file=sys.stderr)
                per_direction_meta[direction] = None
                continue
            dst = args.target_pdb.with_name(f"{stem}_{direction}{suffix}")
            shutil.copy2(src, dst)
            per_direction_meta[direction] = _candidate_meta(
                chosen, src, dst,
                selector="v1_binding_per_direction_ranking_score",
                selector_note="v1.1 long-chain per-direction selection ranked by ranking_score.",
            )

        report = {
            "mode": "per_direction",
            "max_severe_clashes": args.max_severe_clashes,
            "per_direction": per_direction_meta,
        }
        print(json.dumps(report, indent=2))
        if args.meta is not None:
            args.meta.parent.mkdir(parents=True, exist_ok=True)
            args.meta.write_text(json.dumps(report, indent=2), encoding="utf-8")
        return 0

    # Legacy single-best-pose path (kept for backward compatibility).
    chosen = _select(records, args.max_severe_clashes)
    if chosen is None:
        print(
            f"select_v1_binding_candidate: no clash-free candidate "
            f"(severe_clashes ≤ {args.max_severe_clashes}) found among {len(records)} records",
            file=sys.stderr,
        )
        return 12

    src = Path(chosen["complex_pdb"])
    if not src.exists():
        print(f"select_v1_binding_candidate: chosen complex_pdb missing: {src}", file=sys.stderr)
        return 13

    args.target_pdb.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, args.target_pdb)

    selection_meta = _candidate_meta(
        chosen, src, args.target_pdb,
        selector=f"v1_binding_{selector_name}",
        selector_note=selector_note,
    )
    print(json.dumps(selection_meta, indent=2))
    if args.meta is not None:
        args.meta.parent.mkdir(parents=True, exist_ok=True)
        args.meta.write_text(json.dumps(selection_meta, indent=2), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
