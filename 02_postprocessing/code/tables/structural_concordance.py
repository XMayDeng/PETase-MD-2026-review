"""STR-EQ-01 numerical projection, historical Table S2 / current Table S3."""
from __future__ import annotations
import csv
from dataclasses import dataclass
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
SUMMARY_CSV=ROOT/"02_postprocessing/results/aggregate_statistics/protein_model_equivalence_summary.csv"
KEY_RESIDUE_CSV=ROOT/"02_postprocessing/results/aggregate_statistics/named_key_residue_comparison.csv"
EXPECTED_PROTEINS = {'Pro00057': 'TfCut1', 'Pro00062': 'LCC', 'Pro00075': 'FoCut5a', 'Pro00083': 'IsPETase', 'Pro00121': 'HiC', 'Pro00137': 'PHL7'}

PROTEIN_ORDER = tuple(EXPECTED_PROTEINS)

BACKBONE_FLAG_ANGSTROM = 1.0

SIDECHAIN_FLAG_ANGSTROM = 1.5

@dataclass(frozen=True)
class LocalFlag:
    residue_label: str
    metric: str
    value: float

def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))

def optional_float(value: str) -> float | None:
    return float(value) if value.strip() else None

def format_residue_label(label: str) -> str:
    return label[:3].title() + label[3:]

def collect_local_flags(rows: list[dict[str, str]]) -> dict[str, list[LocalFlag]]:
    flags = {protein_id: [] for protein_id in PROTEIN_ORDER}
    for row in rows:
        protein_id = row['protein_id']
        backbone = optional_float(row['backbone_rmsd_after_pocket_fit_angstrom'])
        sidechain = optional_float(row['sidechain_conformation_rmsd_after_residue_backbone_fit_angstrom'])
        if backbone is not None and backbone > BACKBONE_FLAG_ANGSTROM:
            flags[protein_id].append(LocalFlag(format_residue_label(row['production_label']), 'local backbone', backbone))
        if sidechain is not None and sidechain > SIDECHAIN_FLAG_ANGSTROM:
            flags[protein_id].append(LocalFlag(format_residue_label(row['production_label']), 'side-chain conformation', sidechain))
    return flags

def format_flags(flags: list[LocalFlag]) -> str:
    if not flags:
        return 'None'
    grouped: dict[str, list[LocalFlag]] = {}
    for flag in flags:
        grouped.setdefault(flag.metric, []).append(flag)
    parts = []
    for metric, metric_flags in grouped.items():
        residues = '/'.join((flag.residue_label for flag in metric_flags))
        values = '/'.join((f'{flag.value:.2f}' for flag in metric_flags))
        parts.append(f'{residues} {metric} ({values} A)')
    return '; '.join(parts)

def build_rows() -> list[dict[str, str]]:
    summary_rows = read_csv(SUMMARY_CSV)
    summary_by_id = {row['protein_id']: row for row in summary_rows}
    if set(summary_by_id) != set(EXPECTED_PROTEINS):
        raise ValueError('Structural summary does not contain exactly the expected six proteins')
    local_flags = collect_local_flags(read_csv(KEY_RESIDUE_CSV))
    output_rows = []
    for protein_id in PROTEIN_ORDER:
        row = summary_by_id[protein_id]
        if row['protein'] != EXPECTED_PROTEINS[protein_id]:
            raise ValueError(f'Unexpected protein name for {protein_id}: {row['protein']}')
        identity = float(row['sequence_identity_over_aligned_residues'])
        if identity != 1.0:
            raise ValueError(f'Mapped identity is not 1.0 for {protein_id}')
        output_rows.append({'enzyme': row['protein'], 'production_model_id': protein_id, 'experimental_reference': f'{row['experimental_pdb']} chain A', 'experimental_resolution_angstrom': f'{float(row['experimental_resolution_angstrom']):.2f}', 'aligned_residues_n': row['aligned_residue_pairs'], 'aligned_sequence_identity_percent': f'{identity * 100:.1f}', 'production_model_coverage_percent': f'{float(row['production_sequence_coverage']) * 100:.1f}', 'global_calpha_rmsd_angstrom': f'{float(row['global_calpha_rmsd_angstrom']):.3f}', 'pocket_backbone_rmsd_angstrom': f'{float(row['pocket_backbone_local_fit_rmsd_angstrom']):.3f}', 'pocket_sidechain_conformation_rmsd_median_angstrom': f'{float(row['pocket_sidechain_conformation_rmsd_median_angstrom']):.3f}', 'pocket_sidechain_conformation_rmsd_p90_angstrom': f'{float(row['pocket_sidechain_conformation_rmsd_p90_angstrom']):.3f}', 'maximum_triad_distance_delta_angstrom': f'{float(row['maximum_catalytic_geometry_delta_angstrom']):.3f}', 'prespecified_local_flag': format_flags(local_flags[protein_id])})
    return output_rows

def latex_escape(text: str) -> str:
    replacements = {'&': '\\&', '%': '\\%', '_': '\\_', '#': '\\#'}
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text

def latex_flag(text: str) -> str:
    if text == 'None':
        return 'None'
    text = latex_escape(text).replace('/', '/\\allowbreak{}')
    text = text.replace(' A)', ' \\AA)')
    return text
