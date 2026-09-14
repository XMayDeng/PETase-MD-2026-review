from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import math
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
INK = '#20242B'
MUTED = '#626A75'
L10 = '#33658A'
L20 = '#86BBD8'
SCOPE_ANY = '#30343B'
SCOPE_SIDECHAIN = '#4F81BD'
ANGLE_STRICT = '#C66A45'
ZERO = '#9299A2'
BACKGROUND = '#FFFFFF'

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def panel_label(ax: plt.Axes, label: str) -> mpl.text.Text:
    return ax.text(-0.02, 1.03, f'({label.lower()})', transform=ax.transAxes, ha='right', va='bottom', fontsize=8.5, fontweight='bold', clip_on=False, zorder=20)

def clean_axis(ax: plt.Axes) -> None:
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(False)

def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for source in (CASE_SOURCE, PET_SOURCE, SHIFT_SOURCE, ANGLE_SOURCE):
        if not source.is_file():
            raise FileNotFoundError(source)
    case_values = pd.read_csv(CASE_SOURCE)
    pet_summary = pd.read_csv(PET_SOURCE)
    prior_shifts = pd.read_csv(SHIFT_SOURCE)
    angle_data = pd.read_csv(ANGLE_SOURCE)
    native = angle_data[angle_data['criterion'] == 'gromacs_ADH_le_30'].copy()
    reconciliation_errors: list[float] = []
    for row in prior_shifts.itertuples(index=False):
        protein_name = 'CUT1' if row.protein_id == 'Pro00057' else 'HiC'
        for scope, prefix in (('any_residue_atom', ''), ('sidechain_only', 'sidechain_')):
            match = native[(native['protein_name'] == protein_name) & (native['residue_label'] == row.residue_label) & (native['metric'] == scope)]
            if len(match) != 1:
                raise ValueError(f'Missing native reconciliation row: {protein_name} {row.residue_label} {scope}')
            reconciliation_errors.append(abs(float(match.iloc[0]['delta_l20_minus_l10']) - float(getattr(row, f'{prefix}delta_l20_minus_l10'))))
            selected = match.iloc[0]
            if (int(selected['n_l10']), int(selected['n_l20'])) != (
                    int(row.n_l10), int(row.n_l20)) or bool(selected['ci_excludes_zero']) != bool(
                        getattr(row, f'{prefix}ci_excludes_zero')):
                raise ValueError(f'Primary sample-count/resolution mismatch: {row.residue_label} {scope}')
            angle_data.loc[match.index, ['delta_ci95_low', 'delta_ci95_high']] = [
                float(getattr(row, f'{prefix}delta_ci95_low')),
                float(getattr(row, f'{prefix}delta_ci95_high'))]
    max_error = max(reconciliation_errors)
    if max_error >= 1e-06:
        raise ValueError(f'Native H-bond shift reconciliation failed: {max_error}')
    native = angle_data[angle_data['criterion'] == 'gromacs_ADH_le_30'].copy()
    expected_cases = {('Pro00057', 'PET_L10'): 6, ('Pro00057', 'PET_L20'): 6, ('Pro00121', 'PET_L10'): 5, ('Pro00121', 'PET_L20'): 8}
    observed = case_values.groupby(['protein_id', 'pet_kind'])['case'].nunique().to_dict()
    if observed != expected_cases:
        raise ValueError(f'Unexpected trajectory-count map: {observed}')
    return (case_values, pet_summary, native, angle_data)

def draw_occupancy_panel(ax: plt.Axes, case_values: pd.DataFrame, pet_summary: pd.DataFrame, protein_id: str, protein_name: str, residues: list[str], occupancy_max: float, letter: str) -> None:
    styles = {'PET_L10': ('L10', L10, 'o', -0.14), 'PET_L20': ('L20', L20, 's', 0.14)}
    x = np.arange(len(residues), dtype=float)
    for pet_kind, (_, color, marker, offset) in styles.items():
        for xi, residue in zip(x, residues, strict=True):
            rows = case_values[(case_values['protein_id'] == protein_id) & (case_values['pet_kind'] == pet_kind) & (case_values['residue_label'] == residue)].sort_values('case')
            values = rows['any_hbond_occupancy'].to_numpy(float)
            jitter = np.linspace(-0.045, 0.045, len(values)) if len(values) > 1 else np.zeros(1)
            ax.scatter(xi + offset + jitter, values, s=10, marker=marker, color=color, edgecolor='white', linewidth=0.35, zorder=3)
            summary = pet_summary[(pet_summary['protein_id'] == protein_id) & (pet_summary['pet_kind'] == pet_kind) & (pet_summary['residue_label'] == residue)]
            if len(summary) != 1:
                raise ValueError(f'Missing H-bond occupancy summary: {protein_id} {pet_kind} {residue}')
            mean = float(summary.iloc[0]['mean_any_hbond_occupancy'])
            low = float(summary.iloc[0]['mean_any_hbond_occupancy_ci95_low'])
            high = float(summary.iloc[0]['mean_any_hbond_occupancy_ci95_high'])
            ax.errorbar(xi + offset, mean, yerr=np.array([[mean - low], [high - mean]]), fmt=marker, color=color, markerfacecolor=color, markeredgecolor=INK, markeredgewidth=0.5, markersize=4.4, elinewidth=0.9, capsize=1.8, zorder=5)
    counts = case_values[case_values['protein_id'] == protein_id].groupby('pet_kind')['case'].nunique()
    ax.set_title(f'{protein_name}  (L10 n={int(counts['PET_L10'])}; L20 n={int(counts['PET_L20'])})', loc='left', pad=5, fontweight='bold')
    ax.set_xticks(x, residues)
    ax.set_ylim(-0.015, occupancy_max)
    ax.set_xlabel('Candidate residue')
    clean_axis(ax)
    panel_label(ax, letter)

def draw_native_shift_panel(ax: plt.Axes, native: pd.DataFrame, protein_name: str, title: str, residues: list[str], bound: float, letter: str) -> None:
    scopes = [('any_residue_atom', 'Any residue atom', 'o', -0.13, SCOPE_ANY), ('sidechain_only', 'Side-chain only', 'D', 0.13, SCOPE_SIDECHAIN)]
    y = np.arange(len(residues))
    for scope, _, marker, offset, color in scopes:
        rows = native[(native['protein_name'] == protein_name) & (native['metric'] == scope)].set_index('residue_label').reindex(residues)
        values = rows['delta_l20_minus_l10'].to_numpy(float)
        low = rows['delta_ci95_low'].to_numpy(float)
        high = rows['delta_ci95_high'].to_numpy(float)
        ax.errorbar(values, y + offset, xerr=np.vstack([values - low, high - values]), fmt=marker, color=color, markerfacecolor=color, markeredgecolor=color, markeredgewidth=0.8, markersize=4.1, elinewidth=0.85, capsize=1.7, zorder=3)
    ax.axvline(0, color=ZERO, linewidth=0.8, linestyle=(0, (2.2, 2.2)), zorder=0)
    ax.set_yticks(y, residues)
    ax.invert_yaxis()
    ax.set_xlim(-bound, bound)
    ax.set_title(title, loc='left', pad=5, fontweight='bold')
    ax.set_xlabel('L20 − L10 mean trajectory occupancy')
    clean_axis(ax)
    panel_label(ax, letter)

def draw_angle_panel(ax: plt.Axes, angle_data: pd.DataFrame, scope: str, title: str, bound: float, letter: str) -> None:
    residue_order = [('CUT1', 'T62'), ('CUT1', 'T64'), ('CUT1', 'S67'), ('CUT1', 'N213'), ('HiC', 'T164'), ('HiC', 'T166')]
    criteria = [('gromacs_ADH_le_30', 'A–D–H ≤30°', 'o', SCOPE_ANY), ('DHA_ge_150', 'D–H–A ≥150°', 's', SCOPE_SIDECHAIN), ('DHA_ge_160', 'D–H–A ≥160°', '^', ANGLE_STRICT)]
    offsets = [-0.18, 0.0, 0.18]
    subset = angle_data[angle_data['metric'] == scope]
    y = np.arange(len(residue_order))
    for offset, (criterion, _, marker, color) in zip(offsets, criteria, strict=True):
        rows = []
        for protein, residue in residue_order:
            match = subset[(subset['protein_name'] == protein) & (subset['residue_label'] == residue) & (subset['criterion'] == criterion)]
            if len(match) != 1:
                raise ValueError(f'Missing S1 sensitivity row: {scope}, {protein} {residue}, {criterion}')
            rows.append(match.iloc[0])
        values = np.array([float(row['delta_l20_minus_l10']) for row in rows])
        low = np.array([float(row['delta_ci95_low']) for row in rows])
        high = np.array([float(row['delta_ci95_high']) for row in rows])
        ax.errorbar(values, y + offset, xerr=np.vstack([values - low, high - values]), fmt=marker, color=color, markerfacecolor='white', markeredgewidth=0.85, markersize=3.8, elinewidth=0.75, capsize=1.6, zorder=3)
    ax.axvline(0, color=ZERO, linewidth=0.8, linestyle=(0, (2.2, 2.2)), zorder=0)
    ax.set_yticks(y, [f'{('TfCut1' if protein == 'CUT1' else protein)} {residue}' for protein, residue in residue_order])
    ax.invert_yaxis()
    ax.set_xlim(-bound, bound)
    ax.set_title(title, loc='left', pad=5, fontweight='bold')
    ax.set_xlabel('L20 − L10 mean trajectory occupancy')
    clean_axis(ax)
    panel_label(ax, letter)

def build_figure() -> None:
    apply_style()
    case_values, pet_summary, native, angle_data = load_data()
    occupancy_max = max(0.8, float(case_values['any_hbond_occupancy'].max()) * 1.08)
    bound = float(np.nanmax(np.abs(angle_data[['delta_ci95_low', 'delta_ci95_high']].to_numpy())))
    # Preserve the accepted comparison scale across all primary/sensitivity
    # panels, exactly as in the current V9 source renderer.
    if bound >= 0.35:
        raise ValueError(f'Interval exceeds the registered +/-0.35 plotting range: {bound}')
    bound = 0.35
    fig = plt.figure(figsize=(7.0, 6.7))
    grid = fig.add_gridspec(3, 2, height_ratios=[1.12, 1.0, 1.22], left=0.13, right=0.985, top=0.955, bottom=0.085, wspace=0.34, hspace=0.54)
    axes = [fig.add_subplot(grid[row, col]) for row in range(3) for col in range(2)]
    draw_occupancy_panel(axes[0], case_values, pet_summary, 'Pro00057', 'TfCut1', ['T62', 'T64', 'S67', 'N213'], occupancy_max, 'a')
    draw_occupancy_panel(axes[1], case_values, pet_summary, 'Pro00121', 'HiC', ['T164', 'T166'], occupancy_max, 'b')
    axes[0].set_ylabel('Fraction of frames with\na direct H bond')
    draw_native_shift_panel(axes[2], native, 'CUT1', 'TfCut1', ['T62', 'T64', 'S67', 'N213'], bound, 'c')
    draw_native_shift_panel(axes[3], native, 'HiC', 'HiC', ['T164', 'T166'], bound, 'd')
    draw_angle_panel(axes[4], angle_data, 'any_residue_atom', 'Any residue atom', bound, 'e')
    draw_angle_panel(axes[5], angle_data, 'sidechain_only', 'Side-chain only', bound, 'f')
    occupancy_legend = [mlines.Line2D([], [], marker='o', color=L10, linestyle='none', label='L10'), mlines.Line2D([], [], marker='s', color=L20, linestyle='none', label='L20')]
    axes[0].legend(handles=occupancy_legend, loc='upper left', ncol=2, frameon=False, handletextpad=0.35, columnspacing=0.7, borderaxespad=0.2)
    scope_legend = [mlines.Line2D([], [], marker='o', color=SCOPE_ANY, linestyle='none', label='Any residue atom'), mlines.Line2D([], [], marker='D', color=SCOPE_SIDECHAIN, linestyle='none', label='Side-chain only')]
    axes[2].legend(handles=scope_legend, loc='lower left', bbox_to_anchor=(0.0, 1.01), ncol=2, frameon=False, handletextpad=0.4, borderaxespad=0.2)
    for axis in axes[2:4]:
        axis.set_title(axis.get_title(loc='left'), loc='left', y=1.15, pad=5, fontweight='bold')
        for text in axis.texts:
            if text.get_text() in ('(c)', '(d)'):
                text.set_y(1.21)
    criterion_legend = [mlines.Line2D([], [], color=SCOPE_ANY, marker='o', markerfacecolor='white', linestyle='none', label='A–D–H ≤30°'), mlines.Line2D([], [], color=SCOPE_SIDECHAIN, marker='s', markerfacecolor='white', linestyle='none', label='D–H–A ≥150°'), mlines.Line2D([], [], color=ANGLE_STRICT, marker='^', markerfacecolor='white', linestyle='none', label='D–H–A ≥160°')]
    axes[4].legend(handles=criterion_legend, loc='upper left', ncol=1, frameon=False, handletextpad=0.45, borderaxespad=0.2)
    return fig
CASE_SOURCE = ROOT / '02_postprocessing/results/trajectory_summaries/candidate_case_hbond_values_zero_filled.csv'
PET_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/candidate_pet_hbond_summary.csv'
SHIFT_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/candidate_hbond_shift_uncertainty.csv'
ANGLE_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/criterion_residue_shift_statistics.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S08', build_figure)
