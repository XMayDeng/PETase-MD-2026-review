from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
PROTEIN_ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']
PAIR_PROTEINS = ['IsPETase', 'TfCut1', 'LCC', 'PHL7']
PET_ORDER = ['PET_L10', 'PET_L20']
PROTEIN_COLORS = {'IsPETase': '#4472C4', 'TfCut1': '#ED7D31', 'LCC': '#98C97A', 'FoCut5a': '#E45A9D', 'HiC': '#7FAED3', 'PHL7': '#70AD47'}
LINE_STYLES = {'IsPETase': '-', 'TfCut1': '--', 'LCC': '-.', 'FoCut5a': ':', 'HiC': (0, (5, 1)), 'PHL7': (0, (1, 1))}
PAIR_COLORS = {'pair_site_1': '#3E7CB8', 'pair_site_2': '#D47732'}
INK = '#20242B'
BACKGROUND = '#FFFFFF'
MARKER_LINE = '#B69A2A'

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.3, 'axes.labelsize': 6.9, 'xtick.labelsize': 6.1, 'ytick.labelsize': 6.1, 'legend.fontsize': 5.8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def panel_label(ax: plt.Axes, label: str) -> mpl.text.Text:
    return ax.text(0.0, 1.035, f'({label.lower()})', transform=ax.transAxes, ha='left', va='bottom', fontsize=8.5, fontweight='bold', clip_on=False, zorder=20)

def align_panel_label_to_yticks(fig: plt.Figure, ax: plt.Axes, label: mpl.text.Text) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    visible = [tick for tick in ax.get_yticklabels() if tick.get_visible() and tick.get_text().strip()]
    if not visible:
        return
    left_display = min((tick.get_window_extent(renderer=renderer).x0 for tick in visible))
    left_axes = ax.transAxes.inverted().transform((left_display, ax.bbox.y0))[0]
    label.set_x(left_axes)

def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sources = (PROFILE_SOURCE, TRAJECTORY_SOURCE, SUMMARY_SOURCE, PAIR_SOURCE)
    for source in sources:
        if not source.is_file():
            raise FileNotFoundError(source)
    profile = pd.read_csv(PROFILE_SOURCE)
    trajectories = pd.read_csv(TRAJECTORY_SOURCE)
    summary = pd.read_csv(SUMMARY_SOURCE)
    pair = pd.read_csv(PAIR_SOURCE)
    expected_source_rows = {'profile': 66, 'trajectories': 71, 'summary': 12, 'pair': 8}
    observed_source_rows = {'profile': len(profile), 'trajectories': len(trajectories), 'summary': len(summary), 'pair': len(pair)}
    if observed_source_rows != expected_source_rows:
        raise ValueError(f'Unexpected Figure S8 source dimensions: {observed_source_rows}')
    profile_counts = profile.groupby('protein').size().reindex(PROTEIN_ORDER)
    if not profile_counts.eq(11).all():
        raise ValueError('Expected 11 mapped profile positions per protein')
    if sorted(profile['relative_position'].unique()) != list(range(-3, 8)):
        raise ValueError('Unexpected relative-position window')
    display_trajectories = trajectories[trajectories['pet_kind'].isin(PET_ORDER)].copy()
    if len(display_trajectories) != 64:
        raise ValueError(f'Expected 64 retained L10/L20 trajectories, found {len(display_trajectories)}')
    summary_rows = summary[summary['metric'] == 'wloop_mean_rmsf_nm'].set_index('protein').reindex(PROTEIN_ORDER)
    if len(summary_rows) != 6 or summary_rows.isna().any().any():
        raise ValueError('Incomplete protein-level W-loop summary')
    observed_counts = display_trajectories.groupby(['protein', 'pet_kind']).size().unstack('pet_kind').reindex(index=PROTEIN_ORDER, columns=PET_ORDER)
    if not np.array_equal(observed_counts['PET_L10'].to_numpy(int), summary_rows['n_l10'].to_numpy(int)) or not np.array_equal(observed_counts['PET_L20'].to_numpy(int), summary_rows['n_l20'].to_numpy(int)):
        raise ValueError('Trajectory counts disagree with the registered summaries')
    pair_index = set(zip(pair['protein'], pair['site_role'], strict=True))
    expected_pair_index = {(protein, role) for protein in PAIR_PROTEINS for role in ('pair_site_1', 'pair_site_2')}
    if pair_index != expected_pair_index:
        raise ValueError('Incomplete audited pair-site summary')
    return (profile, display_trajectories, summary_rows, pair)

def draw_profile_panel(ax: plt.Axes, profile: pd.DataFrame) -> mpl.text.Text:
    for protein in PROTEIN_ORDER:
        group = profile[profile['protein'] == protein].sort_values('relative_position')
        ax.plot(group['relative_position'], group['length_balanced_mean_rmsf_nm'], color=PROTEIN_COLORS[protein], linestyle=LINE_STYLES[protein], linewidth=1.2, marker='o', markersize=2.5, label=protein)
    ax.axvline(0, color=MARKER_LINE, linestyle='--', linewidth=0.85, zorder=0)
    ax.set_xlim(-3.5, 7.5)
    ax.set_ylim(0.015, 0.395)
    ax.set_xticks(range(-3, 8))
    ax.set_xlabel('Position relative to mapped loop marker')
    ax.set_ylabel('C$\\alpha$ RMSF (nm)')
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(loc='lower center', bbox_to_anchor=(0.56, 1.005), ncol=6, frameon=False, handlelength=1.45, columnspacing=0.75, borderaxespad=0)
    return panel_label(ax, 'a')

def draw_trajectory_panel(ax: plt.Axes, trajectories: pd.DataFrame, summary_rows: pd.DataFrame) -> mpl.text.Text:
    x_positions = np.arange(len(PROTEIN_ORDER), dtype=float)
    for protein_index, protein in enumerate(PROTEIN_ORDER):
        group = trajectories[trajectories['protein'] == protein]
        for pet_kind, marker, offset, filled in (('PET_L10', 'o', -0.08, False), ('PET_L20', '^', 0.08, True)):
            values = group.loc[group['pet_kind'] == pet_kind, 'wloop_mean_rmsf_nm'].to_numpy(float)
            jitter = np.linspace(-0.025, 0.025, len(values)) if len(values) > 1 else np.array([0.0])
            ax.scatter(protein_index + offset + jitter, values, s=17, marker=marker, facecolor=PROTEIN_COLORS[protein] if filled else BACKGROUND, edgecolor=PROTEIN_COLORS[protein], linewidth=0.75, alpha=0.88, zorder=3)
        row = summary_rows.loc[protein]
        center = float(row['length_balanced_mean_nm'])
        ax.errorbar(protein_index, center, yerr=[[center - float(row['bootstrap_ci_low_nm'])], [float(row['bootstrap_ci_high_nm']) - center]], fmt='D', color=INK, markerfacecolor=BACKGROUND, markeredgewidth=0.8, markersize=4.0, elinewidth=0.9, capsize=2.0, zorder=5)
    ax.set_xlim(-0.45, len(PROTEIN_ORDER) - 0.55)
    ax.set_ylim(0.02, 0.52)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(PROTEIN_ORDER, rotation=27, ha='right')
    ax.set_ylabel('Mean loop C$\\alpha$ RMSF (nm)')
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=[mlines.Line2D([], [], marker='o', linestyle='none', markerfacecolor=BACKGROUND, markeredgecolor='#59616B', label='PET L10'), mlines.Line2D([], [], marker='^', linestyle='none', markerfacecolor='#8E959E', markeredgecolor='#59616B', label='PET L20'), mlines.Line2D([], [], marker='D', linestyle='-', color=INK, markerfacecolor=BACKGROUND, label='Balanced mean, 95% CI')], loc='lower center', bbox_to_anchor=(0.59, 1.08), ncol=3, frameon=False, handlelength=1.25, columnspacing=0.7, borderaxespad=0)
    return panel_label(ax, 'b')

def draw_pair_panel(ax: plt.Axes, pair: pd.DataFrame) -> mpl.text.Text:
    x_positions = np.arange(len(PAIR_PROTEINS), dtype=float)
    style_rows = (('pair_site_1', 'S/H-like site', 'o', -0.11), ('pair_site_2', 'I/F-like site', 's', 0.11))
    text_effects = [mpatheffects.Stroke(linewidth=1.5, foreground=BACKGROUND), mpatheffects.Normal()]
    for role, label, marker, offset in style_rows:
        rows = pair[pair['site_role'] == role].set_index('protein').reindex(PAIR_PROTEINS)
        centers = rows['length_balanced_mean_nm'].to_numpy(float)
        low = rows['bootstrap_ci_low_nm'].to_numpy(float)
        high = rows['bootstrap_ci_high_nm'].to_numpy(float)
        ax.errorbar(x_positions + offset, centers, yerr=np.vstack([centers - low, high - centers]), fmt=marker, color=PAIR_COLORS[role], markerfacecolor=BACKGROUND, markeredgewidth=0.9, markersize=4.1, elinewidth=0.85, capsize=2.0, label=label, zorder=3)
        for x_value, row in zip(x_positions + offset, rows.itertuples(), strict=True):
            ax.text(x_value, float(row.bootstrap_ci_high_nm) + 0.009, str(row.residue_label), ha='center', va='bottom', fontsize=5.7, rotation=45, color=PAIR_COLORS[role], path_effects=text_effects, clip_on=False)
    ax.set_xlim(-0.45, len(PAIR_PROTEINS) - 0.55)
    ax.set_ylim(0.02, 0.45)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(PAIR_PROTEINS, rotation=27, ha='right')
    ax.set_ylabel('Pair-site C$\\alpha$ RMSF (nm)')
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(loc='lower center', bbox_to_anchor=(0.6, 1.08), ncol=2, frameon=False, handlelength=1.1, columnspacing=0.8, borderaxespad=0)
    return panel_label(ax, 'c')

def build_figure() -> None:
    apply_style()
    profile, trajectories, summary_rows, pair = load_data()
    fig = plt.figure(figsize=(7.0, 4.45))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], left=0.085, right=0.985, top=0.94, bottom=0.14, hspace=0.5, wspace=0.3)
    ax_profile = fig.add_subplot(grid[0, :])
    ax_trajectory = fig.add_subplot(grid[1, 0])
    ax_pair = fig.add_subplot(grid[1, 1])
    labels = [(ax_profile, draw_profile_panel(ax_profile, profile)), (ax_trajectory, draw_trajectory_panel(ax_trajectory, trajectories, summary_rows)), (ax_pair, draw_pair_panel(ax_pair, pair))]
    for ax, label in labels:
        align_panel_label_to_yticks(fig, ax, label)
    return fig
PROFILE_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/protein_wloop_profile_summary.csv'
TRAJECTORY_SOURCE = ROOT / '02_postprocessing/results/trajectory_summaries/trajectory_wloop_summary.csv'
SUMMARY_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/protein_wloop_metric_summary.csv'
PAIR_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/mapped_pair_site_summary.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S09', build_figure)
