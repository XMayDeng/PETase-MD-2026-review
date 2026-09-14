from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from revised_figure_layout import align_panel_labels_to_yticklabels, contrasting_text_color, panel_label
PROTEIN_ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']
CLASS_ORDER = ['Catalytic/oxyanion', 'W-loop/flexible or rigid pair', 'Cleft/rim patch', 'Flap/binding loop', 'Subsite/hotspot/stability', 'Other/candidate']
CLASS_DISPLAY = {'Catalytic/oxyanion': 'Catalytic/oxyanion', 'W-loop/flexible or rigid pair': 'W-loop/rigid pair', 'Cleft/rim patch': 'Cleft/rim patch', 'Flap/binding loop': 'Flap/binding loop', 'Subsite/hotspot/stability': 'Subsite/hotspot', 'Other/candidate': 'Other/candidate'}
CLASS_COLORS = {'Catalytic/oxyanion': '#E07B54', 'W-loop/flexible or rigid pair': '#E8C547', 'Cleft/rim patch': '#5B9BD5', 'Flap/binding loop': '#D98EC7', 'Subsite/hotspot/stability': '#70AD47', 'Other/candidate': '#A6A6A6'}
CLASS_HATCHES = {'Catalytic/oxyanion': '////', 'W-loop/flexible or rigid pair': '....', 'Cleft/rim patch': '\\\\', 'Flap/binding loop': 'xxxx', 'Subsite/hotspot/stability': '++++', 'Other/candidate': ''}
INK = '#20242B'
MUTED = '#626A75'
BACKGROUND = '#FFFFFF'
from data_checks import require_complete

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'axes.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'hatch.linewidth': 0.35})

def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    for path in (TOP8_PATH, TOP20_PATH):
        if not path.is_file():
            raise FileNotFoundError(path)
    top8 = pd.read_csv(TOP8_PATH)
    top8 = top8[top8['is_top8'].astype(bool)].copy()
    counts = top8.groupby('protein_short').size().reindex(PROTEIN_ORDER)
    if not counts.eq(8).all():
        raise ValueError(f'Expected eight residues per enzyme: {counts.to_dict()}')
    top20 = pd.read_csv(TOP20_PATH)
    composition = require_complete(top20.pivot(index='protein_short', columns='landmark_class', values='contact_mass_percent').reindex(index=PROTEIN_ORDER, columns=CLASS_ORDER))
    if not np.allclose(composition.sum(axis=1).to_numpy(float), 100.0, atol=1e-06):
        raise ValueError('Top-20 contact-mass composition must sum to 100% per enzyme')
    return (top8, top20)

def draw_composition(ax: plt.Axes, top20: pd.DataFrame) -> mpl.text.Text:
    pivot = require_complete(top20.pivot(index='protein_short', columns='landmark_class', values='contact_mass_percent').reindex(index=PROTEIN_ORDER, columns=CLASS_ORDER))
    first = top20.groupby('protein_short').first()
    labels = []
    for protein in PROTEIN_ORDER:
        row = first.loc[protein]
        warn = '*' if int(row['stability_warn_n']) == int(row['retained_n']) else ''
        labels.append(f'{protein}{warn}  (n={int(row['retained_n'])})')
    y = np.arange(len(PROTEIN_ORDER), dtype=float)
    left = np.zeros(len(PROTEIN_ORDER), dtype=float)
    for cls in CLASS_ORDER:
        values = pivot[cls].to_numpy(float)
        bars = ax.barh(y, values, left=left, height=0.61, color=CLASS_COLORS[cls], edgecolor='#4A4F57', linewidth=0.35)
        for bar in bars:
            bar.set_hatch(CLASS_HATCHES[cls])
        for yi, (value, start) in enumerate(zip(values, left)):
            if value >= 10.0:
                text_color = contrasting_text_color(bars[yi].get_facecolor(), dark=INK)
                label = ax.text(start + value / 2.0, yi, f'{value:.0f}', ha='center', va='center', fontsize=6.0, fontweight='bold', color=text_color)
                if bars[yi].get_hatch():
                    label.set_path_effects([mpatheffects.withStroke(linewidth=0.9, foreground=bars[yi].get_facecolor())])
        left += values
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontweight='bold')
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel('Share of pooled top-20 contact mass (%)', labelpad=2)
    ax.spines[['top', 'right']].set_visible(False)
    return panel_label(ax, 'a', y=1.025)

def draw_residue_panel(ax: plt.Axes, group: pd.DataFrame, panel_letter: str, show_xlabel: bool) -> mpl.text.Text:
    group = group.sort_values('rank_within_protein')
    protein = str(group.iloc[0]['protein_short'])
    pdb = str(group.iloc[0]['pdb_id'])
    retained_n = int(group.iloc[0]['n_retained'])
    stability_pass_n = int(group.iloc[0]['stability_pass_n'])
    warn = '*' if stability_pass_n < retained_n else ''
    y = np.arange(len(group), dtype=float)
    bars = ax.barh(y, group['mean_contact_fraction'].to_numpy(float), color=[CLASS_COLORS[cls] for cls in group['landmark_class']], edgecolor='#4A4F57', linewidth=0.35, height=0.68)
    for bar, cls in zip(bars, group['landmark_class']):
        bar.set_hatch(CLASS_HATCHES[cls])
    ax.set_yticks(y)
    ax.set_yticklabels(group['residue_label'], fontsize=7.0)
    ax.invert_yaxis()
    ax.axvline(0.5, color='#737A84', linewidth=0.6, linestyle=(0, (3, 2)), zorder=0)
    ax.set_xlim(0, 1.58)
    ax.spines['bottom'].set_bounds(0, 1)
    ax.set_xticks([0.0, 0.5, 1.0])
    for yi, row in enumerate(group.itertuples(index=False)):
        label = ax.text(1.04, yi, f'{row.mean_contact_fraction:.2f}; {int(row.n_cases_ge_0_5)}/{int(row.n_retained)}', ha='left', va='center', fontsize=7.0)
        label.set_path_effects([mpatheffects.withStroke(linewidth=1.15, foreground=BACKGROUND)])
    header = panel_label(ax, panel_letter, y=1.02)
    ax.text(0.5, 1.02, f'{protein} ({pdb}){warn}', transform=ax.transAxes, ha='center', va='bottom', fontsize=7.5, fontweight='bold')
    if show_xlabel:
        ax.set_xlabel('Mean trajectory\ncontact fraction')
    return header

def build_figure() -> None:
    apply_style()
    top8, top20 = load_inputs()
    fig = plt.figure(figsize=(7.0, 6.7))
    outer = fig.add_gridspec(3, 1, height_ratios=[0.95, 2.27, 0.23], left=0.135, right=0.985, top=0.95, bottom=0.055, hspace=0.44)
    residue_grid = outer[1, 0].subgridspec(2, 3, hspace=0.42, wspace=0.36)
    ax_a = fig.add_subplot(outer[0, 0])
    residue_axes = [fig.add_subplot(residue_grid[row, col]) for row in range(2) for col in range(3)]
    ax_legend = fig.add_subplot(outer[2, 0])
    panel_headers = [(ax_a, draw_composition(ax_a, top20))]
    composition_box = ax_a.get_position()
    ax_a.set_position([0.18, composition_box.y0, composition_box.x1 - 0.18, composition_box.height])
    for index, (ax, protein) in enumerate(zip(residue_axes, PROTEIN_ORDER)):
        header = draw_residue_panel(ax, top8[top8['protein_short'] == protein], chr(ord('b') + index), show_xlabel=False)
        panel_headers.append((ax, header))
    fig.supxlabel('Mean trajectory contact fraction', x=0.5, y=0.145, ha='center', va='center', fontsize=7.5)
    legend_handles = [mpatches.Patch(facecolor=CLASS_COLORS[cls], edgecolor='#4A4F57', linewidth=0.4, hatch=CLASS_HATCHES[cls], label=CLASS_DISPLAY[cls]) for cls in CLASS_ORDER]
    ax_legend.legend(handles=legend_handles, loc='center', bbox_to_anchor=(0.0, 0.0, 1.0, 1.0), bbox_transform=ax_legend.transAxes, mode='expand', ncol=3, frameon=False, fontsize=6.5, columnspacing=0.7, handlelength=0.95, handletextpad=0.35, borderaxespad=0.0)
    ax_legend.set_axis_off()
    fig._publication_center_groups = [residue_axes, [ax_legend]]
    fig._publication_panel_legend_layout = {'key': 'Figure_03', 'panel': ax_a, 'lower': residue_axes, 'legend': ax_legend.get_legend()}
    align_panel_labels_to_yticklabels(fig, panel_headers)
    return fig
TOP8_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/pooled_top8_contacts_unfiltered.csv'
TOP20_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/pooled_top20_contact_class_composition_unfiltered.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_03', build_figure)
