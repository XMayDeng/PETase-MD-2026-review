from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
FOCUS_RESIDUES = [28, 29, 70, 164, 165, 166]
TARGET_ORDER = ['FoCut5a', 'IsPETase', 'TfCut1', 'LCC', 'PHL7']
CONTACT_ORDER = ['HiC', *TARGET_ORDER]
CONDITION_ORDER = ['PET_L10', 'PET_L20']
INK = '#20242B'
MUTED = '#59616B'
BACKGROUND = '#FFFFFF'
PRODUCTION_COLOR = '#3E7CB8'
EXPERIMENTAL_COLOR = '#D47732'
MAPPING_COLORS = {'high_confidence': '#9DBDD7', 'region_supported': '#B8D6A3', 'topology_only': '#EBCB5A', 'ambiguous_or_unstable': '#D9DCE0', 'no_stable_counterpart': '#D9DCE0'}

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.2, 'axes.labelsize': 6.8, 'xtick.labelsize': 6.0, 'ytick.labelsize': 6.0, 'legend.fontsize': 5.8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

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

def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for source in (CE_SOURCE, MAPPING_SOURCE, CONTACT_SOURCE):
        if not source.is_file():
            raise FileNotFoundError(source)
    ce = pd.read_csv(CE_SOURCE)
    mapping = pd.read_csv(MAPPING_SOURCE)
    contact = pd.read_csv(CONTACT_SOURCE)
    if len(ce) != 10:
        raise ValueError(f'Expected 10 pairwise CE rows, found {len(ce)}')
    if len(mapping) != 35:
        raise ValueError(f'Expected 35 mapping rows, found {len(mapping)}')
    if len(contact) != 126:
        raise ValueError(f'Expected 126 contact-summary rows, found {len(contact)}')
    ce_index = set(zip(ce['target_protein'], ce['coordinate_source'], strict=True))
    expected_ce = {(protein, source) for protein in TARGET_ORDER for source in ('production', 'experimental')}
    if ce_index != expected_ce:
        raise ValueError('Pairwise CE table does not match the registered panel')
    focus_mapping = mapping[mapping['hic_residue'].isin(FOCUS_RESIDUES) & mapping['target_protein'].isin(TARGET_ORDER)]
    if len(focus_mapping) != len(FOCUS_RESIDUES) * len(TARGET_ORDER):
        raise ValueError('Incomplete requested-focus structural mapping')
    if not focus_mapping['analysis_role'].eq('requested_focus').all():
        raise ValueError('Figure S6 focus must contain only requested-focus sites')
    focus_contact = contact[contact['hic_residue'].isin(FOCUS_RESIDUES) & contact['protein'].isin(CONTACT_ORDER) & contact['pet_kind'].isin(CONDITION_ORDER)]
    expected_rows = len(FOCUS_RESIDUES) * len(CONTACT_ORDER) * len(CONDITION_ORDER)
    if len(focus_contact) != expected_rows:
        raise ValueError('Incomplete L10/L20 contact summary for requested-focus sites')
    positive_sample_sizes = focus_contact.loc[focus_contact['n_trajectories'] > 0].groupby(['protein', 'pet_kind'])['n_trajectories'].nunique()
    if not positive_sample_sizes.eq(1).all():
        raise ValueError('Inconsistent nonzero retained denominator within a protein/length group')
    return (ce, mapping, contact)

def draw_ce_panel(ax: plt.Axes, ce: pd.DataFrame) -> mpl.text.Text:
    pivot = ce.pivot(index='target_protein', columns='coordinate_source', values='ce_rms_angstrom').reindex(TARGET_ORDER)
    y_positions = np.arange(len(TARGET_ORDER))
    for y, protein in enumerate(TARGET_ORDER):
        ax.plot([pivot.loc[protein, 'production'], pivot.loc[protein, 'experimental']], [y, y], color='#9AA1AA', linewidth=0.75, zorder=1)
    ax.scatter(pivot['production'], y_positions, s=25, marker='o', facecolor=BACKGROUND, edgecolor=PRODUCTION_COLOR, linewidth=1.0, zorder=3)
    ax.scatter(pivot['experimental'], y_positions, s=25, marker='s', facecolor=BACKGROUND, edgecolor=EXPERIMENTAL_COLOR, linewidth=1.0, zorder=3)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(TARGET_ORDER)
    ax.invert_yaxis()
    ax.set_xlim(0.5, 5.05)
    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xlabel('CE RMSD to HiC ($\\AA$)')
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(handles=[mlines.Line2D([], [], marker='o', linestyle='none', markerfacecolor=BACKGROUND, markeredgecolor=PRODUCTION_COLOR, markeredgewidth=1.0, label='Production model'), mlines.Line2D([], [], marker='s', linestyle='none', markerfacecolor=BACKGROUND, markeredgecolor=EXPERIMENTAL_COLOR, markeredgewidth=1.0, label='Experimental PDB')], loc='upper right', ncol=1, handletextpad=0.4, columnspacing=0.8, borderaxespad=0.2, frameon=False)
    return panel_label(ax, 'a')

def draw_mapping_panel(ax: plt.Axes, mapping: pd.DataFrame) -> mpl.text.Text:
    subset = mapping[mapping['hic_residue'].isin(FOCUS_RESIDUES) & mapping['target_protein'].isin(TARGET_ORDER)].copy()
    lookup = subset.set_index(['hic_residue', 'target_protein'])
    for row_index, residue in enumerate(FOCUS_RESIDUES):
        for column_index, protein in enumerate(TARGET_ORDER):
            row = lookup.loc[residue, protein]
            mapping_class = str(row['mapping_class'])
            ax.add_patch(mpatches.Rectangle((column_index - 0.5, row_index - 0.5), 1, 1, facecolor=MAPPING_COLORS[mapping_class], edgecolor=BACKGROUND, linewidth=0.75))
            label = str(row['consensus_label']) if pd.notna(row['consensus_label']) else '—'
            ax.text(column_index, row_index, label, ha='center', va='center', fontsize=6.2, fontweight='bold' if mapping_class == 'high_confidence' else 'normal')
    ax.set_xlim(-0.5, len(TARGET_ORDER) - 0.5)
    ax.set_ylim(len(FOCUS_RESIDUES) - 0.5, -0.5)
    ax.set_xticks(range(len(TARGET_ORDER)))
    ax.set_xticklabels(TARGET_ORDER)
    ax.tick_params(axis='x', top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=2)
    ax.set_yticks(range(len(FOCUS_RESIDUES)))
    labels = subset.drop_duplicates('hic_residue').set_index('hic_residue')['hic_label'].reindex(FOCUS_RESIDUES)
    ax.set_yticklabels([f'HiC {label}' for label in labels])
    ax.tick_params(axis='y', length=0, pad=2)
    legend_handles = [mpatches.Patch(facecolor=MAPPING_COLORS[mapping_class], edgecolor='#59616B', linewidth=0.35, label=label) for mapping_class, label in (('high_confidence', 'High confidence'), ('region_supported', 'Region supported'), ('topology_only', 'Topology only'), ('ambiguous_or_unstable', 'Unresolved'))]
    ax.legend(handles=legend_handles, loc='upper center', bbox_to_anchor=(0.5, -0.075), ncol=4, frameon=False, handlelength=0.9, columnspacing=0.75, borderaxespad=0)
    return panel_label(ax, 'b')

def draw_contact_panel(ax: plt.Axes, contact: pd.DataFrame, mapping: pd.DataFrame, pet_kind: str, panel: str, show_ylabels: bool) -> tuple[mpl.image.AxesImage, mpl.text.Text]:
    subset = contact[contact['hic_residue'].isin(FOCUS_RESIDUES) & contact['protein'].isin(CONTACT_ORDER) & (contact['pet_kind'] == pet_kind)].copy()
    values = subset.pivot(index='hic_residue', columns='protein', values='mean_contact_fraction').reindex(index=FOCUS_RESIDUES, columns=CONTACT_ORDER)
    cmap = mcolors.LinearSegmentedColormap.from_list('contact_use', ['#F6F1F8', '#B9A7C7', '#604572'])
    cmap.set_bad(BACKGROUND)
    image = ax.imshow(values.to_numpy(float), vmin=0.0, vmax=1.0, cmap=cmap, aspect='auto')
    for row_index in range(values.shape[0]):
        for column_index in range(values.shape[1]):
            value = values.iloc[row_index, column_index]
            ax.text(column_index, row_index, '—' if pd.isna(value) else f'{float(value):.2f}', ha='center', va='center', fontsize=5.8, fontweight='bold' if pd.notna(value) else 'normal', color='white' if pd.notna(value) and float(value) >= 0.6 else INK)
    ax.set_xticks(range(len(CONTACT_ORDER)))
    ax.set_xticklabels(CONTACT_ORDER, rotation=25, ha='right', rotation_mode='anchor')
    ax.set_yticks(range(len(FOCUS_RESIDUES)))
    hic_labels = mapping.drop_duplicates('hic_residue').set_index('hic_residue')['hic_label'].reindex(FOCUS_RESIDUES)
    if show_ylabels:
        ax.set_yticklabels([f'HiC {label}' for label in hic_labels])
    else:
        ax.set_yticklabels([])
    ax.tick_params(axis='both', which='major', length=0, pad=2)
    ax.set_xticks(np.arange(-0.5, len(CONTACT_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(FOCUS_RESIDUES), 1), minor=True)
    ax.grid(which='minor', color=BACKGROUND, linewidth=0.75)
    ax.tick_params(which='minor', bottom=False, left=False)
    ax.text(0.5, 1.035, 'PET L10' if pet_kind == 'PET_L10' else 'PET L20', transform=ax.transAxes, ha='center', va='bottom', fontsize=6.5, fontweight='bold', clip_on=False)
    return (image, panel_label(ax, panel))

def build_figure() -> None:
    apply_style()
    ce, mapping, contact = load_data()
    fig = plt.figure(figsize=(7.0, 5.25))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.12], left=0.09, right=0.93, top=0.955, bottom=0.105, hspace=0.4)
    top = outer[0].subgridspec(1, 2, width_ratios=[0.95, 2.25], wspace=0.27)
    bottom = outer[1].subgridspec(1, 3, width_ratios=[1.0, 1.0, 0.045], wspace=0.2)
    ax_ce = fig.add_subplot(top[0, 0])
    ax_mapping = fig.add_subplot(top[0, 1])
    ax_l10 = fig.add_subplot(bottom[0, 0])
    ax_l20 = fig.add_subplot(bottom[0, 1])
    ax_colorbar = fig.add_subplot(bottom[0, 2])
    labels = [(ax_ce, draw_ce_panel(ax_ce, ce)), (ax_mapping, draw_mapping_panel(ax_mapping, mapping))]
    image, l10_label = draw_contact_panel(ax_l10, contact, mapping, 'PET_L10', 'c', True)
    _, l20_label = draw_contact_panel(ax_l20, contact, mapping, 'PET_L20', 'd', False)
    labels.extend([(ax_l10, l10_label), (ax_l20, l20_label)])
    colorbar = fig.colorbar(image, cax=ax_colorbar, ticks=[0.0, 0.5, 1.0])
    colorbar.ax.set_yticklabels(['0.0', '0.5', '1.0'])
    colorbar.set_label('Mean contact fraction', labelpad=3)
    colorbar.outline.set_linewidth(0.5)
    figure_width_pt = fig.get_figwidth() * 72.0
    l10_box = ax_l10.get_position()
    l20_box = ax_l20.get_position()
    bar_box = ax_colorbar.get_position()
    l20_left = l10_box.x1 + 18.0 / figure_width_pt
    ax_l20.set_position([l20_left, l20_box.y0, l20_box.width, l20_box.height])
    ax_colorbar.set_position([l20_left + l20_box.width + 9.0 / figure_width_pt, bar_box.y0, bar_box.width, bar_box.height])
    for ax, label in labels:
        align_panel_label_to_yticks(fig, ax, label)
    return fig
CE_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/hic_pairwise_ce_alignment_summary.csv'
MAPPING_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/hic_structural_equivalent_mapping.csv'
CONTACT_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/hic_equivalent_contact_fraction_by_length.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S07', build_figure)
