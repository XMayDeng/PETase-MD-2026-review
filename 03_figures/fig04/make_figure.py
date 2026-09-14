from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from revised_figure_layout import align_panel_labels_to_yticklabels, contrasting_text_color, panel_label
PROTEIN_ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']
CLASS_ORDER = ['Catalytic/oxyanion', 'W-loop/flexible or rigid pair', 'Cleft/rim patch', 'Flap/binding loop', 'Subsite/hotspot/stability', 'Other/candidate']
CLASS_DISPLAY = {'Catalytic/oxyanion': 'Catalytic/oxyanion', 'W-loop/flexible or rigid pair': 'W-loop/rigid pair', 'Cleft/rim patch': 'Cleft/rim patch', 'Flap/binding loop': 'Flap/binding loop', 'Subsite/hotspot/stability': 'Subsite/hotspot', 'Other/candidate': 'Other/candidate'}
CLASS_TICK_DISPLAY = {'Catalytic/oxyanion': 'Catalytic/\noxyanion', 'W-loop/flexible or rigid pair': 'W-loop/\nrigid pair', 'Cleft/rim patch': 'Cleft/rim\npatch', 'Flap/binding loop': 'Flap/\nbinding loop', 'Subsite/hotspot/stability': 'Subsite/\nhotspot', 'Other/candidate': 'Other/\ncandidate'}
CLASS_COLORS = {'Catalytic/oxyanion': '#E07B54', 'W-loop/flexible or rigid pair': '#E8C547', 'Cleft/rim patch': '#5B9BD5', 'Flap/binding loop': '#D98EC7', 'Subsite/hotspot/stability': '#70AD47', 'Other/candidate': '#A6A6A6'}
CLASS_HATCHES = {'Catalytic/oxyanion': '////', 'W-loop/flexible or rigid pair': '....', 'Cleft/rim patch': '\\\\', 'Flap/binding loop': 'xxxx', 'Subsite/hotspot/stability': '++++', 'Other/candidate': ''}
PROTEIN_COLORS = {'TfCut1': '#C96A25', 'HiC': '#4F81BD'}
INK = '#20242B'
MUTED = '#626A75'
BACKGROUND = '#FFFFFF'
from data_checks import require_complete

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'axes.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'hatch.linewidth': 0.35})

def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    paths = [SHARE_PATH, COVERAGE_PATH, CLASS_BOOTSTRAP_PATH, CANDIDATE_PATH]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    shares = pd.read_csv(SHARE_PATH)
    coverage = pd.read_csv(COVERAGE_PATH)
    class_bootstrap = pd.read_csv(CLASS_BOOTSTRAP_PATH)
    candidates = pd.read_csv(CANDIDATE_PATH)
    for frame in (shares, coverage, class_bootstrap, candidates):
        frame['protein_short'] = frame['protein_short'].replace({'CUT1': 'TfCut1'})
    main_shares = shares[shares['pet_kind'].isin(['PET_L10', 'PET_L20'])].copy()
    sums = main_shares.groupby(['protein_short', 'pet_kind'])['contact_mass_share'].sum()
    if len(sums) != 12 or not np.allclose(sums.to_numpy(float), 1.0, atol=1e-09):
        raise ValueError('Each enzyme-length contact-mass composition must sum to one')
    main_coverage = coverage[coverage['pet_kind'].isin(['PET_L10', 'PET_L20'])].copy()
    if len(main_coverage) != 12:
        raise ValueError('Expected exactly 12 L10/L20 coverage rows')
    expected_counts = {('IsPETase', 'PET_L10'): 4, ('IsPETase', 'PET_L20'): 7, ('TfCut1', 'PET_L10'): 6, ('TfCut1', 'PET_L20'): 6, ('LCC', 'PET_L10'): 3, ('LCC', 'PET_L20'): 8, ('FoCut5a', 'PET_L10'): 5, ('FoCut5a', 'PET_L20'): 7, ('HiC', 'PET_L10'): 5, ('HiC', 'PET_L20'): 8, ('PHL7', 'PET_L10'): 3, ('PHL7', 'PET_L20'): 2}
    observed_counts = {(row.protein_short, row.pet_kind): int(row.n_retained_fingerprint_eligible) for row in main_coverage.itertuples(index=False)}
    if observed_counts != expected_counts:
        raise ValueError(f'Unexpected retained denominators: {observed_counts}')
    class_rows = class_bootstrap[class_bootstrap['landmark_class'].isin(CLASS_ORDER)].copy()
    if len(class_rows) != 36:
        raise ValueError(f'Expected 36 enzyme-class bootstrap rows, found {len(class_rows)}')
    candidate_rows = candidates[candidates['protein_short'].isin(['TfCut1', 'HiC']) & (candidates['shift_class'] == 'L20_enriched_supported_candidate')].copy()
    expected_residues = {('TfCut1', 'T62'), ('TfCut1', 'T64'), ('TfCut1', 'S67'), ('TfCut1', 'F210'), ('TfCut1', 'N213'), ('HiC', 'L138'), ('HiC', 'T164'), ('HiC', 'T166')}
    observed_residues = set(zip(candidate_rows['protein_short'], candidate_rows['residue_label']))
    if observed_residues != expected_residues:
        raise ValueError(f'Unexpected supported candidate set: {observed_residues}')
    return (main_shares, main_coverage, class_rows, candidate_rows)

def draw_panel_a(ax: plt.Axes, shares: pd.DataFrame, coverage: pd.DataFrame) -> mpl.text.Text:
    pair_offset = 0.19
    bar_height = 0.3
    group_y = np.arange(len(PROTEIN_ORDER), dtype=float)
    length_specs = [('PET_L10', 'L10', -pair_offset), ('PET_L20', 'L20', pair_offset)]
    coverage_lookup = {(row.protein_short, row.pet_kind): int(row.n_retained_fingerprint_eligible) for row in coverage.itertuples(index=False)}
    for protein_index, protein in enumerate(PROTEIN_ORDER):
        for pet_kind, pet_label, offset in length_specs:
            y = group_y[protein_index] + offset
            left = 0.0
            for cls in CLASS_ORDER:
                match = shares[(shares['protein_short'] == protein) & (shares['pet_kind'] == pet_kind) & (shares['landmark_class'] == cls)]
                value = 100.0 * float(match['contact_mass_share'].iloc[0]) if len(match) else 0.0
                bar = ax.barh(y, value, left=left, height=bar_height, color=CLASS_COLORS[cls], edgecolor='#4A4F57', linewidth=0.35)[0]
                bar.set_hatch(CLASS_HATCHES[cls])
                if value >= 14.0:
                    text_color = contrasting_text_color(bar.get_facecolor(), dark=INK)
                    label = ax.text(left + value / 2.0, y, f'{value:.0f}', ha='center', va='center', fontsize=6.0, fontweight='bold', color=text_color)
                    if bar.get_hatch():
                        label.set_path_effects([mpatheffects.withStroke(linewidth=0.9, foreground=bar.get_facecolor())])
                left += value
            n = coverage_lookup[protein, pet_kind]
            small_mark = '†' if protein == 'PHL7' else ''
            ax.text(-0.018, y, f'{pet_label}  n={n}{small_mark}', transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=5.8, color=INK, clip_on=False)
    for separator in np.arange(0.5, len(PROTEIN_ORDER) - 0.25, 1.0):
        ax.axhline(separator, color='#E4E7EB', linewidth=0.45, zorder=0)
    ax.set_yticks(group_y)
    ax.set_yticklabels(PROTEIN_ORDER, fontweight='bold')
    ax.tick_params(axis='y', length=0, pad=46)
    ax.invert_yaxis()
    ax.set_xlim(0.0, 100.0)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel('Retained-interface contact-mass share (%)')
    ax.spines[['left', 'top', 'right']].set_visible(False)
    return panel_label(ax, 'a')

def draw_panel_b(ax: plt.Axes, class_rows: pd.DataFrame) -> mpl.text.Text:
    matrix = require_complete(class_rows.pivot(index='protein_short', columns='landmark_class', values='delta_pct_points').reindex(index=PROTEIN_ORDER, columns=CLASS_ORDER))
    vmax = max(22.0, float(np.ceil(np.abs(matrix.to_numpy(float)).max())))
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = mcolors.LinearSegmentedColormap.from_list('loss_gain', ['#C96A4B', '#F7F7F7', '#4F81BD'])
    image = ax.imshow(matrix.to_numpy(float), cmap=cmap, norm=norm, aspect='auto')
    indexed = class_rows.set_index(['protein_short', 'landmark_class'])
    for i, protein in enumerate(PROTEIN_ORDER):
        for j, cls in enumerate(CLASS_ORDER):
            row = indexed.loc[protein, cls]
            value = float(row['delta_pct_points'])
            label = '0.00' if round(value, 2) == 0 else f'{value:+.2f}'
            robust_interval = bool(row['ci95_excludes_zero']) and (not bool(row['small_n_flag']))
            if robust_interval:
                label = f'{label}•'
                ax.add_patch(mpatches.Rectangle((j - 0.48, i - 0.48), 0.96, 0.96, fill=False, edgecolor=INK, linewidth=0.9))
            text_color = contrasting_text_color(image.cmap(image.norm(value)), dark=INK)
            ax.text(j, i, label, ha='center', va='center', fontsize=6.0, fontweight='normal', color=text_color)
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_TICK_DISPLAY[cls] for cls in CLASS_ORDER], rotation=38, ha='right', rotation_mode='anchor', linespacing=0.92)
    ax.set_yticks(range(len(PROTEIN_ORDER)))
    ax.set_yticklabels(['PHL7†' if value == 'PHL7' else value for value in PROTEIN_ORDER])
    ax.tick_params(length=0, pad=1.5)
    ax.set_xticks(np.arange(-0.5, len(CLASS_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(PROTEIN_ORDER), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.65)
    ax.tick_params(which='minor', bottom=False, left=False)
    header = panel_label(ax, 'b')
    cbar = ax.figure.colorbar(image, ax=ax, orientation='vertical', fraction=0.045, pad=0.07, aspect=24, ticks=[-vmax, 0, vmax])
    cbar.ax.yaxis.set_label_position('right')
    cbar.set_label('Percentage-point change', labelpad=1)
    cbar.outline.set_linewidth(0.5)
    return header

def draw_panel_c(ax: plt.Axes, candidates: pd.DataFrame) -> mpl.text.Text:
    candidate_order = [('TfCut1', 'F210'), ('TfCut1', 'S67'), ('TfCut1', 'T62'), ('TfCut1', 'N213'), ('TfCut1', 'T64'), ('HiC', 'T164'), ('HiC', 'L138'), ('HiC', 'T166')]
    chemistry_follow_up = {('TfCut1', 'F210'), ('TfCut1', 'S67'), ('TfCut1', 'T64'), ('HiC', 'L138'), ('HiC', 'T166')}
    indexed = candidates.set_index(['protein_short', 'residue_label'])
    y = np.arange(len(candidate_order), dtype=float)
    for yi, key in zip(y, candidate_order):
        row = indexed.loc[key]
        value = float(row['mean_delta_l20_minus_l10'])
        low = float(row['delta_ci95_low'])
        high = float(row['delta_ci95_high'])
        protein, _ = key
        color = PROTEIN_COLORS[protein]
        marker = 's' if protein == 'TfCut1' else 'o'
        filled = key in chemistry_follow_up
        ax.errorbar(value, yi, xerr=np.array([[value - low], [high - value]]), fmt=marker, color=color, markerfacecolor=color if filled else 'white', markeredgecolor=color, markeredgewidth=0.9, markersize=4.6, elinewidth=0.9, capsize=1.8, zorder=3)
    ax.axvline(0, color='#59616B', linewidth=0.7)
    ax.axhline(4.5, color='#D5DAE0', linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([f'{protein} {residue}' for protein, residue in candidate_order])
    ax.invert_yaxis()
    ax.set_xlim(-0.05, 0.86)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
    ax.set_xlabel('L20 − L10 mean contact fraction')
    return panel_label(ax, 'c')

def build_figure() -> None:
    apply_style()
    shares, coverage, class_rows, candidates = load_inputs()
    fig = plt.figure(figsize=(7.0, 5.8))
    grid = fig.add_gridspec(3, 2, height_ratios=[1.4, 0.25, 1.0], width_ratios=[1.1, 0.9], left=0.16, right=0.985, top=0.965, bottom=0.14, hspace=0.5, wspace=0.62)
    ax_a = fig.add_subplot(grid[0, :])
    ax_legend = fig.add_subplot(grid[1, :])
    ax_b = fig.add_subplot(grid[2, 0])
    ax_c = fig.add_subplot(grid[2, 1])
    panel_headers = [(ax_a, draw_panel_a(ax_a, shares, coverage)), (ax_b, draw_panel_b(ax_b, class_rows)), (ax_c, draw_panel_c(ax_c, candidates))]
    composition_box = ax_a.get_position()
    ax_a.set_position([0.21, composition_box.y0, composition_box.x1 - 0.21, composition_box.height])
    heatmap_box = ax_b.get_position()
    ax_b.set_position([0.08, heatmap_box.y0, heatmap_box.x1 - 0.08, heatmap_box.height])
    legend_handles = [mpatches.Patch(facecolor=CLASS_COLORS[cls], edgecolor='#4A4F57', linewidth=0.4, hatch=CLASS_HATCHES[cls], label=CLASS_DISPLAY[cls]) for cls in CLASS_ORDER]
    ax_legend.legend(handles=legend_handles, loc='center', bbox_to_anchor=(0.0, 0.0, 1.0, 1.0), bbox_transform=ax_legend.transAxes, mode='expand', ncol=3, frameon=False, fontsize=6.5, columnspacing=0.7, handlelength=0.95, handletextpad=0.35, handleheight=0.9, borderaxespad=0.0)
    ax_legend.set_axis_off()
    legend_box = ax_legend.get_position()
    ax_legend.set_position([legend_box.x0, legend_box.y0 - 0.018, legend_box.width, legend_box.height])
    fig._publication_panel_legend_layout = {'key': 'Figure_04', 'panel': ax_a, 'lower': [ax_b, ax_c], 'legend': ax_legend.get_legend(), 'legend_gap_pt': 10.0, 'lower_row': fig.axes[2:], 'lower_row_gap_pt': 12.0}
    align_panel_labels_to_yticklabels(fig, panel_headers)
    return fig
SHARE_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/pet_length_landmark_share_filtered_contact_ge_0_2.csv'
COVERAGE_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/pet_length_case_coverage.csv'
CLASS_BOOTSTRAP_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/class_contact_shift_bootstrap.csv'
CANDIDATE_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/manuscript_candidate_residues.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_04', build_figure)
