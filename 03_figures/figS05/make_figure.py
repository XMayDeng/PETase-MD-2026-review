from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
PROTEIN_ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']
PET_ORDER = ['PET_L4', 'PET_L10', 'PET_L20']
PET_COLORS = {'PET_L4': '#2F4858', 'PET_L10': '#33658A', 'PET_L20': '#86BBD8'}
PET_MARKERS = {'PET_L4': 's', 'PET_L10': 'o', 'PET_L20': '^'}
PET_OFFSETS = {'PET_L4': -0.22, 'PET_L10': 0.0, 'PET_L20': 0.22}
INK = '#20242B'
BACKGROUND = '#FFFFFF'

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def panel_label(ax: plt.Axes, label: str) -> mpl.text.Text:
    return ax.text(0.0, 1.025, f'({label.lower()})', transform=ax.transAxes, ha='left', va='bottom', fontsize=8.5, fontweight='bold', clip_on=False, zorder=20)

def align_panel_label_to_yticks(fig: plt.Figure, ax: plt.Axes, label: mpl.text.Text) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    visible = [tick for tick in ax.get_yticklabels() if tick.get_visible() and tick.get_text().strip()]
    if not visible:
        return
    left_display = min((tick.get_window_extent(renderer=renderer).x0 for tick in visible))
    left_axes = ax.transAxes.inverted().transform((left_display, ax.bbox.y0))[0]
    label.set_x(left_axes)

def load_data() -> pd.DataFrame:
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    data = pd.read_csv(SOURCE)
    if len(data) != 71:
        raise ValueError(f'Expected 71 retained trajectories, found {len(data)}')
    if not data['retention_class'].eq('retained').all():
        raise ValueError('Figure S5 input must contain retained trajectories only')
    if set(data['protein']) != set(PROTEIN_ORDER):
        raise ValueError('Unexpected protein set in Figure S5 input')
    if not set(data['pet_kind']).issubset(PET_ORDER):
        raise ValueError('Unexpected PET length in Figure S5 input')
    expected_counts = {'IsPETase': {'PET_L4': 2, 'PET_L10': 4, 'PET_L20': 7}, 'TfCut1': {'PET_L4': 2, 'PET_L10': 6, 'PET_L20': 6}, 'LCC': {'PET_L4': 2, 'PET_L10': 3, 'PET_L20': 8}, 'FoCut5a': {'PET_L4': 0, 'PET_L10': 5, 'PET_L20': 7}, 'HiC': {'PET_L4': 0, 'PET_L10': 5, 'PET_L20': 8}, 'PHL7': {'PET_L4': 1, 'PET_L10': 3, 'PET_L20': 2}}
    observed = data.groupby(['protein', 'pet_kind']).size().to_dict()
    expected_observed = {(protein, pet_kind): count for protein, by_length in expected_counts.items() for pet_kind, count in by_length.items() if count > 0}
    if observed != expected_observed:
        raise ValueError(f'Unexpected retained-count map: {observed}')
    return data

def draw_panel(ax: plt.Axes, data: pd.DataFrame, value_column: str, title: str, letter: str) -> mpl.text.Text:
    x = np.arange(len(PROTEIN_ORDER), dtype=float)
    for pet_kind in PET_ORDER:
        for protein_index, protein in enumerate(PROTEIN_ORDER):
            values = 10.0 * data.loc[(data['protein'] == protein) & (data['pet_kind'] == pet_kind), value_column].to_numpy(float)
            position = protein_index + PET_OFFSETS[pet_kind]
            if len(values) >= 3 and np.ptp(values) > 1e-08:
                violin = ax.violinplot([values], positions=[position], widths=0.19, showextrema=False, bw_method=0.5)
                for body in violin['bodies']:
                    body.set_facecolor(PET_COLORS[pet_kind])
                    body.set_edgecolor('#4A4F57')
                    body.set_linewidth(0.45)
                    body.set_alpha(0.28)
            if len(values):
                jitter = np.linspace(-0.035, 0.035, len(values)) if len(values) > 1 else np.array([0.0])
                ax.scatter(position + jitter, values, marker=PET_MARKERS[pet_kind], s=17, facecolor='white', edgecolor=PET_COLORS[pet_kind], linewidth=0.9, zorder=3)
                mean = float(values.mean())
                ax.plot([position - 0.055, position + 0.055], [mean, mean], color=INK, linewidth=1.05, zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(PROTEIN_ORDER, rotation=25, ha='right')
    ax.set_ylim(1.7, 40.5)
    ax.set_title(title, loc='left', x=0.03, pad=5, fontweight='bold')
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    return panel_label(ax, letter)

def build_figure() -> None:
    apply_style()
    data = load_data()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.35), sharey=True, gridspec_kw={'left': 0.09, 'right': 0.985, 'top': 0.9, 'bottom': 0.235, 'wspace': 0.2})
    left_label = draw_panel(axes[0], data, 'mean_ser_og_to_any_pet_heavy_min_nm', 'Any PET heavy atom', 'a')
    right_label = draw_panel(axes[1], data, 'mean_ser_og_to_ester_carbonyl_c_min_nm', 'PET ester carbonyl carbon', 'b')
    axes[0].set_ylabel('Catalytic Ser O$_{\\gamma}$–PET minimum distance (Å)')
    axes[1].tick_params(labelleft=False)
    handles = [mlines.Line2D([], [], color=PET_COLORS[pet_kind], marker=PET_MARKERS[pet_kind], markerfacecolor='white', markeredgewidth=0.9, linestyle='none', label=pet_kind.replace('PET_', '')) for pet_kind in PET_ORDER]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.025), columnspacing=1.35, handletextpad=0.4)
    align_panel_label_to_yticks(fig, axes[0], left_label)
    right_label.set_x(-0.1)
    return fig
SOURCE = ROOT / '02_postprocessing/results/trajectory_summaries/catalytic_ser_pet_trajectory_summary.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S05', build_figure)
