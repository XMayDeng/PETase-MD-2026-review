from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
PROTEIN_ORDER = ['Pro00083', 'Pro00057', 'Pro00062', 'Pro00075', 'Pro00121', 'Pro00137']
PROTEIN_NAMES = {'Pro00083': 'IsPETase', 'Pro00057': 'TfCut1', 'Pro00062': 'LCC', 'Pro00075': 'FoCut5a', 'Pro00121': 'HiC', 'Pro00137': 'PHL7'}
CLASS_ORDER = ['Catalytic/oxyanion', 'W-loop/flexible or rigid pair', 'Cleft/rim patch', 'Flap/binding loop', 'Subsite/hotspot/stability', 'Other/candidate']
CLASS_LABELS = {'Catalytic/oxyanion': 'Catalytic/oxyanion', 'W-loop/flexible or rigid pair': 'W-loop/rigid pair', 'Cleft/rim patch': 'Cleft/rim patch', 'Flap/binding loop': 'Flap/binding loop', 'Subsite/hotspot/stability': 'Subsite/hotspot', 'Other/candidate': 'Other/candidate'}
CLASS_COLORS = {'Catalytic/oxyanion': '#E07B54', 'W-loop/flexible or rigid pair': '#E8C547', 'Cleft/rim patch': '#5B9BD5', 'Flap/binding loop': '#D98EC7', 'Subsite/hotspot/stability': '#70AD47', 'Other/candidate': '#8E939A'}
STABILITY_WARNING_PROTEINS = {'Pro00083', 'Pro00062', 'Pro00075'}
SMALL_N_PROTEIN = 'Pro00137'
INK = '#20242B'
BACKGROUND = '#FFFFFF'
ZERO_LINE = '#737B85'

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.2, 'axes.labelsize': 6.8, 'xtick.labelsize': 6.0, 'ytick.labelsize': 6.0, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

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

def load_data() -> pd.DataFrame:
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    data = pd.read_csv(SOURCE)
    if len(data) != 36:
        raise ValueError(f'Expected 36 class comparisons, found {len(data)}')
    observed_index = set(zip(data['protein_id'], data['landmark_class'], strict=True))
    expected_index = {(protein_id, landmark_class) for protein_id in PROTEIN_ORDER for landmark_class in CLASS_ORDER}
    if observed_index != expected_index:
        raise ValueError('Class-comparison grid does not match the registered panel')
    if not (data['bootstrap_ci95_low_pct_points'] <= data['delta_pct_points']).all() or not (data['delta_pct_points'] <= data['bootstrap_ci95_high_pct_points']).all():
        raise ValueError('A bootstrap interval does not contain its point estimate')
    if data['figure7_match_error_pct_points'].abs().max() > 1e-10:
        raise ValueError('Registered main-figure point estimates are not reproduced')
    expected_counts = {'Pro00083': (4, 7), 'Pro00057': (6, 6), 'Pro00062': (3, 8), 'Pro00075': (5, 7), 'Pro00121': (5, 8), 'Pro00137': (3, 2)}
    for protein_id, (n_l10, n_l20) in expected_counts.items():
        group = data[data['protein_id'] == protein_id]
        if not group['n_l10'].eq(n_l10).all() or not group['n_l20'].eq(n_l20).all():
            raise ValueError(f'Unexpected retained denominator for {protein_id}')
        if not group['small_n_flag'].eq(protein_id == SMALL_N_PROTEIN).all():
            raise ValueError(f'Unexpected small-n flag for {protein_id}')
        if not group['stability_warning'].eq(protein_id in STABILITY_WARNING_PROTEINS).all():
            raise ValueError(f'Unexpected stability-WARN flag for {protein_id}')
    return data

def draw_panel(ax: plt.Axes, data: pd.DataFrame, protein_id: str, panel_index: int) -> mpl.text.Text:
    group = data[data['protein_id'] == protein_id].set_index('landmark_class').reindex(CLASS_ORDER)
    y_positions = np.arange(len(CLASS_ORDER))
    for y, landmark_class in enumerate(CLASS_ORDER):
        row = group.loc[landmark_class]
        value = float(row['delta_pct_points'])
        low = float(row['bootstrap_ci95_low_pct_points'])
        high = float(row['bootstrap_ci95_high_pct_points'])
        color = CLASS_COLORS[landmark_class]
        ax.errorbar(value, y, xerr=np.array([[value - low], [high - value]]), fmt='o', markersize=4.1, markerfacecolor=BACKGROUND if bool(row['small_n_flag']) else color, markeredgecolor=INK, markeredgewidth=0.55, ecolor=color, elinewidth=1.0, capsize=2.0, capthick=0.8, zorder=3)
    ax.axvline(0, color=ZERO_LINE, linewidth=0.8, linestyle='--', zorder=1)
    ax.set_xlim(-60, 60)
    ax.set_xticks([-50, -25, 0, 25, 50])
    ax.set_yticks(y_positions)
    ax.set_yticklabels([CLASS_LABELS[value] for value in CLASS_ORDER])
    ax.invert_yaxis()
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='x', labelbottom=panel_index >= 4)
    n_l10 = int(group['n_l10'].iloc[0])
    n_l20 = int(group['n_l20'].iloc[0])
    symbol = ''
    if protein_id == SMALL_N_PROTEIN:
        symbol += '†'
    ax.text(0.5, 1.035, f'{PROTEIN_NAMES[protein_id]}{symbol}  (n={n_l10}/{n_l20})', transform=ax.transAxes, ha='center', va='bottom', fontsize=6.8, fontweight='bold', clip_on=False)
    return panel_label(ax, chr(ord('a') + panel_index))

def build_figure() -> None:
    apply_style()
    data = load_data()
    fig, axes = plt.subplots(3, 2, figsize=(7.0, 5.45), sharex=True, gridspec_kw={'left': 0.2, 'right': 0.985, 'top': 0.955, 'bottom': 0.105, 'wspace': 0.72, 'hspace': 0.34})
    labels = []
    for panel_index, (ax, protein_id) in enumerate(zip(axes.flat, PROTEIN_ORDER, strict=True)):
        labels.append((ax, draw_panel(ax, data, protein_id, panel_index)))
    fig.supxlabel('L20 − L10 contact-class share (percentage points)', x=0.59, y=0.025, fontsize=6.8)
    for ax, label in labels:
        align_panel_label_to_yticks(fig, ax, label)
    fig._supxlabel._publication_role = 'normal'
    return fig
SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/class_contact_shift_bootstrap.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S03', build_figure)
