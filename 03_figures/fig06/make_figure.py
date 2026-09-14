from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
from PIL import Image
from revised_figure_layout import panel_label, align_panel_labels_to_yticklabels
ROOT = Path(__file__).resolve().parents[2]
HBOND = '#2589A8'
AROMATIC = '#8965A8'
NONPOLAR = '#737B84'
NONPOLAR_LIGHT = '#A6ADB5'
INK = '#20242B'
MUTED = '#626A75'
BACKGROUND = '#FFFFFF'
SNAPSHOT_CROPS = {'T64': (580, 130, 1820, 1250), 'S67': (280, 390, 1520, 1510), 'T166': (550, 360, 1790, 1480)}

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def load_chemistry_effects() -> pd.DataFrame:
    hbond = pd.read_csv(HBOND_PATH)
    local = pd.read_csv(LOCAL_PATH)
    requested_local = [('T64_any_native_ADH_le_30', 'TfCut1 T64 · any direct H bond', 'H bond', False), ('T64_backbone_native_ADH_le_30', 'TfCut1 T64 · backbone N', 'H bond', False), ('T64_sidechain_native_ADH_le_30', 'TfCut1 T64 · side-chain OG1', 'H bond', False), ('F210_aromatic_proximity_d0p55_any', 'TfCut1 F210 · aromatic proximity', 'Aromatic', False), ('F210_any_pi_d0p55_a30_o0p20_any', 'TfCut1 F210 · pi-like geometry', 'Aromatic', False), ('L138_nonpolar_d0p40_any', 'HiC L138 · nonpolar <0.40 nm', 'Nonpolar', False), ('L138_nonpolar_d0p45_any', 'HiC L138 · nonpolar <0.45 nm†', 'Nonpolar', True), ('L138_nonpolar_d0p50_any', 'HiC L138 · nonpolar <0.50 nm†', 'Nonpolar', True)]
    local_rows: dict[str, dict[str, object]] = {}
    for metric_id, label, chemistry, sensitivity in requested_local:
        selected = local[local['metric_id'] == metric_id]
        if len(selected) != 1:
            raise ValueError(f'Expected one local-interaction row for {metric_id}')
        row = selected.iloc[0]
        local_rows[metric_id] = {'label': label, 'chemistry': chemistry, 'sensitivity': sensitivity, 'estimate': float(row['delta_l20_minus_l10']), 'low': float(row['delta_ci95_low']), 'high': float(row['delta_ci95_high']), 'resolved': bool(row['ci_excludes_zero']), 'n_l10': int(row['n_l10']), 'n_l20': int(row['n_l20'])}

    def hbond_row(protein: str, residue: int, sidechain: bool) -> dict[str, object]:
        selected = hbond[(hbond['protein_name'] == protein) & (hbond['residue_number'] == residue)]
        if len(selected) != 1:
            raise ValueError(f'Expected one H-bond row for {protein} {residue}')
        row = selected.iloc[0]
        prefix = 'sidechain_' if sidechain else ''
        estimate_key = f'{prefix}delta_l20_minus_l10'
        low_key = f'{prefix}delta_ci95_low'
        high_key = f'{prefix}delta_ci95_high'
        resolved_key = f'{prefix}ci_excludes_zero'
        return {'chemistry': 'H bond', 'sensitivity': False, 'estimate': float(row[estimate_key]), 'low': float(row[low_key]), 'high': float(row[high_key]), 'resolved': bool(row[resolved_key]), 'n_l10': int(row['n_l10']), 'n_l20': int(row['n_l20'])}
    # HB01 is the registered primary direct-H-bond interval source. LOC01
    # retains the backbone-specific row and other chemistry-specific metrics.
    for metric_id, sidechain in (
            ('T64_any_native_ADH_le_30', False),
            ('T64_sidechain_native_ADH_le_30', True)):
        target = local_rows[metric_id]
        primary = hbond_row('CUT1', 64, sidechain=sidechain)
        if (target['n_l10'], target['n_l20'], target['resolved']) != (
                primary['n_l10'], primary['n_l20'], primary['resolved']) or not np.isclose(
                    target['estimate'], primary['estimate'], atol=1e-6, rtol=0):
            raise ValueError('HB01/LOC01 primary-estimate reconciliation failed: ' + metric_id)
        target['low'], target['high'] = primary['low'], primary['high']
    s67 = hbond_row('CUT1', 67, sidechain=False)
    s67['label'] = 'TfCut1 S67 · direct H bond'
    t166 = hbond_row('HiC', 166, sidechain=True)
    t166['label'] = 'HiC T166 · side-chain H bond'
    order = [local_rows['T64_any_native_ADH_le_30'], local_rows['T64_backbone_native_ADH_le_30'], local_rows['T64_sidechain_native_ADH_le_30'], s67, t166, local_rows['F210_aromatic_proximity_d0p55_any'], local_rows['F210_any_pi_d0p55_a30_o0p20_any'], local_rows['L138_nonpolar_d0p40_any'], local_rows['L138_nonpolar_d0p45_any'], local_rows['L138_nonpolar_d0p50_any']]
    result = pd.DataFrame(order)
    if len(result) != 10:
        raise ValueError('Figure 6 chemistry panel must contain ten rows')
    expected_resolved = [True, True, True, True, True, True, False, False, False, True]
    if result['resolved'].tolist() != expected_resolved:
        raise ValueError(f'Unexpected interval-resolution pattern: {result['resolved'].tolist()}')
    return result

def validate_frames() -> pd.DataFrame:
    registry = pd.read_csv(FRAME_REGISTRY_PATH)
    expected = {('CUT1', 64): ('PET_L20', 75490.0, 'OG1', 'O4'), ('CUT1', 67): ('PET_L20', 56190.0, 'OG', 'O2'), ('HiC', 166): ('PET_L20', 52440.0, 'OG1', 'O4')}
    selected_rows = []
    for (protein, residue), values in expected.items():
        selected = registry[(registry['protein_name'] == protein) & (registry['residue_number'] == residue)]
        if len(selected) != 1:
            raise ValueError(f'Expected one representative frame for {protein} {residue}')
        row = selected.iloc[0]
        observed = (row['pet_kind'], float(row['representative_time_ps']), row['protein_atom_name'], row['pet_atom_name'])
        if observed != values:
            raise ValueError(f'Unexpected representative frame metadata: {observed}')
        selected_rows.append(row)
    for path in SNAPSHOTS.values():
        if not path.is_file():
            raise FileNotFoundError(path)
        if Image.open(path).size != (2400, 1800):
            raise ValueError(f'Unexpected snapshot dimensions for {path}')
    return pd.DataFrame(selected_rows)

def draw_chemistry_panel(ax: plt.Axes, data: pd.DataFrame) -> mpl.text.Text:
    """Draw the chemistry-specific effect panel using the Figure 4c axis style."""
    y = np.arange(len(data))[::-1]
    ax.axvline(0, color='#9299A2', lw=0.8, ls=(0, (2.2, 2.2)), zorder=0)
    xlim = (-0.06, 0.68)
    for position, row in zip(y, data.itertuples(index=False), strict=True):
        color = {'H bond': HBOND, 'Aromatic': AROMATIC, 'Nonpolar': NONPOLAR_LIGHT if row.sensitivity else NONPOLAR}[row.chemistry]
        marker = '^' if row.sensitivity else {'H bond': 'o', 'Aromatic': 'D', 'Nonpolar': 's'}[row.chemistry]
        marker_face = color if row.resolved else BACKGROUND
        ax.errorbar(row.estimate, position, xerr=np.array([[row.estimate - row.low], [row.high - row.estimate]]), fmt=marker, markersize=4.6, markerfacecolor=marker_face, markeredgecolor=color, markeredgewidth=0.9, ecolor=color, elinewidth=1.05, capsize=2.0, capthick=0.8, zorder=3)
        text_x = min(row.high + 0.018, xlim[1] - 0.07)
        ax.text(text_x, position, f'{row.estimate:+.3f}', va='center', ha='left', fontsize=6.0, color=MUTED)
    ax.set_yticks(y, data['label'])
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.75, len(data) - 0.25)
    ax.set_xticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    ax.set_xlabel('L20 − L10 mean occupancy')
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='y', length=2.5, width=0.55, pad=3)
    for boundary in (6.5, 4.5, 2.5):
        ax.axhline(boundary, color='#E7E9ED', lw=0.7, zorder=0)
    chemistry_legend = [mlines.Line2D([], [], marker='o', linestyle='none', color=HBOND, label='Direct H bond'), mlines.Line2D([], [], marker='D', linestyle='none', color=AROMATIC, label='Aromatic'), mlines.Line2D([], [], marker='s', linestyle='none', color=NONPOLAR, label='Nonpolar'), mlines.Line2D([], [], marker='^', linestyle='none', color=NONPOLAR_LIGHT, label='Cutoff sensitivity')]
    ax.legend(handles=chemistry_legend, bbox_to_anchor=(0.995, 0.975), loc='upper right', frameon=False, ncol=2, handletextpad=0.25, columnspacing=0.8, labelspacing=0.35, borderaxespad=0.15)
    header = panel_label(ax, 'a')
    return header

def draw_snapshot_panel(ax: plt.Axes, path: Path, crop_box: tuple[int, int, int, int], letter: str, title: str, subtitle: str, distance_label: str, distance_xy: tuple[float, float]) -> None:
    with Image.open(path) as source:
        image = source.convert('RGB').crop(crop_box)
    ax.imshow(image, interpolation='lanczos')
    ax.set_axis_off()
    panel_label(ax, letter, x=0.01, y=1.015, va='bottom')
    ax.text(0.5, 1.015, title, transform=ax.transAxes, fontsize=7.5, fontweight='bold', va='bottom', ha='center')
    ax.text(0.5, -0.035, subtitle, transform=ax.transAxes, fontsize=5.8, color=MUTED, va='top', ha='center', clip_on=False)
    ax.text(distance_xy[0], distance_xy[1], distance_label, transform=ax.transAxes, fontsize=6.0, color=INK, va='center', ha='center', zorder=20)

def build_figure(rendered) -> None:
    global SNAPSHOTS
    SNAPSHOTS = {r: rendered / f'Figure_07_{r}_no_distance_label.png' for r in ['T64', 'S67', 'T166']}
    apply_style()
    chemistry = load_chemistry_effects()
    frames = validate_frames()
    fig = plt.figure(figsize=(7.0, 5.5))
    top = fig.add_gridspec(1, 1, left=0.315, right=0.985, top=0.955, bottom=0.515)
    bottom = fig.add_gridspec(1, 3, left=0.06, right=0.985, top=0.355, bottom=0.045, wspace=0.06)
    ax_a = fig.add_subplot(top[0, 0])
    ax_b = fig.add_subplot(bottom[0, 0])
    ax_c = fig.add_subplot(bottom[0, 1])
    ax_d = fig.add_subplot(bottom[0, 2])
    header_a = draw_chemistry_panel(ax_a, chemistry)
    align_panel_labels_to_yticklabels(fig, [(ax_a, header_a)])
    draw_snapshot_panel(ax_b, SNAPSHOTS['T64'], SNAPSHOT_CROPS['T64'], 'b', 'TfCut1 T64 (OG1)', 'L20, 75.49 ns', '2.86 Å', (0.3, 0.55))
    draw_snapshot_panel(ax_c, SNAPSHOTS['S67'], SNAPSHOT_CROPS['S67'], 'c', 'TfCut1 S67', 'L20, 56.19 ns', '2.59 Å', (0.35, 0.58))
    draw_snapshot_panel(ax_d, SNAPSHOTS['T166'], SNAPSHOT_CROPS['T166'], 'd', 'HiC T166', 'L20, 52.44 ns', '2.58 Å', (0.7, 0.67))
    return fig
HBOND_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/candidate_hbond_shift_uncertainty.csv'
LOCAL_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/interaction_shift_statistics.csv'
FRAME_REGISTRY_PATH = ROOT / '02_postprocessing/results/structural_snapshots/representative_hbond_frames.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_06', build_figure)
