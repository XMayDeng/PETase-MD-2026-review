from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import pandas as pd
ENZYME_ORDER = ['TfCut1', 'LCC']
COLORS = {'TfCut1': '#ED7D31', 'LCC': '#9DCB84'}
LINE_STYLES = {1: '-', 2: '--', 3: ':'}
MARKERS = {1: 'o', 2: 's', 3: '^'}
INK = '#20242B'
BACKGROUND = '#FFFFFF'

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

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

def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    for source in (FRAME_SOURCE, TRAJECTORY_SOURCE, COMPARISON_SOURCE):
        if not source.is_file():
            raise FileNotFoundError(source)
    frames = pd.read_csv(FRAME_SOURCE)
    trajectories = pd.read_csv(TRAJECTORY_SOURCE)
    comparison_table = pd.read_csv(COMPARISON_SOURCE)
    if len(comparison_table) != 1:
        raise ValueError('Expected one registered LCC-minus-TfCut1 comparison')
    comparison = comparison_table.iloc[0]
    if len(frames) != 48006:
        raise ValueError(f'Expected 48,006 frame rows, found {len(frames)}')
    counts = frames.groupby(['enzyme', 'replica_id']).size()
    expected_index = pd.MultiIndex.from_product([ENZYME_ORDER, [1, 2, 3]], names=['enzyme', 'replica_id'])
    if not counts.reindex(expected_index).eq(8001).all():
        raise ValueError('Expected 8,001 frames for each enzyme/replica trajectory')
    if len(trajectories) != 6:
        raise ValueError(f'Expected six trajectory summaries, found {len(trajectories)}')
    if set(frames['enzyme']) != set(ENZYME_ORDER):
        raise ValueError('Unexpected enzyme in Figure S4 frame source')
    if not frames['pet_kind'].eq('PET_L20').all():
        raise ValueError('Figure S4 must use PET L20 trajectories only')
    if not frames['direction'].eq('head-side').all():
        raise ValueError('Figure S4 must use head-side trajectories only')
    return (frames, trajectories, comparison)

def build_figure() -> None:
    apply_style()
    frames, trajectories, _ = load_data()
    fig, (ax_time, ax_dist) = plt.subplots(1, 2, figsize=(7.0, 3.25), sharey=True, gridspec_kw={'width_ratios': [1.82, 1.0], 'left': 0.09, 'right': 0.985, 'top': 0.86, 'bottom': 0.17, 'wspace': 0.25})
    for (enzyme, replica), group in frames.groupby(['enzyme', 'replica_id'], sort=False):
        group = group.sort_values('time_ns')
        smoothed = group['pet_rg_nm'].rolling(101, center=True, min_periods=1).mean()
        ax_time.plot(group['time_ns'], group['pet_rg_nm'], color=COLORS[enzyme], alpha=0.09, linewidth=0.35)
        ax_time.plot(group['time_ns'], smoothed, color=COLORS[enzyme], linestyle=LINE_STYLES[int(replica)], linewidth=1.1)
    ax_time.set_xlim(20, 100)
    ax_time.set_ylim(0.83, 2.22)
    ax_time.set_xlabel('Time (ns)')
    ax_time.set_ylabel('PET radius of gyration, $R_g$ (nm)')
    ax_time.grid(False)
    ax_time.spines[['top', 'right']].set_visible(False)
    left_label = panel_label(ax_time, 'a')
    enzyme_handles = [mlines.Line2D([], [], color=COLORS[enzyme], linewidth=1.6, label=enzyme) for enzyme in ENZYME_ORDER]
    replica_handles = [mlines.Line2D([], [], color='#545B65', linestyle=LINE_STYLES[replica], linewidth=1.2, marker=MARKERS[replica], markersize=4.2, markerfacecolor=BACKGROUND, markeredgewidth=0.8, label=f'r{replica}') for replica in (1, 2, 3)]
    ax_time.legend(handles=enzyme_handles, loc='lower left', bbox_to_anchor=(0.06, 1.01), frameon=False, ncol=2, handlelength=1.35, columnspacing=1.1, borderaxespad=0)
    positions = {'TfCut1': 0, 'LCC': 1}
    frame_values = [frames.loc[frames['enzyme'] == enzyme, 'pet_rg_nm'].to_numpy(float) for enzyme in ENZYME_ORDER]
    violins = ax_dist.violinplot(frame_values, positions=[positions[enzyme] for enzyme in ENZYME_ORDER], widths=0.72, showextrema=False)
    for body, enzyme in zip(violins['bodies'], ENZYME_ORDER, strict=True):
        body.set_facecolor(COLORS[enzyme])
        body.set_edgecolor('#4A4F57')
        body.set_alpha(0.24)
        body.set_linewidth(0.5)
    for enzyme, xpos in positions.items():
        group = trajectories[trajectories['enzyme'] == enzyme].sort_values('replica_id')
        if len(group) != 3:
            raise ValueError(f'Expected three trajectory summaries for {enzyme}')
        for row in group.itertuples(index=False):
            replica = int(row.replica_id)
            ax_dist.scatter(xpos + (replica - 2) * 0.08, row.mean_rg_nm, marker=MARKERS[replica], s=26, facecolor='white', edgecolor=COLORS[enzyme], linewidth=1.0, zorder=4)
        values = group['mean_rg_nm'].to_numpy(float)
        ax_dist.plot([xpos - 0.2, xpos + 0.2], [values.mean(), values.mean()], color=INK, linewidth=1.25, zorder=5)
    ax_dist.set_xticks([0, 1], ENZYME_ORDER)
    ax_dist.tick_params(labelleft=True)
    ax_dist.grid(False)
    ax_dist.spines[['top', 'right']].set_visible(False)
    right_label = panel_label(ax_dist, 'b')
    ax_dist.legend(handles=replica_handles, loc='lower center', bbox_to_anchor=(0.5, 1.01), frameon=False, ncol=3, handlelength=1.4, columnspacing=0.8, borderaxespad=0)
    align_panel_label_to_yticks(fig, ax_time, left_label)
    align_panel_label_to_yticks(fig, ax_dist, right_label)
    return fig
FRAME_SOURCE = ROOT / '02_postprocessing/results/frame_observables/pet_rg_per_frame.csv.gz'
TRAJECTORY_SOURCE = ROOT / '02_postprocessing/results/trajectory_summaries/pet_rg_trajectory_summary.csv'
COMPARISON_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/pet_rg_group_comparison.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S06', build_figure)
