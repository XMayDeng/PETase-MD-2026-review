from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
ROOT = Path(__file__).resolve().parents[2]
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
LENGTH_ORDER = ('PET_L10', 'PET_L20')
COLORS = {'PET_L10': '#2684B7', 'PET_L20': '#D66A1F'}
MARKERS = {'PET_L10': 'o', 'PET_L20': '^'}
LINESTYLES = {'PET_L10': '-', 'PET_L20': '--'}
DISPLAY_NAMES = {'PET_L10': 'L10', 'PET_L20': 'L20'}
INK = '#20242B'
GRAY = '#737B85'
LIGHT_GRAY = '#D4D9DF'
BACKGROUND = '#FFFFFF'

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 7.6, 'axes.labelsize': 8.2, 'xtick.labelsize': 7.6, 'ytick.labelsize': 7.6, 'legend.fontsize': 7.6, 'axes.linewidth': 0.7, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def panel_label(ax: plt.Axes, label: str) -> mpl.text.Text:
    return ax.text(0.0, 1.035, f'({label.lower()})', transform=ax.transAxes, ha='left', va='bottom', fontsize=10.0, fontweight='bold', clip_on=False, zorder=20)

def align_panel_label_to_yticks(fig: plt.Figure, ax: plt.Axes, label: mpl.text.Text) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    visible = [tick for tick in ax.get_yticklabels() if tick.get_visible() and tick.get_text().strip()]
    if not visible:
        return
    left_display = min((tick.get_window_extent(renderer=renderer).x0 for tick in visible))
    left_axes = ax.transAxes.inverted().transform((left_display, ax.bbox.y0))[0]
    label.set_x(left_axes)

def load_data() -> tuple[pd.DataFrame, ...]:
    sources = [FRAME_SOURCE, TRAJECTORY_SOURCE, GROUP_SOURCE, REFERENCE_SOURCE, MIXTURE_SOURCE]
    for source in sources:
        if not source.is_file():
            raise FileNotFoundError(source)
    frames = pd.read_csv(FRAME_SOURCE, compression='gzip')
    trajectories = pd.read_csv(TRAJECTORY_SOURCE)
    group_summary = pd.read_csv(GROUP_SOURCE)
    references = pd.read_csv(REFERENCE_SOURCE)
    mixtures = pd.read_csv(MIXTURE_SOURCE)
    if len(frames) != 40005 or frames['case'].nunique() != 5:
        raise ValueError('Unexpected PHL7 frame population')
    if not frames.groupby('case').size().eq(8001).all():
        raise ValueError('Every retained trajectory must contribute 8,001 frames')
    if set(frames['pet_kind']) != set(LENGTH_ORDER):
        raise ValueError('Unexpected PET-length labels')
    expected_trajectory_counts = {'PET_L10': 3, 'PET_L20': 2}
    observed_counts = trajectories.groupby('pet_kind').size().to_dict()
    if observed_counts != expected_trajectory_counts:
        raise ValueError(f'Unexpected trajectory counts: {observed_counts}')
    if len(trajectories) != 5 or len(references) != 8:
        raise ValueError('Unexpected trajectory or structural-reference count')
    selected = mixtures.loc[mixtures['selected_by_bic'].astype(bool)]
    if len(selected) != 15 or not selected['components'].eq(1).all():
        raise ValueError('The registered one-component sensitivity result changed')
    frames = frames.copy()
    frames['chi1_display_deg'] = -180.0 + frames['chi1_delta_from_trans_deg']
    if not np.allclose(np.mod(frames['chi1_display_deg'] + 180.0, 360.0) - 180.0, np.mod(frames['chi1_deg'] + 180.0, 360.0) - 180.0, atol=1e-08):
        raise ValueError('Periodic chi1 display mapping is not angle-equivalent')
    return (frames, trajectories, group_summary, references)

def highest_density_levels(density: np.ndarray) -> list[float]:
    values = density[np.isfinite(density) & (density > 0)]
    ordered = np.sort(values)[::-1]
    cumulative = np.cumsum(ordered) / ordered.sum()
    thresholds = []
    for mass in (0.95, 0.8, 0.5):
        index = min(int(np.searchsorted(cumulative, mass)), len(ordered) - 1)
        thresholds.append(float(ordered[index]))
    return sorted(set(thresholds))

def draw_density_contours(ax: plt.Axes, frame: pd.DataFrame, color: str, linestyle: str) -> None:
    x_edges = np.linspace(-205.0, -120.0, 171)
    y_edges = np.linspace(0.0, 135.0, 181)
    histogram, _, _ = np.histogram2d(frame['chi1_display_deg'].to_numpy(dtype=float), frame['chi2_deg'].to_numpy(dtype=float), bins=(x_edges, y_edges))
    density = gaussian_filter(histogram.T, sigma=2.0, mode='nearest')
    x_centers = (x_edges[:-1] + x_edges[1:]) / 2.0
    y_centers = (y_edges[:-1] + y_edges[1:]) / 2.0
    levels = highest_density_levels(density)
    ax.contour(x_centers, y_centers, density, levels=levels, colors=[color], linewidths=np.linspace(0.75, 1.15, len(levels)), linestyles=linestyle, alpha=0.96, zorder=2)

def draw_reference_space(ax: plt.Axes, frames: pd.DataFrame, references: pd.DataFrame) -> mpl.text.Text:
    for pet_kind in LENGTH_ORDER:
        draw_density_contours(ax, frames.loc[frames['pet_kind'].eq(pet_kind)], COLORS[pet_kind], LINESTYLES[pet_kind])
    ispetase = references.loc[references['reference'].eq('IsPETase_5XG0')]
    label_offsets = {'A': (8, 0), 'B': (-12, -18), 'C': (14, 17)}
    for row in ispetase.itertuples(index=False):
        ax.scatter(row.chi1_deg, row.chi2_deg, s=25, marker='s', facecolor=BACKGROUND, edgecolor=GRAY, linewidth=0.8, zorder=5)
        ax.annotate(str(row.chain), (row.chi1_deg, row.chi2_deg), xytext=label_offsets[str(row.chain)], textcoords='offset points', fontsize=8.0, fontweight='semibold', color=INK, ha='right' if row.chain == 'B' else 'left', va='center', arrowprops={'arrowstyle': '-', 'color': INK, 'lw': 0.6, 'shrinkB': 4} if row.chain != 'A' else None, zorder=10)
    phl7 = references.loc[references['reference'].isin(['PHL7_7NEI', 'PHL7_8BRA'])]
    ax.scatter(phl7['chi1_deg'], phl7['chi2_deg'], s=22, marker='D', facecolor=INK, edgecolor=BACKGROUND, linewidth=0.6, zorder=6)
    ax.annotate('PHL7 references', (float(phl7['chi1_deg'].mean()), float(phl7['chi2_deg'].mean())), xytext=(0.97, 0.59), textcoords='axes fraction', fontsize=7.6, color=INK, ha='right', va='center', arrowprops={'arrowstyle': '-', 'color': INK, 'lw': 0.6, 'shrinkB': 4}, zorder=10)
    ax.legend(handles=[Line2D([], [], color=COLORS[pet_kind], linestyle=LINESTYLES[pet_kind], linewidth=1.0, label=DISPLAY_NAMES[pet_kind]) for pet_kind in LENGTH_ORDER], loc='upper right', bbox_to_anchor=(0.97, 0.4), frameon=False, handlelength=2.0, handletextpad=0.6, borderaxespad=0)
    ax.set_xlim(-205, -120)
    ax.set_ylim(0, 135)
    ax.set_xticks([-200, -180, -160, -140, -120])
    ax.set_yticks([20, 50, 80, 110])
    ax.set_xlabel('W156 $\\chi_1$ (°)')
    ax.set_ylabel('W156 $\\chi_2$ (°)')
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    return panel_label(ax, 'a')

def short_case_label(case: str) -> str:
    text = case.replace('Pro00137_PET_', '')
    text = text.replace('_side', '').replace('bidirectional', 'bi')
    length, direction, replica = text.rsplit('_', 2)
    return f'{length} {direction} {replica}'

def draw_time_strips(axes: list[plt.Axes], frames: pd.DataFrame, references: pd.DataFrame) -> mpl.text.Text:
    cases = frames[['case', 'pet_kind']].drop_duplicates().sort_values(['pet_kind', 'case'])
    ispetase = references.loc[references['reference'].eq('IsPETase_5XG0')].set_index('chain')['chi2_deg'].to_dict()
    phl7 = references.loc[references['reference'].isin(['PHL7_7NEI', 'PHL7_8BRA'])]['chi2_deg']
    phl7_low = float(phl7.min())
    phl7_high = float(phl7.max())
    for index, (ax, row) in enumerate(zip(axes, cases.itertuples(index=False))):
        group = frames.loc[frames['case'].eq(row.case)].iloc[::10]
        ax.axhspan(phl7_low, phl7_high, color=LIGHT_GRAY, alpha=0.7, lw=0)
        for state in ('A', 'B', 'C'):
            ax.axhline(ispetase[state], color=GRAY, lw=0.6, ls=':' if state == 'A' else '--', alpha=0.9, zorder=0)
        ax.plot(group['time_ns'], group['chi2_deg'], color=COLORS[row.pet_kind], lw=0.6, alpha=0.92, rasterized=True)
        ax.text(0.01, 1.04, short_case_label(row.case), transform=ax.transAxes, fontsize=8.0, fontweight='semibold', color=INK, ha='left', va='bottom', clip_on=False, zorder=6)
        ax.set_xlim(20, 100)
        ax.set_ylim(0, 135)
        ax.set_yticks([20, 100])
        ax.tick_params(axis='y', labelleft=True, pad=2)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines['bottom'].set_visible(index == len(axes) - 1)
        ax.tick_params(axis='x', labelbottom=index == len(axes) - 1)
        if index != len(axes) - 1:
            ax.tick_params(axis='x', length=0)
    axes[-1].set_xlabel('Time (ns)', labelpad=2)
    axes[len(axes) // 2].set_ylabel('W156 $\\chi_2$ (°)', labelpad=5)
    return panel_label(axes[0], 'b')

def summary_lookup(group_summary: pd.DataFrame, metric: str, pet_kind: str) -> pd.Series:
    subset = group_summary.loc[group_summary['metric'].eq(metric) & group_summary['pet_kind'].eq(pet_kind)]
    if len(subset) != 1:
        raise ValueError(f'Missing group summary for {metric}/{pet_kind}')
    return subset.iloc[0]

def draw_metric_panel(ax: plt.Axes, trajectories: pd.DataFrame, group_summary: pd.DataFrame, metrics: list[tuple[str, str]], ylabel: str, label: str, ylim: tuple[float, float]) -> mpl.text.Text:
    offsets = {'PET_L10': -0.14, 'PET_L20': 0.14}
    for metric_index, (metric, _) in enumerate(metrics):
        for pet_kind in LENGTH_ORDER:
            subset = trajectories.loc[trajectories['pet_kind'].eq(pet_kind), metric].to_numpy(dtype=float)
            jitter = np.linspace(-0.075, 0.005, len(subset))
            x_values = metric_index + offsets[pet_kind] + jitter
            ax.scatter(x_values, subset, s=22, marker=MARKERS[pet_kind], facecolor=COLORS[pet_kind], edgecolor=BACKGROUND, linewidth=0.6, alpha=0.95, zorder=4)
            summary = summary_lookup(group_summary, metric, pet_kind)
            x_mean = metric_index + offsets[pet_kind] + 0.065
            ax.errorbar(x_mean, float(summary['mean']), yerr=np.asarray([[float(summary['mean'] - summary['bootstrap_95_low'])], [float(summary['bootstrap_95_high'] - summary['mean'])]]), fmt='D', markersize=2.8, markerfacecolor=INK, markeredgecolor=INK, color=INK, ecolor=INK, elinewidth=0.8, capsize=2.6, capthick=0.75, zorder=5)
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([item[1] for item in metrics])
    ax.set_ylabel(ylabel)
    ax.set_xlim(-0.2425, 1.2425)
    ax.set_ylim(*ylim)
    ax.grid(False)
    ax.spines[['top', 'right']].set_visible(False)
    return panel_label(ax, label)

def build_figure() -> None:
    apply_style()
    frames, trajectories, group_summary, references = load_data()
    fig = plt.figure(figsize=(7.0, 6.0))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.08, 1.0], height_ratios=[1.28, 0.8], left=0.1, right=0.985, bottom=0.15, top=0.94, wspace=0.42, hspace=0.4)
    ax_a = fig.add_subplot(grid[0, 0])
    labels: list[tuple[plt.Axes, mpl.text.Text]] = [(ax_a, draw_reference_space(ax_a, frames, references))]
    strip_grid = grid[0, 1].subgridspec(5, 1, hspace=0.6)
    strip_axes = [fig.add_subplot(strip_grid[index, 0]) for index in range(5)]
    labels.append((strip_axes[0], draw_time_strips(strip_axes, frames, references)))
    ax_c = fig.add_subplot(grid[1, 0])
    labels.append((ax_c, draw_metric_panel(ax_c, trajectories, group_summary, [('w156_contact_fraction', 'W156 contact'), ('same_unit_coengagement_fraction', 'Same-repeat\nco-engagement')], 'Trajectory occupancy', 'c', (0.0, 1.08))))
    ax_d = fig.add_subplot(grid[1, 1])
    labels.append((ax_d, draw_metric_panel(ax_d, trajectories, group_summary, [('chi2_circular_sd_deg', '$\\chi_2$ circular SD'), ('chi2_q90_span_deg', '$\\chi_2$ 5th–95th span')], 'Degrees', 'd', (0.0, 34.0))))
    legend_handles = [Line2D([0], [0], marker=MARKERS[pet_kind], color=COLORS[pet_kind], markerfacecolor=COLORS[pet_kind], markeredgecolor=BACKGROUND, markeredgewidth=0.6, linestyle='none', label=f'{DISPLAY_NAMES[pet_kind]} (n={int(trajectories['pet_kind'].eq(pet_kind).sum())})') for pet_kind in LENGTH_ORDER]
    legend_handles.append(ax_c.containers[0])
    fig.legend(handles=legend_handles, labels=[*(handle.get_label() for handle in legend_handles[:2]), 'Mean and 95% CI'], loc='lower center', bbox_to_anchor=(0.5, 0.035), ncol=3, handlelength=1.6, columnspacing=1.3, borderaxespad=0, frameon=False)
    for ax, label in labels:
        align_panel_label_to_yticks(fig, ax, label)
    return fig
FRAME_SOURCE = ROOT / '02_postprocessing/results/frame_observables/frame_level_w156_rotamer_contacts.csv.gz'
TRAJECTORY_SOURCE = ROOT / '02_postprocessing/results/trajectory_summaries/trajectory_rotamer_summary.csv'
GROUP_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/length_group_summary.csv'
REFERENCE_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/reference_structure_angles.csv'
MIXTURE_SOURCE = ROOT / '02_postprocessing/results/aggregate_statistics/mixture_model_sensitivity.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S04', build_figure)
