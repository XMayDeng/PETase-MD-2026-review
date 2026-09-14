from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patheffects as mpatheffects
from s1_helpers import CATEGORIES, CATEGORY_COLORS, CATEGORY_HATCHES, CONDITIONS, INK, PROTEIN_ORDER, align_panel_label_to_yticks, apply_style, panel_label
from data_adapter import s1_tables as load_tables
from publication_style import apply_final_typography

def apply_trial_style() -> None:
    """Expose the inherited publication settings in this standalone source."""
    apply_style()
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def align_heatmap_to_bar_bounds(fig, ax_left, ax_right, cax, right_label):
    """Apply the author-approved S1 alignment at final publication size."""
    apply_final_typography(fig, "Figure_S01")
    fig.canvas.draw()
    left_bounds = ax_left.get_position().frozen()
    right_bounds = ax_right.get_position().frozen()
    colorbar_bounds = cax.get_position().frozen()
    label_position = fig.transFigure.inverted().transform(
        right_label.get_transform().transform(right_label.get_position()))
    first_bar = ax_left.containers[0].patches[0]
    top = first_bar.get_window_extent(fig.canvas.get_renderer()).transformed(
        fig.transFigure.inverted()).y1
    bottom = left_bounds.y0
    height = top - bottom
    ax_right.set_position([right_bounds.x0, bottom, right_bounds.width, height])
    cax.set_box_aspect(height * fig.get_figheight() /
                       (colorbar_bounds.width * fig.get_figwidth()))
    cax.set_position([colorbar_bounds.x0, bottom, colorbar_bounds.width, height])
    right_label.set_y((label_position[1] - bottom) / height)
    fig.canvas.draw()
    np.testing.assert_allclose(ax_left.get_position().bounds, left_bounds.bounds,
                               atol=1e-12, rtol=0)
    for axis in (ax_right, cax):
        bounds = axis.get_position()
        np.testing.assert_allclose([bounds.y0, bounds.y1], [bottom, top],
                                   atol=1e-12, rtol=0)
    return {
        "panel_a_bounds": list(left_bounds.bounds),
        "panel_b_bounds": list(ax_right.get_position().bounds),
        "colorbar_bounds": list(cax.get_position().bounds),
        "panel_a_first_bar_top": top,
        "panel_b_top_matches_first_bar": True,
        "panel_a_b_x_axis_baselines_aligned": True,
    }


def build_figure() -> None:
    apply_trial_style()
    _, counts, retained = load_tables()
    fig = plt.figure(figsize=(7.0, 3.55))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], left=0.105, right=0.945, top=0.94, bottom=0.285, wspace=0.34)
    ax_left = fig.add_subplot(grid[0, 0])
    ax_right = fig.add_subplot(grid[0, 1])
    y = np.arange(len(PROTEIN_ORDER))
    left = np.zeros(len(PROTEIN_ORDER))
    legend_handles = []
    for category, color, hatch in zip(CATEGORIES, CATEGORY_COLORS, CATEGORY_HATCHES, strict=True):
        values = counts[category].to_numpy(int)
        bars = ax_left.barh(y, values, left=left, height=0.64, color=color, edgecolor='#4A4F57', linewidth=0.45, label=category)
        for bar in bars:
            bar.set_hatch(hatch)
        legend_handles.append(bars[0])
        for yi, (value, start) in enumerate(zip(values, left, strict=True)):
            if value >= 3:
                label = ax_left.text(start + value / 2, yi, str(value), ha='center', va='center', fontsize=6.5, fontweight='bold')
                label.set_path_effects([mpatheffects.withStroke(linewidth=1.0, foreground=color)])
        left += values
    ax_left.set_yticks(y, PROTEIN_ORDER)
    ax_left.invert_yaxis()
    ax_left.set_xlim(0, 21)
    ax_left.set_xticks([0, 7, 14, 21])
    ax_left.set_xlabel('Trajectory count (21 per enzyme)')
    ax_left.grid(False)
    ax_left.spines[['top', 'right']].set_visible(False)
    left_label = panel_label(ax_left, 'a')
    cmap = mcolors.ListedColormap(['#F3F4F5', '#C9D9EA', '#7FA6CC', '#356B9B'])
    norm = mcolors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    image = ax_right.imshow(retained.to_numpy(), cmap=cmap, norm=norm, aspect='auto')
    for i in range(retained.shape[0]):
        for j in range(retained.shape[1]):
            value = int(retained.iloc[i, j])
            ax_right.text(j, i, str(value), ha='center', va='center', fontsize=6.5, fontweight='bold', color='white' if value == 3 else INK)
    ax_right.set_xticks(range(len(CONDITIONS)))
    ax_right.set_xticklabels(CONDITIONS, rotation=35, ha='right', rotation_mode='anchor')
    ax_right.set_yticks(range(len(PROTEIN_ORDER)), PROTEIN_ORDER)
    ax_right.tick_params(length=0, pad=1.8)
    ax_right.set_xticks(np.arange(-0.5, len(CONDITIONS), 1), minor=True)
    ax_right.set_yticks(np.arange(-0.5, len(PROTEIN_ORDER), 1), minor=True)
    ax_right.grid(which='minor', color='white', linewidth=0.7)
    ax_right.tick_params(which='minor', bottom=False, left=False)
    right_label = panel_label(ax_right, 'b')
    cbar = fig.colorbar(image, ax=ax_right, fraction=0.05, pad=0.035, ticks=[0, 1, 2, 3])
    cbar.set_label('Retained trajectories (of 3)', labelpad=3)
    cbar.outline.set_linewidth(0.5)
    ax_left.legend(handles=legend_handles, labels=CATEGORIES, loc='upper center', bbox_to_anchor=(0.5, -0.27), ncol=2, frameon=False, columnspacing=1.0, handlelength=1.1, handletextpad=0.4, borderaxespad=0)
    align_panel_label_to_yticks(fig, ax_left, left_label)
    align_panel_label_to_yticks(fig, ax_right, right_label)
    fig._publication_relative_layout = align_heatmap_to_bar_bounds(
        fig, ax_left, ax_right, cbar.ax, right_label)
    return fig
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S01', build_figure)
