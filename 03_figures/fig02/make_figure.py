from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpatheffects
from main_style import INK, PROTEIN_ORDER, apply_style, TYPOGRAPHY
from revised_figure_layout import contrasting_text_color
from data_adapter import retention_tables as _retention_tables

def build_figure() -> plt.Figure:
    apply_style()
    _, condition = _retention_tables()
    conditions = ['L10/bi', 'L10/head', 'L10/tail', 'L20/bi', 'L20/head', 'L20/tail', 'L4/single']
    matrix = condition.pivot(index='protein_short', columns='condition', values='mean_retention_fraction').reindex(index=PROTEIN_ORDER, columns=conditions)
    fig, ax = plt.subplots(figsize=(3.55, 2.35))
    cmap = mcolors.LinearSegmentedColormap.from_list('retention_green', ['#F7FAF7', '#B8D8B8', '#2F6B45'])
    image = ax.imshow(matrix.to_numpy(), cmap=cmap, vmin=0, vmax=1, aspect='auto')
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = float(matrix.iloc[i, j])
            ax.text(j, i, f'{value:.2f}', ha='center', va='center', fontsize=TYPOGRAPHY.annotation, color=contrasting_text_color(image.cmap(image.norm(value)), dark=INK), fontweight='bold' if value >= 0.995 else 'normal')
    anomaly_i = PROTEIN_ORDER.index('PHL7')
    anomaly_j = conditions.index('L20/head')
    ax.add_patch(mpatches.Rectangle((anomaly_j - 0.5, anomaly_i - 0.5), 1, 1, fill=False, edgecolor='#B21F2D', linewidth=1.5, zorder=10, clip_on=False, path_effects=[mpatheffects.Stroke(linewidth=2.3, foreground='white'), mpatheffects.Normal()]))
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=35, ha='right', rotation_mode='anchor')
    ax.set_yticks(range(len(PROTEIN_ORDER)))
    ax.set_yticklabels(PROTEIN_ORDER)
    ax.tick_params(length=0, pad=1.5)
    ax.set_xticks(np.arange(-0.5, len(conditions), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(PROTEIN_ORDER), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.55)
    ax.tick_params(which='minor', bottom=False, left=False)
    for boundary in (2.5, 5.5):
        ax.axvline(boundary, color='white', linewidth=1.25)
    cbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.035, ticks=[0, 0.5, 1])
    cbar.set_label('Mean retention fraction', labelpad=3)
    cbar.outline.set_linewidth(0.5)
    fig.subplots_adjust(left=0.205, right=0.85, top=0.965, bottom=0.245)
    return fig
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_02', build_figure)
