from __future__ import annotations
import matplotlib as mpl
import matplotlib.pyplot as plt
PROTEIN_ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']
CONDITIONS = ['L10/bi', 'L10/head', 'L10/tail', 'L20/bi', 'L20/head', 'L20/tail', 'L4/single']
CATEGORIES = ['Detached', 'Partial', 'Retained + WARN', 'Retained + PASS']
CATEGORY_COLORS = ['#C7CCD2', '#E8C547', '#D98EC7', '#5B9BD5']
CATEGORY_HATCHES = ['////', '....', 'xxxx', '\\\\\\\\']
INK = '#20242B'
BACKGROUND = '#FFFFFF'

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

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
