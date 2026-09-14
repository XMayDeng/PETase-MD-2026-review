from __future__ import annotations
import matplotlib as mpl
from base_typography import TYPOGRAPHY, matplotlib_typography_rc

INK = '#20242B'

BACKGROUND = '#FFFFFF'

PROTEIN_ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']

def apply_style(base_font_size: float | None=None) -> None:
    typography = matplotlib_typography_rc()
    if base_font_size is not None:
        typography.update({'font.size': base_font_size, 'axes.titlesize': base_font_size + 1.0, 'axes.labelsize': base_font_size, 'xtick.labelsize': max(6.0, base_font_size - 0.5), 'ytick.labelsize': max(6.0, base_font_size - 0.5), 'legend.fontsize': max(6.0, base_font_size - 0.5)})
    mpl.rcParams.update({**typography, 'axes.linewidth': 0.6, 'xtick.major.width': 0.55, 'ytick.major.width': 0.55, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'axes.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'hatch.linewidth': 0.35})
