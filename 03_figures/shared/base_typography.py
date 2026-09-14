from dataclasses import dataclass

@dataclass(frozen=True)
class FigureTypography:
    """Target text hierarchy in points at final publication size."""
    panel_label: float = 8.5
    panel_title: float = 7.5
    axis_label: float = 7.0
    category_label: float = 6.5
    tick_label: float = 6.5
    legend: float = 6.5
    compact_legend: float = 5.8
    annotation: float = 6.0
    minor_annotation: float = 5.8
    structure_residue: float = 7.5
TYPOGRAPHY = FigureTypography()
FONT_FAMILY = {'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans']}

def matplotlib_typography_rc() -> dict[str, object]:
    """Return shared matplotlib defaults for final-size manuscript artwork."""
    return {**FONT_FAMILY, 'font.size': TYPOGRAPHY.category_label, 'axes.titlesize': TYPOGRAPHY.panel_title, 'axes.labelsize': TYPOGRAPHY.axis_label, 'xtick.labelsize': TYPOGRAPHY.tick_label, 'ytick.labelsize': TYPOGRAPHY.tick_label, 'legend.fontsize': TYPOGRAPHY.legend, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42}
