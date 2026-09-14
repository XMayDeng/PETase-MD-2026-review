#!/usr/bin/env python3
"""Build the JCIM table-of-contents graphic for the PET hydrolase panel.

The artwork is an original, evidence-led schematic. It is intentionally not a
quantitative plot and does not reproduce a manuscript figure. The visual
sequence is:

    distinct retained interfaces -> L10/L20 remodeling -> local targets

Outputs are fixed at the ACS TOC maximum size (3.25 x 1.75 inches) and 300 dpi.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, PathPatch, Polygon
from matplotlib.path import Path as MplPath
from PIL import Image


HERE = Path(__file__).resolve().parent
PREVIEW = HERE / "TOC_graphic_schematic.png"
FINAL_TIFF = HERE / "TOC_graphic_schematic.tif"

WIDTH_MM = 82.55
HEIGHT_MM = 44.45
WIDTH_IN = WIDTH_MM / 25.4
HEIGHT_IN = HEIGHT_MM / 25.4
DPI = 300
EXPECTED_PIXELS = (975, 525)

INK = "#20242B"
MUTED = "#66717D"
PROTEIN = "#E5EBF0"
PROTEIN_EDGE = "#8393A3"
PET_L10 = "#2684B7"
PET_L20 = "#D66A1F"
PET_GOLD = "#E9A321"
HBOND = "#168CA8"
AROMATIC = "#8265A8"

# Interface-class colors are shared with the manuscript figures.
CATALYTIC = "#F07F4F"
W_LOOP = "#F2C94C"
CLEFT = "#5B9BD5"
FLAP = "#D98BC3"
SUBSITE = "#70AD47"
OTHER = "#A6A6A6"


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 7.0,
            "font.weight": "regular",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "text.color": INK,
            "axes.edgecolor": INK,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "lines.solid_capstyle": "round",
            "lines.solid_joinstyle": "round",
        }
    )


def protein_path(cx: float, cy: float, rx: float, ry: float, phase: float) -> MplPath:
    """Return a smooth, slightly irregular closed path for a protein surface."""

    theta = np.linspace(0, 2 * np.pi, 97)
    modulation = (
        1.0
        + 0.10 * np.sin(3 * theta + phase)
        + 0.055 * np.sin(5 * theta - 0.7 * phase)
        + 0.025 * np.cos(7 * theta + phase)
    )
    x = cx + rx * modulation * np.cos(theta)
    y = cy + ry * modulation * np.sin(theta)
    vertices = np.column_stack([x, y])
    codes = np.full(len(vertices), MplPath.LINETO, dtype=np.uint8)
    codes[0] = MplPath.MOVETO
    codes[-1] = MplPath.CLOSEPOLY
    vertices[-1] = vertices[0]
    return MplPath(vertices, codes)


def draw_protein(
    ax: plt.Axes,
    cx: float,
    cy: float,
    rx: float,
    ry: float,
    phase: float,
    hotspots: list[tuple[float, str]],
    *,
    hotspot_size: float = 0.010,
    alpha: float = 1.0,
) -> None:
    patch = PathPatch(
        protein_path(cx, cy, rx, ry, phase),
        facecolor=PROTEIN,
        edgecolor=PROTEIN_EDGE,
        linewidth=0.85,
        alpha=alpha,
        zorder=2,
    )
    ax.add_patch(patch)

    # A pale inner contour suggests a cleft without claiming atomistic geometry.
    inner = PathPatch(
        protein_path(cx + 0.006 * np.cos(phase), cy, rx * 0.57, ry * 0.50, phase + 0.8),
        facecolor="none",
        edgecolor="white",
        linewidth=1.0,
        alpha=0.92,
        zorder=3,
    )
    ax.add_patch(inner)

    for angle, color in hotspots:
        rad = np.deg2rad(angle)
        hx = cx + rx * 0.92 * np.cos(rad)
        hy = cy + ry * 0.92 * np.sin(rad)
        ax.add_patch(
            Circle(
                (hx, hy),
                hotspot_size,
                facecolor=color,
                edgecolor=INK,
                linewidth=0.55,
                zorder=5,
            )
        )


def draw_pet_chain(
    ax: plt.Axes,
    x0: float,
    y0: float,
    units: int,
    *,
    color: str,
    spacing: float,
    scale: float,
    tilt: float = 0.0,
    zorder: int = 8,
) -> None:
    """Draw a compact PET-chain icon using linked aromatic repeat motifs."""

    angle = np.deg2rad(tilt)
    direction = np.array([np.cos(angle), np.sin(angle)])
    normal = np.array([-np.sin(angle), np.cos(angle)])
    centers: list[np.ndarray] = []
    for idx in range(units):
        wave = 0.010 * np.sin(idx * 1.25)
        center = np.array([x0, y0]) + idx * spacing * direction + wave * normal
        centers.append(center)
        angles = np.linspace(0, 2 * np.pi, 7)[:-1] + np.pi / 6 + angle
        points = np.column_stack(
            [center[0] + scale * np.cos(angles), center[1] + scale * 0.72 * np.sin(angles)]
        )
        ax.add_patch(
            Polygon(
                points,
                closed=True,
                fill=False,
                edgecolor=color,
                linewidth=1.05,
                zorder=zorder,
            )
        )
        if idx:
            p0 = centers[idx - 1] + 0.80 * scale * direction
            p1 = center - 0.80 * scale * direction
            ax.plot(
                [p0[0], p1[0]],
                [p0[1], p1[1]],
                color=color,
                linewidth=1.05,
                zorder=zorder,
            )


def draw_left_panel(ax: plt.Axes) -> None:
    ax.text(
        0.135,
        0.953,
        "Distinct\ninterfaces",
        ha="center",
        va="top",
        fontsize=6.6,
        fontweight="bold",
        linespacing=0.82,
    )

    icons = [
        ("IsPETase", 0.066, 0.705, 0.1, [(25, CATALYTIC), (126, W_LOOP), (210, CLEFT)]),
        ("TfCut1", 0.198, 0.705, 0.8, [(10, CATALYTIC), (110, W_LOOP), (225, CLEFT)]),
        ("LCC", 0.066, 0.420, 1.5, [(45, CATALYTIC), (145, W_LOOP), (245, CLEFT)]),
        ("FoCut5a", 0.198, 0.420, 2.2, [(60, FLAP), (95, FLAP), (135, FLAP)]),
        ("HiC", 0.066, 0.145, 2.9, [(315, CLEFT), (350, CLEFT), (25, CLEFT)]),
        ("PHL7", 0.198, 0.145, 3.6, [(175, SUBSITE), (215, SUBSITE), (255, SUBSITE)]),
    ]
    for name, x, y, phase, hotspots in icons:
        draw_protein(ax, x, y + 0.028, 0.046, 0.087, phase, hotspots, hotspot_size=0.008)
        ax.text(x, y - 0.082, name, ha="center", va="top", fontsize=6.25)


def draw_center_panel(ax: plt.Axes) -> None:
    # L10 state.
    draw_protein(
        ax,
        0.382,
        0.525,
        0.076,
        0.235,
        0.55,
        [(80, W_LOOP), (215, CLEFT)],
        hotspot_size=0.011,
    )
    draw_pet_chain(
        ax,
        0.329,
        0.565,
        3,
        color=PET_L10,
        spacing=0.027,
        scale=0.014,
        tilt=-8,
    )
    ax.text(0.382, 0.855, "L10", color=PET_L10, ha="center", fontsize=7.2, fontweight="bold")

    # L20 state. Additional/repositioned hotspots deliberately differ from L10.
    draw_protein(
        ax,
        0.602,
        0.525,
        0.076,
        0.235,
        1.15,
        [(32, CATALYTIC), (122, W_LOOP), (205, CLEFT), (286, OTHER)],
        hotspot_size=0.011,
    )
    draw_pet_chain(
        ax,
        0.524,
        0.586,
        6,
        color=PET_L20,
        spacing=0.025,
        scale=0.013,
        tilt=-10,
    )
    ax.text(0.602, 0.855, "L20", color=PET_L20, ha="center", fontsize=7.2, fontweight="bold")

    ax.add_patch(
        FancyArrowPatch(
            (0.466, 0.525),
            (0.515, 0.525),
            arrowstyle="-|>",
            mutation_scale=8.5,
            linewidth=0.9,
            color=INK,
            zorder=10,
        )
    )
    ax.text(
        0.492,
        0.186,
        "Enzyme-specific\ninterface remodeling",
        ha="center",
        va="center",
        fontsize=7.2,
        fontweight="bold",
        linespacing=0.95,
    )


def draw_hbond_icon(ax: plt.Axes, x: float, y: float, label: str) -> None:
    ax.plot([x, x + 0.045], [y, y + 0.018], color=HBOND, linewidth=1.15, zorder=8)
    ax.plot(
        [x + 0.045, x + 0.087],
        [y + 0.018, y - 0.004],
        color=HBOND,
        linewidth=0.9,
        linestyle=(0, (2.0, 1.6)),
        zorder=8,
    )
    ax.add_patch(Circle((x, y), 0.008, facecolor=HBOND, edgecolor=INK, linewidth=0.55, zorder=9))
    ax.add_patch(
        Circle((x + 0.091, y - 0.006), 0.008, facecolor=PET_GOLD, edgecolor=INK, linewidth=0.55, zorder=9)
    )
    ax.text(x - 0.004, y + 0.045, label, ha="left", va="bottom", fontsize=6.2, fontweight="bold")


def draw_aromatic_icon(ax: plt.Axes, x: float, y: float, label: str) -> None:
    angles = np.linspace(0, 2 * np.pi, 7)[:-1] + np.pi / 6
    for shift, color in [(0.0, AROMATIC), (0.054, PET_GOLD)]:
        points = np.column_stack(
            [x + shift + 0.020 * np.cos(angles), y + 0.016 * np.sin(angles)]
        )
        ax.add_patch(
            Polygon(points, closed=True, fill=False, edgecolor=color, linewidth=1.05, zorder=8)
        )
    ax.plot(
        [x + 0.020, x + 0.034],
        [y, y],
        color=AROMATIC,
        linewidth=0.9,
        linestyle=(0, (1.8, 1.5)),
        zorder=8,
    )
    ax.text(x - 0.020, y + 0.041, label, ha="left", va="bottom", fontsize=6.2, fontweight="bold")


def draw_right_panel(ax: plt.Axes) -> None:
    ax.text(
        0.850,
        0.935,
        "Local targets",
        ha="center",
        va="center",
        fontsize=7.4,
        fontweight="bold",
    )

    for y0, height in [(0.485, 0.305), (0.115, 0.285)]:
        ax.add_patch(
            FancyBboxPatch(
                (0.725, y0),
                0.250,
                height,
                boxstyle="round,pad=0.008,rounding_size=0.018",
                facecolor="#F7F9FB",
                edgecolor="#B6C0CA",
                linewidth=0.75,
                zorder=1,
            )
        )

    ax.text(0.742, 0.748, "TfCut1", ha="left", fontsize=7.0, fontweight="bold")
    draw_hbond_icon(ax, 0.750, 0.635, "T64/S67")
    draw_aromatic_icon(ax, 0.875, 0.530, "F210")

    ax.text(0.742, 0.360, "HiC", ha="left", fontsize=7.0, fontweight="bold")
    draw_hbond_icon(ax, 0.804, 0.245, "T166")


def build_figure() -> plt.Figure:
    configure_style()
    fig = plt.Figure(figsize=(WIDTH_IN, HEIGHT_IN), dpi=DPI, facecolor="white")
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    draw_left_panel(ax)
    draw_center_panel(ax)
    draw_right_panel(ax)

    # Lightweight progression arrows keep the reading order unambiguous.
    for x0, x1 in [(0.268, 0.301), (0.686, 0.714)]:
        ax.add_patch(
            FancyArrowPatch(
                (x0, 0.52),
                (x1, 0.52),
                arrowstyle="-|>",
                mutation_scale=7.5,
                linewidth=0.8,
                color=MUTED,
                zorder=12,
            )
        )
    return fig


def export(fig: plt.Figure) -> None:
    fig.canvas.draw()
    rgba = np.asarray(fig.canvas.buffer_rgba())
    if (rgba.shape[1], rgba.shape[0]) != EXPECTED_PIXELS:
        raise RuntimeError(
            f"Unexpected canvas size {(rgba.shape[1], rgba.shape[0])}; expected {EXPECTED_PIXELS}"
        )

    rgb = Image.fromarray(rgba, mode="RGBA").convert("RGB")
    rgb.save(PREVIEW, format="PNG", dpi=(DPI, DPI), optimize=True)
    rgb.save(FINAL_TIFF, format="TIFF", dpi=(DPI, DPI), compression="tiff_lzw")

    with Image.open(PREVIEW) as preview, Image.open(FINAL_TIFF) as final:
        if preview.size != EXPECTED_PIXELS or final.size != EXPECTED_PIXELS:
            raise RuntimeError("Exported image dimensions changed unexpectedly")
        if preview.mode != "RGB" or final.mode != "RGB":
            raise RuntimeError("ACS deliverables must be opaque RGB images")


def main() -> None:
    fig = build_figure()
    export(fig)
    plt.close(fig)
    print(f"Wrote {PREVIEW}")
    print(f"Wrote {FINAL_TIFF}")


if __name__ == "__main__":
    raise SystemExit("Drawing source only. Use make_figure.py --output NEW_DIR/TOC_graphic.png")
