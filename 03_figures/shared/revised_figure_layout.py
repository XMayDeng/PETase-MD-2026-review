"""Shared layout and annotation-contrast helpers for manuscript figures."""

from __future__ import annotations

from collections.abc import Iterable

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


PANEL_LABEL_SIZE = 8.5
ANNOTATION_COLOR_POLICY = "maximum sRGB luminance contrast against the rendered fill"


def relative_luminance(color, *, background="white") -> float:
    """Return linear-sRGB luminance, compositing transparency over the canvas."""
    red, green, blue, alpha = mcolors.to_rgba(color)
    canvas = mcolors.to_rgb(background)
    rgb = [alpha * value + (1.0 - alpha) * base
           for value, base in zip((red, green, blue), canvas)]
    linear = [value / 12.92 if value <= 0.04045
              else ((value + 0.055) / 1.055) ** 2.4 for value in rgb]
    return sum(weight * value for weight, value in zip((0.2126, 0.7152, 0.0722), linear))


def contrasting_text_color(fill, *, dark="#20242B", light="white", background="white"):
    """Choose the more legible text color from the actual fill, not its data value."""
    fill_luminance = relative_luminance(fill, background=background)

    def contrast(color):
        text_luminance = relative_luminance(color, background=background)
        return ((max(fill_luminance, text_luminance) + 0.05)
                / (min(fill_luminance, text_luminance) + 0.05))

    return light if contrast(light) > contrast(dark) else dark


def panel_label(
    ax: plt.Axes,
    label: str,
    *,
    x: float = 0.0,
    y: float = 1.02,
    va: str = "bottom",
) -> mpl.text.Text:
    """Place a parenthesized lower-case panel label at a panel's upper left."""
    normalized = label.strip().strip("()").lower()
    return ax.text(
        x,
        y,
        f"({normalized})",
        transform=ax.transAxes,
        ha="left",
        va=va,
        fontsize=PANEL_LABEL_SIZE,
        fontweight="bold",
        zorder=30,
    )


def align_panel_labels_to_yticklabels(
    fig: plt.Figure,
    panel_headers: Iterable[tuple[plt.Axes, mpl.text.Text]],
) -> None:
    """Align each panel label with the left edge of its y-tick label column."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax, label in panel_headers:
        yticklabels = [
            tick
            for tick in ax.get_yticklabels()
            if tick.get_visible() and tick.get_text().strip()
        ]
        if not yticklabels:
            continue
        left_display = min(
            tick.get_window_extent(renderer=renderer).x0 for tick in yticklabels
        )
        left_axes = ax.transAxes.inverted().transform((left_display, ax.bbox.y0))[0]
        label.set_x(left_axes)
