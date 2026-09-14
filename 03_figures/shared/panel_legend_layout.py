"""Place legends relative to their owning panels at final physical size."""

import re

import numpy as np
from matplotlib.transforms import Bbox


def _reanchor_header(ax, renderer):
    ticks = [text for text in ax.get_yticklabels()
             if text.get_visible() and text.get_text().strip()]
    if not ticks:
        return
    left = min(text.get_window_extent(renderer).x0 for text in ticks)
    x = ax.transAxes.inverted().transform((left, ax.bbox.y0))[0]
    for text in ax.texts:
        if re.fullmatch(r"\([a-z]\)", text.get_text().strip()):
            text.set_x(x)
            text.set_ha("left")


def _fit_panel(fig, panel, lower, right_inset=0.0):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    reference = Bbox.union([ax.get_tightbbox(renderer) for ax in lower])
    target_left = reference.x0
    target_right = reference.x1 - right_inset * fig.bbox.width
    for _ in range(3):
        renderer = fig.canvas.get_renderer()
        bounds = panel.get_tightbbox(renderer)
        left = panel.bbox.x0 + target_left - bounds.x0
        right = panel.bbox.x1 + target_right - bounds.x1
        position = panel.get_position()
        panel.set_position([left / fig.bbox.width, position.y0,
                            (right - left) / fig.bbox.width, position.height])
        _reanchor_header(panel, renderer)
        fig.canvas.draw()


def refine_panel_legends(fig):
    """Keep panel geometry and apply only explicitly requested legend changes."""
    spec = fig._publication_panel_legend_layout
    panel = spec["panel"]
    if spec["key"] in ("Figure_03", "Figure_04"):
        _fit_panel(fig, panel, spec["lower"],
                   right_inset=0.025 if spec["key"] == "Figure_04" else 0.0)
    fig._publication_relative_layout = {
        "panel_data_width_fraction": panel.get_position().width,
        "legend_layout": "restored to the pre-relocation position and expanded arrangement",
    }
    if "legend_gap_pt" not in spec:
        return
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend = spec.get("legend", panel.get_legend())
    gap = spec["legend_gap_pt"]
    if spec["key"] == "Figure_04":
        box = panel.get_position()
        center = (box.x0 + box.x1) / 2
        reference_bottom = panel.xaxis.label.get_window_extent(renderer).y0
        top = (reference_bottom - gap * fig.dpi / 72) / fig.bbox.height
        anchor = (box.x0, top, box.width, 0.0)
    else:
        image_axes = panel.child_axes[0]
        pixels = np.asarray(image_axes.images[0].get_array())
        rows = np.where(np.any(pixels[:, :, :3] != 255, axis=(1, 2)))[0]
        reference_bottom = image_axes.transData.transform((0, rows[-1] + 0.5))[1]
        center = (image_axes.bbox.x0 + image_axes.bbox.x1) / (2 * fig.bbox.width)
        anchor = (center, (reference_bottom - gap * fig.dpi / 72) / fig.bbox.height)
    legend.set_loc("upper center")
    legend.borderaxespad = 0.0
    legend.set_bbox_to_anchor(anchor, transform=fig.transFigure)
    fig.canvas.draw()
    bounds = legend.get_window_extent(fig.canvas.get_renderer())
    fig._publication_relative_layout.update({
        "legend_layout": "existing arrangement moved closer and centered on its own panel",
        "panel_center_x": center,
        "legend_center_x": (bounds.x0 + bounds.x1) / (2 * fig.bbox.width),
        "gap_pt": (reference_bottom - bounds.y1) * 72 / fig.dpi,
    })
    if "lower_row_gap_pt" in spec:
        renderer = fig.canvas.get_renderer()
        row_top = max(ax.get_tightbbox(renderer).y1 for ax in spec["lower"])
        desired_top = bounds.y0 - spec["lower_row_gap_pt"] * fig.dpi / 72
        shift = (desired_top - row_top) / fig.bbox.height
        for ax in spec["lower_row"]:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + shift, box.width, box.height])
        fig.canvas.draw()
        row_top = max(ax.get_tightbbox(fig.canvas.get_renderer()).y1 for ax in spec["lower"])
        fig._publication_relative_layout["legend_to_lower_row_pt"] = (bounds.y0 - row_top) * 72 / fig.dpi
