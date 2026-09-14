"""Shared final-physical-size typography and raster export for manuscript figures."""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path
import matplotlib as mpl
from matplotlib.collections import Collection
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator
from matplotlib.transforms import Bbox
import numpy as np
DPI = 600
SIZES = {'Figure_01': (6.5, 4.1), 'Figure_02': (3.33, 2.55), 'Figure_03': (6.5, 5.85), 'Figure_S02': (6.5, 4.65), 'Figure_04': (6.5, 4.95), 'Figure_05': (6.5, 3.36), 'Figure_06': (6.5, 5.3), 'Figure_S01': (6.5, 3.7), 'Figure_S03': (6.5, 5.2), 'Figure_S04': (6.5, 5.5), 'Figure_S06': (6.5, 3.15), 'Figure_S05': (6.5, 3.4), 'Figure_S07': (6.5, 5.35), 'Figure_S08': (6.5, 6.0), 'Figure_S09': (6.5, 4.9)}
FONTS = {'normal': 7.0, 'axis': 7.5, 'title': 8.0, 'panel': 9.0, 'ancillary': 6.5}
LINE_FLOOR = 0.6
BOTTOM_WHITESPACE_TRIM_IN = {'Figure_05': 0.24}
VISIBLE_BOTTOM_PADDING_PT = {'Figure_04': 7.5, 'Figure_05': 7.5}

def _array_hash(values) -> str:
    data = np.ma.asarray(values)
    if data.dtype.kind in 'OUS':
        payload = repr(data.tolist()).encode()
    else:
        payload = np.asarray(data.filled(np.nan) if np.ma.is_masked(data) else data, dtype='<f8').tobytes()
    return hashlib.sha256(payload).hexdigest()

def evidence_signature(fig) -> str:
    """Hash numerical artists, not presentation widths, fonts, or canvas geometry."""
    records = []
    for ax in fig.axes:
        records.append(['limits', list(ax.get_xlim()), list(ax.get_ylim())])
        for line in ax.lines:
            records.append(['line', _array_hash(line.get_xdata()), _array_hash(line.get_ydata()), str(line.get_marker()), str(line.get_color()), str(line.get_markerfacecolor())])
        for item in ax.collections:
            records.append(['collection', _array_hash(item.get_offsets()), [_array_hash(path.vertices) for path in item.get_paths()], _array_hash(item.get_array()) if item.get_array() is not None else None])
        for item in ax.patches:
            if isinstance(item, Rectangle):
                records.append(['rectangle', list(item.get_xy()), item.get_width(), item.get_height(), str(item.get_facecolor()), item.get_hatch()])
        for item in ax.images:
            records.append(['image', _array_hash(item.get_array()), list(item.get_extent()), list(item.get_clim())])
        records.extend((['annotation', ' '.join(text.get_text().split())] for text in ax.texts if text.get_text().strip()))
    return hashlib.sha256(json.dumps(records, sort_keys=True, default=float).encode()).hexdigest()

def _visible_texts(fig):
    hidden = set()
    for ax in fig.axes:
        for axis, limits in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            low, high = sorted(limits)
            for tick in [*axis.get_major_ticks(), *axis.get_minor_ticks()]:
                if tick.get_loc() < low - 1e-10 or tick.get_loc() > high + 1e-10:
                    hidden.update([tick.label1, tick.label2])
        if not ax.axison:
            hidden.update(ax.get_xticklabels(which='both'))
            hidden.update(ax.get_yticklabels(which='both'))
            hidden.update([ax.xaxis.label, ax.yaxis.label])
    return [text for text in fig.findobj(match=Text) if text.get_visible() and text.get_text().strip() and (text not in hidden)]

def apply_final_typography(fig, key: str) -> None:
    width, height = SIZES[key]
    height += BOTTOM_WHITESPACE_TRIM_IN.get(key, 0.0)
    fig.set_size_inches(width, height, forward=True)
    fig.set_dpi(DPI)
    axis_labels = {label for ax in fig.axes for label in (ax.xaxis.label, ax.yaxis.label)}
    axis_labels.update((label for label in (getattr(fig, '_supxlabel', None), getattr(fig, '_supylabel', None)) if label is not None))
    titles = {title for ax in fig.axes for title in (ax.title, ax._left_title, ax._right_title)}
    for ax in fig.axes:
        titles.update((text for text in ax.texts if text.get_fontweight() in ('bold', 'semibold', 700) and text.get_position()[1] > 1.0 and (text.get_transform() == ax.transAxes or text.get_transform() == ax.get_xaxis_transform())))
    for text in fig.findobj(match=Text):
        value = text.get_text().strip()
        role = 'normal'
        if re.fullmatch('\\([a-z]\\)', value):
            role = 'panel'
        elif text in axis_labels:
            role = 'axis'
        elif text in titles or (text.get_fontweight() in ('bold', 'semibold', 700) and text.get_fontsize() >= 7.0 and (not re.match('^[+−\\-\\d]', value))):
            role = 'title'
        elif value in ('Whole fold', 'Active-site detail') or re.fullmatch('L\\d+, [\\d.]+ ns', value):
            role = 'ancillary'
        role = getattr(text, '_publication_role', role)
        text.set_fontsize(FONTS[role])
        text.set_fontfamily('DejaVu Sans')
        text._publication_role = role
    for legend in fig.findobj(match=Legend):
        for text in legend.get_texts():
            text.set_fontsize(FONTS['normal'])
            text._publication_role = 'normal'
        legend._legend_box.sep = legend.labelspacing * FONTS['normal']
    for item in fig.findobj():
        if isinstance(item, Line2D):
            if item.get_linewidth() > 0:
                item.set_linewidth(max(LINE_FLOOR, item.get_linewidth()))
            if item.get_markeredgewidth() > 0:
                item.set_markeredgewidth(max(LINE_FLOOR, item.get_markeredgewidth()))
        elif isinstance(item, Patch):
            if item.get_linewidth() > 0:
                item.set_linewidth(max(LINE_FLOOR, item.get_linewidth()))
        elif isinstance(item, Collection):
            widths = np.asarray(item.get_linewidths())
            item.set_linewidths(np.where(widths > 0, np.maximum(widths, LINE_FLOOR), widths))
    mpl.rcParams['hatch.linewidth'] = LINE_FLOOR
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        if not ax.axison:
            continue
        ticks = [tick for tick in ax.get_yticklabels() if tick.get_visible() and tick.get_text().strip()]
        if not ticks:
            continue
        x = min((tick.get_window_extent(renderer).x0 for tick in ticks))
        x_axes = ax.transAxes.inverted().transform((x, ax.bbox.y0))[0]
        for text in ax.texts:
            if re.fullmatch('\\([a-z]\\)', text.get_text().strip()):
                text.set_x(x_axes)
                text.set_ha('left')
    for axes_group in getattr(fig, '_publication_center_groups', []):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        bounds = Bbox.union([ax.get_tightbbox(renderer) for ax in axes_group])
        shift = 0.5 - (bounds.x0 + bounds.x1) / (2 * fig.bbox.width)
        for ax in axes_group:
            position = ax.get_position()
            ax.set_position([position.x0 + shift, position.y0, position.width, position.height])
    if hasattr(fig, '_publication_panel_legend_layout'):
        from panel_legend_layout import refine_panel_legends
        refine_panel_legends(fig)
    fig.canvas.draw()

def save_publication_figure(fig, output, **kwargs) -> None:
    output = Path(output)
    match = re.match('(Figure_(?:S)?\\d{2})_', output.name)
    if match is None or match.group(1) not in SIZES:
        raise ValueError(f'Unregistered publication figure: {output}')
    key = match.group(1)
    fig.canvas.draw()
    for ax in fig.axes:
        ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
        ax.yaxis.set_major_locator(FixedLocator(ax.get_yticks()))
        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
    before = evidence_signature(fig)
    apply_final_typography(fig, key)
    after = evidence_signature(fig)
    if before != after:
        raise RuntimeError(f'{key}: numerical artists changed during typography normalization')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_records = []
    outside = []
    bottom_trim = BOTTOM_WHITESPACE_TRIM_IN.get(key, 0.0)
    width, height = SIZES[key]
    source_size = fig.get_size_inches().tolist()
    if key in VISIBLE_BOTTOM_PADDING_PT:
        pixels = np.asarray(fig.canvas.buffer_rgba())
        ink_rows = np.where(np.any(pixels[:, :, :3] != 255, axis=(1, 2)))[0]
        white_rows = pixels.shape[0] - 1 - int(ink_rows[-1])
        keep_rows = int(np.ceil(VISIBLE_BOTTOM_PADDING_PT[key] * DPI / 72))
        trim_rows = max(0, white_rows - keep_rows)
        bottom_trim = trim_rows / DPI
        height = (pixels.shape[0] - trim_rows) / DPI
    for text in _visible_texts(fig):
        bounds = text.get_window_extent(renderer)
        record = {'text': text.get_text(), 'size_pt': text.get_fontsize(), 'role': getattr(text, '_publication_role', 'normal')}
        text_records.append(record)
        if bounds.x0 < -0.5 or bounds.y0 < bottom_trim * DPI - 0.5 or bounds.x1 > fig.bbox.width + 0.5 or (bounds.y1 > fig.bbox.height + 0.5):
            outside.append(record)
    if bottom_trim:
        rows = round(bottom_trim * DPI)
        pixels = np.asarray(fig.canvas.buffer_rgba())
        if not np.all(pixels[-rows:, :, :3] == 255):
            raise RuntimeError(f'{key}: requested bottom trim contains visible artwork')
    export_height = np.nextafter(height, np.inf) if bottom_trim else height
    fig.savefig(output, dpi=DPI, facecolor='white', format='png', bbox_inches=Bbox.from_bounds(0, bottom_trim, width, export_height), pad_inches=0)
    report = {'figure': key, 'output': str(output), 'dimensions_inches': [width, height], 'source_canvas_inches': source_size, 'bottom_whitespace_trim_inches': bottom_trim, 'visible_bottom_padding_pt': VISIBLE_BOTTOM_PADDING_PT.get(key), 'dpi': DPI, 'expected_pdf_ppi': DPI, 'font_roles_pt': FONTS, 'minimum_visible_vector_line_pt': LINE_FLOOR, 'embedded_structure_text': 'none; 22 residue labels are Matplotlib callouts with exact ChimeraX-projected anchors' if key == 'Figure_01' else 'none', 'numerical_artist_signature_before': before, 'numerical_artist_signature_after': after, 'evidence_preserved': before == after, 'minimum_matplotlib_font_pt': min((record['size_pt'] for record in text_records), default=None), 'text_outside_canvas': outside, 'texts': text_records, 'relative_panel_layout': getattr(fig, '_publication_relative_layout', None), 'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest()}
    print(f'{key}: {width:g} × {height:g} in; {len(text_records)} text labels; outside canvas: {len(outside)}')
    return report
