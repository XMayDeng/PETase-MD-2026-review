from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.patheffects as mpatheffects
import numpy as np
import pandas as pd
from PIL import Image
from revised_figure_layout import panel_label
ROOT = Path(__file__).resolve().parents[2]
INK = '#20242B'
MUTED = '#626A75'
BACKGROUND = '#FFFFFF'
PET_ORANGE = '#C8782A'
PROTEIN_ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']
CLASS_ORDER = ['Catalytic/oxyanion', 'W-loop/flexible or rigid pair', 'Cleft/rim patch', 'Flap/binding loop', 'Subsite/hotspot/stability', 'Other/candidate']
CLASS_DISPLAY = {'Catalytic/oxyanion': 'Catalytic/oxyanion', 'W-loop/flexible or rigid pair': 'W-loop/rigid pair', 'Cleft/rim patch': 'Cleft/rim patch', 'Flap/binding loop': 'Flap/binding loop', 'Subsite/hotspot/stability': 'Subsite/hotspot', 'Other/candidate': 'Other/candidate'}
CLASS_COLORS = {'Catalytic/oxyanion': '#E07B54', 'W-loop/flexible or rigid pair': '#E8C547', 'Cleft/rim patch': '#5B9BD5', 'Flap/binding loop': '#D98EC7', 'Subsite/hotspot/stability': '#70AD47', 'Other/candidate': '#A6A6A6'}
ANCHOR_MARKER_COLORS = ('#FF0000', '#00B000', '#0050FF', '#C000C0')
LABEL_OFFSETS_AXES = ((0.075, 0.055), (-0.075, 0.055), (0.075, -0.055), (-0.075, -0.055))

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'text.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def extract_label_anchors(path: Path, count: int) -> list[tuple[float, float]]:
    image = np.asarray(Image.open(path).convert('RGB'), dtype=float) / 255.0
    hsv = mcolors.rgb_to_hsv(image)
    anchors: list[tuple[float, float]] = []
    for marker_color in ANCHOR_MARKER_COLORS[:count]:
        target_rgb = np.asarray(mcolors.to_rgb(marker_color), dtype=float)
        target_hsv = mcolors.rgb_to_hsv(target_rgb.reshape(1, 1, 3))[0, 0]
        hue_distance = np.abs(hsv[:, :, 0] - target_hsv[0])
        hue_distance = np.minimum(hue_distance, 1.0 - hue_distance)
        mask = (hue_distance <= 0.035) & (hsv[:, :, 1] >= 0.45) & (hsv[:, :, 2] >= 0.18)
        rows, columns = np.where(mask)
        if len(rows) < 20:
            raise ValueError(f'Could not recover anchor color {marker_color} from {path}')
        anchors.append((float(np.median(columns) / (image.shape[1] - 1)), float(1.0 - np.median(rows) / (image.shape[0] - 1))))
    return anchors

def assemble_figure(top8: pd.DataFrame, selected: pd.DataFrame, panel_paths: dict[str, Path], anchor_paths: dict[str, Path]) -> None:
    selection_lookup = selected.set_index('protein')
    anchor_lookup: dict[str, dict[int, tuple[float, float]]] = {}
    for protein in PROTEIN_ORDER:
        protein_rows = top8[top8['protein_short'] == protein].sort_values('rank_within_protein')
        frame_contacts = {int(value) for value in str(selection_lookup.loc[protein, 'frame_contact_hotspots_4p5A']).split(';') if value and value != 'nan'}
        label_rows = protein_rows[protein_rows['residue'].astype(int).isin(frame_contacts)].head(4)
        anchors = extract_label_anchors(anchor_paths[protein], len(label_rows))
        anchor_lookup[protein] = {int(label_row.residue): anchor for label_row, anchor in zip(label_rows.itertuples(), anchors)}
    fig, axes = plt.subplots(2, 3, figsize=(7.0, 4.9))
    for panel_index, (ax, protein) in enumerate(zip(axes.flat, PROTEIN_ORDER)):
        protein_rows = top8[top8['protein_short'] == protein].sort_values('rank_within_protein')
        row = selection_lookup.loc[protein]
        frame_contacts = {int(value) for value in str(row['frame_contact_hotspots_4p5A']).split(';') if value and value != 'nan'}
        label_rows = protein_rows[protein_rows['residue'].astype(int).isin(frame_contacts)].head(4)
        ax.set_axis_off()
        structure_ax = ax.inset_axes([0.0, 0.055, 1.0, 0.82])
        structure_ax.imshow(Image.open(panel_paths[protein]).convert('RGB'))
        structure_ax.set_axis_off()
        for label_index, label_row in enumerate(label_rows.itertuples()):
            anchor_x, anchor_y = anchor_lookup[protein][int(label_row.residue)]
            offset_x, offset_y = LABEL_OFFSETS_AXES[label_index]
            label_x = float(np.clip(anchor_x + offset_x, 0.06, 0.94))
            label_y = float(np.clip(anchor_y + offset_y, 0.08, 0.92))
            label = structure_ax.annotate(label_row.residue_label, xy=(anchor_x, anchor_y), xycoords=structure_ax.transAxes, xytext=(label_x, label_y), textcoords=structure_ax.transAxes, ha='left' if offset_x > 0 else 'right', va='center', fontsize=6.0, color='#161A20', zorder=30)
            label.set_path_effects([mpatheffects.withStroke(linewidth=1.15, foreground=BACKGROUND)])
        panel_label(ax, chr(ord('a') + panel_index), x=0.012, y=0.985, va='top')
        ax.text(0.5, 0.975, protein, transform=ax.transAxes, ha='center', va='top', fontsize=7.5, fontweight='bold', zorder=30)
        condition = f'{str(row['pet_kind']).replace('PET_', '')}, {float(row['frame_time_ps']) / 1000.0:.2f} ns'
        ax.text(0.5, 0.012, condition, transform=ax.transAxes, ha='center', va='bottom', fontsize=5.8, color=MUTED)
    legend_handles = [mpatches.Patch(facecolor=CLASS_COLORS[cls], edgecolor='#4A4F57', label=CLASS_DISPLAY[cls]) for cls in CLASS_ORDER]
    legend_handles.extend([mpatches.Patch(facecolor=PET_ORANGE, edgecolor='#4A4F57', label='Local PET contact units'), mpatches.Patch(facecolor='#7B838D', edgecolor='#4A4F57', label='Opaque: contact in frame'), mpatches.Patch(facecolor='#C8CDD3', edgecolor='#4A4F57', label='Pale: pooled hotspot only')])
    fig.legend(handles=legend_handles, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.005), columnspacing=0.8, handlelength=1.1)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.13, wspace=0.018, hspace=0.075)
    return fig

def build_figure(rendered):
    apply_style()
    top8 = pd.read_csv(TOP8_PATH)
    top8 = top8[top8.is_top8.astype(bool)].copy()
    selected = pd.read_csv(SELECTION_PATH)
    if len(selected) != 6 or set(selected.protein) != set(PROTEIN_ORDER):
        raise ValueError('Expected one declared frame per enzyme')
    if not top8.groupby('protein_short').size().reindex(PROTEIN_ORDER).eq(8).all():
        raise ValueError('Expected eight hotspots per enzyme')
    if selected.frame_contact_hotspots_4p5A.isna().any():
        raise ValueError('Missing frame-contact labels')
    panels = {p: rendered / f'Figure_04_{p}_closeup.png' for p in PROTEIN_ORDER}
    anchors = {p: rendered / f'Figure_04_{p}_label_anchor_map.png' for p in PROTEIN_ORDER}
    return assemble_figure(top8, selected, panels, anchors)
TOP8_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/pooled_top8_contacts_unfiltered.csv'
SELECTION_PATH = ROOT / '02_postprocessing/results/structural_snapshots/representative_interface_frames.csv'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_S02', build_figure)
