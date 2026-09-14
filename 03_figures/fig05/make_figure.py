from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea, VPacker
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from PIL import Image
from revised_figure_layout import panel_label, align_panel_labels_to_yticklabels, contrasting_text_color
ROOT = Path(__file__).resolve().parents[2]
HIC_COLOR = '#3D7EBB'
HIC_DARK = '#1F5B91'
FOCUT5A_COLOR = '#D07832'
FOCUT5A_DARK = '#A95318'
INK = '#20242B'
MUTED = '#626A75'
GRID = '#D9DEE5'
BACKGROUND = '#FFFFFF'
ANCHOR_COLORS = ('#FF0000', '#00B000', '#0050FF')
FOCUS_PAIRS = [(28, 'S28', 44, 'S44'), (29, 'T29', 45, 'T45'), (70, 'F70', 87, 'A87'), (164, 'T164', 181, 'T181'), (165, 'G165', 182, 'G182'), (166, 'T166', 183, 'S183')]
STRUCTURE_ZOOM = 1.12
STRUCTURE_DOWN_SHIFT = 0.025
HEATMAP_HEIGHT_FRACTION = 0.81
STRUCTURAL_LABELS = (('S28', 'S44', (0.41, 0.38)), ('T29', 'T45', (0.27, 0.275)), ('F70', 'A87', (0.5, 0.525)), ('T164', 'T181', (0.81, 0.71)), ('G165', 'G182', (0.765, 0.59)), ('T166', 'S183', (0.82, 0.49)))

def apply_style() -> None:
    mpl.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7.0, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.5, 'axes.linewidth': 0.6, 'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'figure.facecolor': BACKGROUND, 'savefig.facecolor': BACKGROUND, 'svg.fonttype': 'none', 'pdf.fonttype': 42})

def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in (MAPPING_PATH, CONTACT_PATH, CE_PATH, HIC_STRUCTURE, FOCUT5A_STRUCTURE):
        if not path.is_file():
            raise FileNotFoundError(path)
    mapping = pd.read_csv(MAPPING_PATH)
    contact = pd.read_csv(CONTACT_PATH)
    ce = pd.read_csv(CE_PATH)
    focus = mapping[(mapping['target_protein'] == 'FoCut5a') & mapping['hic_residue'].isin([pair[0] for pair in FOCUS_PAIRS])].copy()
    observed = {(int(row.hic_residue), str(row.hic_label), int(row.consensus_simulation_residue), str(row.consensus_label), str(row.mapping_class)) for row in focus.itertuples(index=False)}
    expected = {(hic_res, hic_label, fo_res, fo_label, 'high_confidence') for hic_res, hic_label, fo_res, fo_label in FOCUS_PAIRS}
    if observed != expected:
        raise ValueError(f'Unexpected HiC--FoCut5a mapping: {observed}')
    contact_focus = contact[contact['protein'].isin(['HiC', 'FoCut5a']) & contact['pet_kind'].isin(['PET_L10', 'PET_L20']) & contact['hic_residue'].isin([pair[0] for pair in FOCUS_PAIRS])].copy()
    if len(contact_focus) != 24:
        raise ValueError(f'Expected 24 focused contact rows, found {len(contact_focus)}')
    expected_n = {('HiC', 'PET_L10'): 5, ('HiC', 'PET_L20'): 8, ('FoCut5a', 'PET_L10'): 5, ('FoCut5a', 'PET_L20'): 7}
    observed_n = {(protein, pet_kind): int(group['n_trajectories'].iloc[0]) for (protein, pet_kind), group in contact_focus.groupby(['protein', 'pet_kind'])}
    if observed_n != expected_n:
        raise ValueError(f'Unexpected focused retained denominators: {observed_n}')
    return (focus, contact_focus, ce)

def extract_anchors(path: Path) -> list[tuple[float, float]]:
    image = np.asarray(Image.open(path).convert('RGB'), dtype=float) / 255.0
    hsv = mcolors.rgb_to_hsv(image)
    anchors: list[tuple[float, float]] = []
    for marker_color in ANCHOR_COLORS:
        target_rgb = np.asarray(mcolors.to_rgb(marker_color), dtype=float)
        target_hsv = mcolors.rgb_to_hsv(target_rgb.reshape(1, 1, 3))[0, 0]
        hue_distance = np.abs(hsv[:, :, 0] - target_hsv[0])
        hue_distance = np.minimum(hue_distance, 1.0 - hue_distance)
        mask = (hue_distance <= 0.035) & (hsv[:, :, 1] >= 0.45) & (hsv[:, :, 2] >= 0.18)
        rows, columns = np.where(mask)
        if len(rows) < 20:
            raise ValueError(f'Could not recover Figure 5 anchor {marker_color}')
        anchors.append((float(np.median(columns) / (image.shape[1] - 1)), float(1.0 - np.median(rows) / (image.shape[0] - 1))))
    return anchors

def paired_residue_label(hic_text: str, focut5a_text: str, *, multiline: bool=False) -> HPacker | VPacker:
    """Build a mixed-color HiC/FoCut5a structural-pair label."""
    common = {'fontsize': 6.0}
    hic_area = TextArea(hic_text, textprops={**common, 'color': HIC_DARK})
    arrow_area = TextArea(' ↔ ', textprops={**common, 'color': MUTED})
    focut5a_area = TextArea(focut5a_text, textprops={**common, 'color': FOCUT5A_DARK})
    if not multiline:
        return HPacker(children=[hic_area, arrow_area, focut5a_area], align='center', pad=0, sep=0)
    top_line = HPacker(children=[hic_area, arrow_area], align='center', pad=0, sep=0)
    return VPacker(children=[top_line, focut5a_area], align='right', pad=0, sep=0)

def crop_contract() -> tuple[float, float, float, float]:
    crop_width = 1.0 / STRUCTURE_ZOOM
    crop_height = 1.0 / STRUCTURE_ZOOM
    crop_left = (1.0 - crop_width) / 2.0
    centered_bottom = (1.0 - crop_height) / 2.0
    crop_bottom = centered_bottom + STRUCTURE_DOWN_SHIFT * crop_height
    crop_top = 1.0 - crop_bottom - crop_height
    return (crop_left, crop_top, crop_width, crop_height)

def transform_position(position: tuple[float, float]) -> tuple[float, float]:
    crop_left, crop_top, crop_width, crop_height = crop_contract()
    crop_bottom = 1.0 - crop_top - crop_height
    x, y = position
    return ((x - crop_left) / crop_width, (y - crop_bottom) / crop_height)

def crop_structure(image: Image.Image) -> Image.Image:
    crop_left, crop_top, crop_width, crop_height = crop_contract()
    width, height = image.size
    left = round(crop_left * width)
    top = round(crop_top * height)
    right = round((crop_left + crop_width) * width)
    bottom = round((crop_top + crop_height) * height)
    return image.crop((left, top, right, bottom))

def add_halo_pair_label(ax: plt.Axes, hic_text: str, focut5a_text: str, position: tuple[float, float]) -> None:
    label = paired_residue_label(hic_text, focut5a_text)
    for text_area in label.get_children():
        text_area._text.set_path_effects([path_effects.withStroke(linewidth=1.05, foreground='white')])
    ax.add_artist(AnnotationBbox(label, xy=position, xycoords=ax.transAxes, box_alignment=(0.5, 0.5), frameon=False, pad=0, annotation_clip=False, zorder=16))

def draw_structural_panel(ax: plt.Axes, close_path: Path, anchor_path: Path) -> mpl.text.Text:
    del anchor_path
    ax.set_axis_off()
    image_ax = ax.inset_axes([0.0, 0.055, 1.0, 0.945])
    structure = Image.open(close_path).convert('RGB')
    image_ax.imshow(crop_structure(structure))
    image_ax.set_axis_off()
    for hic_text, focut5a_text, position in STRUCTURAL_LABELS:
        add_halo_pair_label(image_ax, hic_text, focut5a_text, transform_position(position))
    header = panel_label(ax, 'a')
    handles = [mlines.Line2D([], [], color=HIC_COLOR, linewidth=3.2, label='HiC'), mlines.Line2D([], [], color=FOCUT5A_COLOR, linewidth=3.2, label='FoCut5a')]
    ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, -0.025), ncol=2, frameon=False, handlelength=1.6, columnspacing=1.2)
    return header

def draw_contact_panel(ax: plt.Axes, contact: pd.DataFrame) -> mpl.text.Text:
    columns = [('HiC', 'PET_L10'), ('HiC', 'PET_L20'), ('FoCut5a', 'PET_L10'), ('FoCut5a', 'PET_L20')]
    matrix = np.zeros((len(FOCUS_PAIRS), len(columns)), dtype=float)
    for row_index, (hic_res, _, _, _) in enumerate(FOCUS_PAIRS):
        for column_index, (protein, pet_kind) in enumerate(columns):
            match = contact[(contact['hic_residue'] == hic_res) & (contact['protein'] == protein) & (contact['pet_kind'] == pet_kind)]
            if len(match) != 1:
                raise ValueError(f'Missing mapped-site contact: {hic_res}, {protein}, {pet_kind}')
            matrix[row_index, column_index] = float(match.iloc[0]['mean_contact_fraction'])
    cmap = mcolors.LinearSegmentedColormap.from_list('contact_use', ['#F6F1F8', '#B9A7C7', '#604572'])
    image = ax.imshow(matrix, vmin=0.0, vmax=1.0, cmap=cmap, aspect='auto')
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            ax.text(column_index, row_index, f'{value:.2f}', ha='center', va='center', fontsize=6.0, fontweight='bold', color=contrasting_text_color(image.cmap(image.norm(value)), dark=INK))
    ax.axvline(1.5, color=INK, linewidth=1.1)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels(['L10\n(n=5)', 'L20\n(n=8)', 'L10\n(n=5)', 'L20\n(n=7)'], rotation=0, ha='center')
    for tick in ax.get_xticklabels():
        tick.set_color(INK)
    axis_transform = ax.get_xaxis_transform()
    ax.text(0.5, 1.035, 'HiC', transform=axis_transform, ha='center', va='bottom', fontsize=6.5, fontweight='bold', color=HIC_DARK, clip_on=False)
    ax.text(2.5, 1.035, 'FoCut5a', transform=axis_transform, ha='center', va='bottom', fontsize=6.5, fontweight='bold', color=FOCUT5A_DARK, clip_on=False)
    ax.set_yticks(range(len(FOCUS_PAIRS)))
    native_pair_labels = [f'{hic_label} ↔ {fo_label}' for _, hic_label, _, fo_label in FOCUS_PAIRS]
    ax.set_yticklabels(native_pair_labels)
    for tick in ax.get_yticklabels():
        tick.set_alpha(0.0)
    for row_index, (_, hic_label, _, fo_label) in enumerate(FOCUS_PAIRS):
        ax.add_artist(AnnotationBbox(paired_residue_label(hic_label, fo_label), xy=(-0.55, row_index), xycoords=ax.transData, box_alignment=(1.0, 0.5), frameon=False, pad=0, annotation_clip=False))
    ax.tick_params(axis='x', length=0, pad=3)
    ax.tick_params(axis='y', length=0, pad=2)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(FOCUS_PAIRS), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.75)
    ax.tick_params(which='minor', bottom=False, left=False)
    header = panel_label(ax, 'b')
    cbar = ax.figure.colorbar(image, ax=ax, orientation='vertical', fraction=0.055, pad=0.04, ticks=[0.0, 0.5, 1.0], use_gridspec=False)
    cbar.ax.set_box_aspect(None)
    cbar.ax.set_aspect('auto')

    def matched_bounds(colorbar_axes, renderer):
        bounds = ax.get_position()
        return Bbox.from_bounds(bounds.x1 + bounds.width * (0.04 / 0.905), bounds.y0, bounds.width * (0.055 / 0.905), bounds.height)
    cbar.ax.set_axes_locator(matched_bounds)
    cbar.set_label('Mean contact fraction', labelpad=3)
    cbar.outline.set_linewidth(0.5)
    return header

def build_figure(rendered) -> None:
    apply_style()
    mapping, contact, ce = load_inputs()
    if not mapping.mapping_class.eq('high_confidence').all():
        raise ValueError('Focused mappings must be high confidence')
    close_path = rendered / 'Figure_06_HiC_FoCut5a_overlay.png'
    anchor_path = rendered / 'Figure_06_HiC_FoCut5a_anchor_map.png'
    fig = plt.figure(figsize=(7.0, 3.85))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.05, 0.95], left=0.035, right=0.915, top=0.91, bottom=0.12, wspace=0.3)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    heatmap_box = ax_b.get_position()
    heatmap_height = heatmap_box.height * HEATMAP_HEIGHT_FRACTION
    ax_b.set_position([heatmap_box.x0, heatmap_box.y1 - heatmap_height, heatmap_box.width, heatmap_height])
    panel_headers = [(ax_a, draw_structural_panel(ax_a, close_path, anchor_path)), (ax_b, draw_contact_panel(ax_b, contact))]
    fig._publication_panel_legend_layout = {'key': 'Figure_05', 'panel': ax_a, 'legend_gap_pt': 10.0}
    align_panel_labels_to_yticklabels(fig, panel_headers)
    return fig
MAPPING_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/hic_structural_equivalent_mapping.csv'
CONTACT_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/hic_equivalent_contact_fraction_by_length.csv'
CE_PATH = ROOT / '02_postprocessing/results/aggregate_statistics/hic_pairwise_ce_alignment_summary.csv'
HIC_STRUCTURE = ROOT / '01_simulation/inputs/production_structures/00121.pdb'
FOCUT5A_STRUCTURE = ROOT / '01_simulation/inputs/production_structures/00075.pdb'
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_05', build_figure)
