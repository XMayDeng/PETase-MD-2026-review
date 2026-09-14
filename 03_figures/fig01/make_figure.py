from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'shared'))
import json
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as mpatheffects
import numpy as np
from PIL import Image
ROOT = Path(__file__).resolve().parents[2]
PROTEINS = [('IsPETase', '6EQD'), ('TfCut1', '7QJR'), ('LCC', '4EB0'), ('FoCut5a', '5AJH'), ('HiC', '4OYY'), ('PHL7', '7NEI')]
LABEL_LAYOUT = {'IsPETase': {'S160': (0.46, 0.68), 'D206': (0.65, 0.18), 'H237': (0.12, 0.37), 'W185': (0.65, 0.9)}, 'TfCut1': {'S131': (0.4, 0.7), 'D177': (0.6, 0.24), 'H209': (0.16, 0.4), 'W156': (0.84, 0.78)}, 'LCC': {'S165': (0.42, 0.77), 'D210': (0.62, 0.31), 'H242': (0.16, 0.44), 'W190': (0.84, 0.86)}, 'FoCut5a': {'S122': (0.56, 0.8), 'D177': (0.64, 0.28), 'H190': (0.14, 0.5)}, 'HiC': {'S105': (0.55, 0.8), 'D160': (0.62, 0.29), 'H173': (0.14, 0.5)}, 'PHL7': {'S131': (0.4, 0.73), 'D177': (0.6, 0.26), 'H209': (0.16, 0.41), 'W156': (0.84, 0.84)}}

def crop_context(path: Path) -> Image.Image:
    """Apply the original whole-fold whitespace-only assembly crop."""
    picture = Image.open(path).convert('RGB')
    rows, columns = np.where(np.any(np.asarray(picture) < 248, axis=2))
    padding = int(max(np.ptp(rows), np.ptp(columns)) * 0.07)
    return picture.crop((max(0, int(columns.min()) - padding), max(0, int(rows.min()) - padding), min(picture.width, int(columns.max()) + padding + 1), min(picture.height, int(rows.max()) + padding + 1)))

def build_figure(rendered) -> None:
    WORK = rendered
    ARCHIVE = rendered
    plt.rcParams.update({'font.family': ['DejaVu Sans'], 'font.sans-serif': ['DejaVu Sans'], 'font.size': 7.0, 'text.color': '#20242B', 'figure.facecolor': 'white'})
    fig, axes = plt.subplots(2, 3, figsize=(6.95, 4.35))
    callout_panels = []
    for index, (ax, (name, pdb)) in enumerate(zip(axes.flat, PROTEINS)):
        ax.set_axis_off()
        context = ax.inset_axes([0.005, 0.1, 0.385, 0.76])
        context.imshow(crop_context(ARCHIVE / f'Figure_01_{name}_context.png'))
        context.set_axis_off()
        pocket = ax.inset_axes([0.39, 0.095, 0.61, 0.8])
        pocket.imshow(Image.open(WORK / f'Figure_01_{name}_active_site.png').convert('RGB'))
        pocket.set_axis_off()
        projection = json.loads((WORK / f'Figure_01_{name}_label_projection.json').read_text())
        labels = []
        for record in projection['records']:
            label = pocket.annotate(record['label'], xy=record['anchor_screen_fraction'], xytext=LABEL_LAYOUT[name][record['label']], xycoords='axes fraction', textcoords='axes fraction', ha='center', va='center', fontsize=7.0, color='#161A20', arrowprops={'arrowstyle': '-', 'color': '#626A75', 'linewidth': 0.6, 'shrinkA': 3.0, 'shrinkB': 2.0}, annotation_clip=False)
            label.set_path_effects([mpatheffects.withStroke(linewidth=1.4, foreground='white')])
            labels.append(label)
        callout_panels.append((name, pocket, projection, labels))
        ax.text(0.008, 0.99, f'({chr(97 + index)})', transform=ax.transAxes, ha='left', va='top', fontsize=8.0, fontweight='bold')
        ax.text(0.53, 0.985, f'{name} ({pdb})', transform=ax.transAxes, ha='center', va='top', fontsize=7.5, fontweight='bold')
        for x, label in [(0.1975, 'Whole fold'), (0.695, 'Active-site detail')]:
            ax.text(x, 0.065, label, transform=ax.transAxes, ha='center', va='top', fontsize=6.3, color='#626A75')
    fig.subplots_adjust(left=0.012, right=0.988, top=0.995, bottom=0.12, wspace=0.018, hspace=0.045)
    fig.legend(handles=[mpatches.Patch(facecolor='#E07B54', label='Catalytic triad'), mpatches.Patch(facecolor='#E8C547', label='Assigned W-loop Trp')], loc='lower center', bbox_to_anchor=(0.5, 0.018), ncol=2, frameon=False, fontsize=7.0, handlelength=1.1, columnspacing=1.5)
    return fig
if __name__ == '__main__':
    from runtime import run_figure
    run_figure('Figure_01', build_figure)
