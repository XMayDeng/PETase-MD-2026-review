"""Strict, package-local execution for registered manuscript figures."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import tempfile

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def registry():
    return json.loads((ROOT/'metadata/figures.json').read_text())


def check_environment():
    if sys.flags.optimize:
        raise RuntimeError('Optimized Python is not supported: validation assertions must remain enabled')
    contract = json.loads((ROOT/'environment/versions.json').read_text())
    actual_python = '.'.join(str(x) for x in sys.version_info[:3])
    if actual_python != contract['python']:
        raise RuntimeError(f"Python version mismatch: {actual_python} != {contract['python']}")
    for package, expected in contract['packages'].items():
        actual = importlib.metadata.version(package)
        if actual != expected:
            raise RuntimeError(f'{package} version mismatch: {actual} != {expected}')
    for face, expected in contract['fonts'].items():
        path = Path(mpl.get_data_path())/'fonts/ttf'/face
        if not path.is_file() or sha256(path) != expected:
            raise RuntimeError(f'Required bundled font missing or changed: {face}')
    mpl.rcParams['font.family'] = ['DejaVu Sans']
    mpl.rcParams['font.sans-serif'] = ['DejaVu Sans']
    resolved = Path(fm.findfont(fm.FontProperties(family=['DejaVu Sans']), fallback_to_default=False))
    required = Path(mpl.get_data_path())/'fonts/ttf/DejaVuSans.ttf'
    if resolved.resolve() != required.resolve():
        raise RuntimeError('Unexpected DejaVu Sans font resolution')


def check_inputs(figure_id):
    item = registry()[figure_id]
    for relative, expected in item['inputs'].items():
        path = ROOT/relative
        if not path.is_file():
            raise FileNotFoundError(f'{figure_id}: required input missing: {relative}')
        if sha256(path) != expected:
            raise RuntimeError(f'{figure_id}: input hash mismatch: {relative}')
    return item


def run_figure(figure_id, builder):
    parser = argparse.ArgumentParser(description=f'Reproduce {figure_id} from declared package inputs')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--replace', action='store_true', help='Replace this package-generated output only')
    parser.add_argument('--chimerax',type=Path,help='Required explicit executable for molecular-structure figures')
    parser.add_argument('--verbose', action='store_true', help='Show checksum and layout diagnostics')
    args = parser.parse_args()
    item = check_inputs(figure_id)
    output = args.output if args.output is not None else ROOT/item['output']
    if output.suffix.lower() != '.png':
        parser.error('Only PNG output is supported')
    if output.name != Path(item['output']).name:
        parser.error('Use the registered PNG filename: '+Path(item['output']).name)
    if output.exists() and not args.replace:
        parser.error('Output exists; choose a new output or explicitly pass --replace')
    check_environment()
    if 'render_scenes' in item and args.chimerax is None:
        parser.error('This figure requires --chimerax with the exact tested executable; no raster-cache substitution is allowed')
    with tempfile.TemporaryDirectory(prefix='petase_molecular_render_') as temporary:
        work=Path(temporary)
        if 'render_scenes' in item:
            from molecular_render import render_scenes
            render_scenes(ROOT,item,args.chimerax,work)
            fig=builder(work)
        else:
            fig = builder()
        from publication_style import save_publication_figure
        staged=work/output.name
        report = save_publication_figure(fig, staged)
        plt.close(fig)
        if report['text_outside_canvas'] or not report['evidence_preserved']:
            raise RuntimeError(f'{figure_id}: layout or numeric-artist QA failed')
        if sha256(staged) != item['reference_sha256']:
            raise RuntimeError(f'{figure_id}: generated PNG differs from the declared reference; inspect the cause')
        # Publish only a validated image. A mismatch leaves no misleading PNG.
        output.parent.mkdir(parents=True, exist_ok=True)
        mode='wb' if args.replace else 'xb'
        with output.open(mode) as handle:
            handle.write(staged.read_bytes())
    if args.verbose:
        print(json.dumps(dict(figure=figure_id,status='PASS',png_sha256=sha256(output),
                             numeric_artists_preserved=True,text_outside_canvas=0),indent=2))
    else:
        print(f'PASS: {figure_id} -> {output}')
