#!/usr/bin/env python3
"""Redraw the accepted TOC schematic into a new external PNG."""
import argparse
from contextlib import redirect_stdout
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path,
                        help='New external file named TOC_graphic.png; no overwrite')
    parser.add_argument('--verbose', action='store_true', help='Show checksum and drawing diagnostics')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Use Python without -O/-OO; pixel-check assertions must stay enabled')
    output = args.output.resolve()
    protected_root = HERE.parents[1] if HERE.parent.name == '03_figures' else HERE
    if output.name != 'TOC_graphic.png':
        parser.error('Use the filename TOC_graphic.png')
    if args.output.is_symlink() or output.exists():
        parser.error('Output already exists; choose a new destination')
    if output == protected_root or protected_root in output.parents:
        parser.error('Output must be outside the package/source folder')
    contract = json.loads((HERE / 'source_manifest.json').read_text())
    for name, expected in contract['inputs'].items():
        if digest(HERE / name) != expected:
            raise RuntimeError('TOC input hash mismatch: ' + name)
    environment = contract['environment']
    if '.'.join(map(str, sys.version_info[:3])) != environment['python']:
        raise RuntimeError('Use the Python version recorded in source_manifest.json')
    for name, expected in environment['packages'].items():
        if importlib.metadata.version(name) != expected:
            raise RuntimeError('Dependency version mismatch: ' + name)
    sys.dont_write_bytecode = True
    with tempfile.TemporaryDirectory(prefix='petase_toc_render_') as temporary:
        scratch = Path(temporary)
        os.environ['MPLCONFIGDIR'] = str(scratch / 'matplotlib_cache')
        os.environ['MPLBACKEND'] = 'Agg'
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.font_manager as fm
        import numpy as np
        from PIL import Image
        import make_toc_graphic_schematic as base
        import make_toc_graphic_schematic_v1_refined as refined

        for name, expected in environment['fonts'].items():
            if digest(Path(mpl.get_data_path()) / 'fonts/ttf' / name) != expected:
                raise RuntimeError('Bundled font hash mismatch: ' + name)
        configure_original = base.configure_style

        def configure_exact_style():
            configure_original()
            mpl.rcParams['font.family'] = ['DejaVu Sans']
            mpl.rcParams['font.sans-serif'] = ['DejaVu Sans']

        base.configure_style = configure_exact_style
        resolved_font = Path(fm.findfont(fm.FontProperties(family=['DejaVu Sans']),
                                         fallback_to_default=False))
        required_font = Path(mpl.get_data_path()) / 'fonts/ttf/DejaVuSans.ttf'
        if resolved_font.resolve() != required_font.resolve():
            raise RuntimeError('Unexpected font resolution')
        staged = scratch / 'TOC_graphic.png'
        refined.OUTPUT = staged
        captured = io.StringIO()
        with redirect_stdout(captured):
            refined.main()
        drawing_checks = json.loads(captured.getvalue())
        drawing_checks.pop('output', None)
        image_hash = digest(staged)
        if image_hash != contract['reference_sha256']:
            raise RuntimeError('Redrawn PNG differs from the accepted TOC; no output published')
        with Image.open(staged) as picture:
            assert list(picture.size) == contract['reference_pixels']
            assert picture.mode in ('RGB', 'RGBA')
            if picture.mode == 'RGBA':
                assert np.all(np.asarray(picture)[:, :, 3] == 255)
        for name, expected in contract['inputs'].items():
            assert digest(HERE / name) == expected, 'Drawing modified its input: ' + name
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('xb') as handle:
            handle.write(staged.read_bytes())
        assert digest(output) == image_hash
    if args.verbose:
        print(json.dumps(dict(status='PASS', output=str(output), sha256=image_hash,
            reference_sha256=contract['reference_sha256'], exact_reference_match=True,
            drawing_checks=drawing_checks, full_md_trajectories_required=False,
            source_files_unchanged=True), indent=2))
    else:
        print(f'PASS: TOC graphic (975 x 525 pixels) -> {output}')


if __name__ == '__main__':
    main()
