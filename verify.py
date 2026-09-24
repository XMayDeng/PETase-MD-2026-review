"""Verify the review repository and its declared numerical reproductions."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
IGNORED_DIRS = {'.git', '__pycache__', '.pytest_cache', '.venv'}
LOCAL_OUTPUTS = {('work',), ('validation_runs',)}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inventory():
    with (ROOT / 'metadata/file_manifest.csv').open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    expected = {row['path'] for row in rows}
    if len(expected) != len(rows):
        raise ValueError('Duplicate manifest path')
    actual = set()
    for path in ROOT.rglob('*'):
        parts = path.relative_to(ROOT).parts
        if any(part in IGNORED_DIRS for part in parts):
            continue
        if any(parts[:len(prefix)] == prefix for prefix in LOCAL_OUTPUTS):
            continue
        if path.is_symlink():
            raise ValueError('Symlink not allowed in the review package: ' + str(path))
        if path.is_file() and path.relative_to(ROOT).as_posix() != 'metadata/file_manifest.csv':
            actual.add(path.relative_to(ROOT).as_posix())
    if actual != expected:
        raise ValueError(f'File-set mismatch: missing={sorted(expected-actual)}, extra={sorted(actual-expected)}')
    for row in rows:
        relative = Path(row['path'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Unsafe manifest path')
        path = ROOT / relative
        if path.stat().st_nlink != 1 or path.stat().st_size != int(row['bytes']) or digest(path) != row['sha256']:
            raise ValueError('Manifest mismatch: ' + row['path'])
    with (ROOT / '01_simulation/configurations/stages.csv').open(newline='') as handle:
        stages = list(csv.DictReader(handle))
    stage_paths = {r['tpr_path'] for r in stages}
    if len(stages) != 126 or len(stage_paths) != 126 or any(r['stage'] != 'MD' for r in stages):
        raise ValueError('Expected exactly 126 production MD stage records')
    if stage_paths != {row['path'] for row in rows if row['path'].endswith('.tpr')}:
        raise ValueError('Stage index and delivered TPR files disagree')
    indexed = {r['path']: r for r in rows}
    for stage in stages:
        if stage['sha256'] != indexed[stage['tpr_path']]['sha256']:
            raise ValueError('TPR stage/manifest hash mismatch')
    for registry_name in ('analyses.json', 'figures.json', 'tables.json'):
        registry = json.loads((ROOT / 'metadata' / registry_name).read_text())
        for key, record in registry.items():
            for field in ('inputs', 'references'):
                values = record.get(field, {})
                if not isinstance(values, dict):
                    continue
                for relative, expected_hash in values.items():
                    path = ROOT / relative
                    if not path.is_file() or digest(path) != expected_hash:
                        raise ValueError(f'{key}: required {field} mismatch: {relative}')
    toc = json.loads((ROOT / '03_figures/toc/source_manifest.json').read_text())
    for name, expected_hash in toc['inputs'].items():
        if digest(ROOT / '03_figures/toc' / name) != expected_hash:
            raise ValueError('TOC input mismatch: ' + name)
    if digest(ROOT / '03_figures/toc/TOC_graphic.png') != toc['reference_sha256']:
        raise ValueError('TOC reference image mismatch')
    return rows


def check_generated_tables(output):
    registry = json.loads((ROOT / 'metadata/tables.json').read_text())
    for key, record in registry.items():
        generated = output / 'tables' / Path(record['output']).name
        if not generated.is_file() or digest(generated) != record['reference_sha256']:
            raise ValueError('Regenerated table differs: ' + key)
    pairs = {
        'tables/Supplementary_Table_S5_per_ring_contact_geometry.tex':
            'Supplementary_Table_S5_per_ring_contact_geometry.tex',
        'Supplementary_Data_S1.csv': 'Supplementary_Data_S1_residue_class_assignments.csv',
        'Supplementary_Data_S2.csv': 'Supplementary_Data_S2_residue_contact_shift_candidates.csv',
    }
    for generated, reference in pairs.items():
        path = output / generated
        if not path.is_file() or digest(path) != digest(ROOT / '02_postprocessing/results/tables' / reference):
            raise ValueError('Regenerated supplementary artifact differs: ' + generated)


def run_check(label, command, *, verbose=False, cwd=None):
    print(f'Checking {label}...', flush=True)
    result = subprocess.run(command, cwd=cwd, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if verbose or result.returncode:
        print(result.stdout, end='', flush=True)
    result.check_returncode()
    print(f'PASS: {label}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checksums-only', action='store_true', help='Check all files and declared input dependencies without running analyses')
    parser.add_argument('--gmx', help='Also parse all 126 TPRs with this GROMACS executable')
    parser.add_argument('--render', action='store_true', help='Also redraw all 15 numbered figures')
    parser.add_argument('--chimerax', type=Path)
    parser.add_argument('--output-dir', type=Path, help='New external figure output directory with --render')
    parser.add_argument('--verbose', action='store_true', help='Show detailed output from each check')
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error('Do not use Python -O/-OO')
    if args.checksums_only and (args.gmx or args.render or args.chimerax or args.output_dir):
        parser.error('--checksums-only cannot be combined with execution flags')
    if args.render != bool(args.chimerax and args.output_dir):
        parser.error('--render requires both --chimerax and --output-dir')
    if not args.render and (args.chimerax or args.output_dir):
        parser.error('Rendering paths require --render')
    rows = inventory()
    if not args.checksums_only:
        sys.path.insert(0, str(ROOT / '03_figures/shared'))
        from runtime import check_environment, check_inputs
        check_environment()
        figures = json.loads((ROOT / 'metadata/figures.json').read_text())
        for key, record in figures.items():
            check_inputs(key)
            if digest(ROOT / record['output']) != record['reference_sha256']:
                raise ValueError('Reference PNG mismatch: ' + key)
        commands = [
            ('preparation metadata', '01_simulation/verify_preparation.py'),
            ('supplementary-data exports', '02_postprocessing/code/validation/verify_projections.py'),
        ]
        for label, command in commands:
            run_check(label, [sys.executable, '-B', str(ROOT / command)], verbose=args.verbose)
        with tempfile.TemporaryDirectory(prefix='petase_review_check_') as temporary:
            work = Path(temporary)
            run_check('statistics and tables',
                      [sys.executable, '-B', str(ROOT / '02_postprocessing/run.py'),
                       '--output-dir', str(work / 'analysis')], verbose=args.verbose, cwd=work)
            check_generated_tables(work / 'analysis')
            if args.gmx:
                run_check('126 production TPRs',
                          [sys.executable, '-B', str(ROOT / '01_simulation/run.py'),
                           'check', '--gmx', args.gmx, '--output-dir', str(work / 'tpr_check')],
                          verbose=args.verbose, cwd=work)
        if args.render:
            run_check('15 numbered figures',
                      [sys.executable, '-B', str(ROOT / 'reproduce.py'), '--figures', 'all',
                       '--chimerax', str(args.chimerax.resolve()), '--output-dir', str(args.output_dir.resolve())],
                      verbose=args.verbose)
        # Supported execution must not change any delivered source or data.
        for row in rows:
            if digest(ROOT / row['path']) != row['sha256']:
                raise ValueError('Execution changed a delivered file: ' + row['path'])
    print(json.dumps(dict(status='PASS', manifest_files=len(rows), production_MD_TPRs=126,
                         numerical_recalculations=not args.checksums_only,
                         parsed_TPRs=126 if args.gmx else 0,
                         rendered_figures=15 if args.render else 0,
                         scope='Delivered-data statistics and artifact reconstruction; not recovery of original trajectories'), indent=2))


if __name__ == '__main__':
    main()
