"""Run explicitly requested figures using package-local code and data."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figures',nargs='+',required=True,help='Exact IDs, for example Figure_02 Figure_S01, or all')
    parser.add_argument('--output-dir',type=Path,required=True,help='Write PNGs to a new directory outside the package')
    parser.add_argument('--chimerax',type=Path,help='Explicit tested executable for structural figures')
    args=parser.parse_args()
    # Resolve user paths before subprocesses switch to the package root.
    if args.output_dir is not None:
        args.output_dir=args.output_dir.resolve()
    if args.output_dir==ROOT or ROOT in args.output_dir.parents or args.output_dir.exists():
        parser.error('Use a new output directory outside the package; reference PNGs are read-only')
    if args.chimerax is not None:
        args.chimerax=args.chimerax.resolve()
    registry=json.loads((ROOT/'metadata/figures.json').read_text())
    selected=args.figures
    if selected==['all']:
        selected=[f'Figure_{i:02d}' for i in range(1,7)]+[f'Figure_S{i:02d}' for i in range(1,10)]
    missing=[key for key in selected if key not in registry]
    if missing:parser.error('Unknown figure ID: '+', '.join(missing))
    if len(set(selected)) != len(selected):parser.error('Duplicate figure IDs are not allowed')
    if any('render_scenes' in registry[key] for key in selected) and args.chimerax is None:
        parser.error('Structural figures require an explicit --chimerax executable')
    # Reject output conflicts before producing any part of a requested batch.
    for key in selected:
        record=registry[key]
        output=(args.output_dir/Path(record['output']).name
                if args.output_dir is not None else ROOT/record['output'])
        if output.exists():
            parser.error('Output exists; use a new directory: '+str(output))
    for key in selected:
        record=registry[key]
        command=[sys.executable,'-B',str(ROOT/record['entrypoint'])]
        if args.output_dir is not None:
            command+=['--output',str(args.output_dir/Path(record['output']).name)]
        if args.chimerax is not None:command+=['--chimerax',str(args.chimerax)]
        subprocess.run(command,cwd=ROOT,check=True)


if __name__=='__main__':main()
