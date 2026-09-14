#!/usr/bin/env python3
"""Recreate all 15 figures and six table environments in a new output folder."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent.parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--chimerax',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve()
    if out.exists() or out==ROOT or ROOT in out.parents:
        parser.error('Use a new output directory outside the package')
    out.mkdir(parents=True)
    commands=[
        ['reproduce.py','--figures','all','--output-dir',str(out/'figures'),'--chimerax',str(args.chimerax.resolve())],
        ['02_postprocessing/code/tables/build_tables.py','--output-dir',str(out/'tables')],
        ['02_postprocessing/code/tables/build_supplementary_table_s5.py','--output',
         str(out/'tables/Supplementary_Table_S5_per_ring_contact_geometry.tex')],
    ]
    checks=[]
    for i,args in enumerate(commands,1):
        with (out/f'step_{i:02d}.log').open('w') as log:
            result=subprocess.run([sys.executable,'-B',str(ROOT/args[0]),*args[1:]],
                                  stdout=log,stderr=subprocess.STDOUT)
        checks.append(dict(step=i,entry=args[0],status='PASS' if result.returncode==0 else 'FAIL'))
        (out/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        print(json.dumps(checks[-1]),flush=True)
        if result.returncode:raise SystemExit(result.returncode)
    assert len(list((out/'figures').glob('*.png')))==15
    assert len(list((out/'tables').glob('*.tex')))==6


if __name__=='__main__':main()
