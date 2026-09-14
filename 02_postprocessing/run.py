#!/usr/bin/env python3
"""One entry: recompute supported statistics and rebuild all delivered tables."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

BASE=Path(__file__).resolve().parent
ROOT=BASE.parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();out=args.output_dir.resolve()
    if out.exists() or out==ROOT or ROOT in out.parents:p.error('Use a new output directory outside the package')
    out.mkdir(parents=True)
    commands=[
        ['code/recompute_analysis.py','--output-dir',str(out/'statistics')],
        ['code/statistics/recompute_structural_sensitivity.py','--output-dir',str(out/'structural_sensitivity')],
        ['code/validation/verify_retention.py'],
        ['code/validation/verify_statistics.py'],
        ['code/tables/build_tables.py','--output-dir',str(out/'tables')],
        ['code/tables/build_supplementary_data_s1.py','--output',str(out/'Supplementary_Data_S1.csv')],
        ['code/tables/build_supplementary_data_s2.py','--output',str(out/'Supplementary_Data_S2.csv')],
        ['code/tables/build_supplementary_table_s5.py','--output',str(out/'tables/Supplementary_Table_S5_per_ring_contact_geometry.tex')],
    ]
    completed=[]
    for index,command in enumerate(commands,1):
        cmd=[sys.executable,'-B',str(BASE/command[0]),*command[1:]]
        with (out/f'step_{index:02d}.log').open('w') as log:
            result=subprocess.run(cmd,cwd=out,stdout=log,stderr=subprocess.STDOUT)
        completed.append(dict(step=index,entry=command[0],status='PASS' if result.returncode==0 else 'FAIL'))
        (out/'checks.json').write_text(json.dumps(completed,indent=2)+'\n')
        print(json.dumps(completed[-1]),flush=True)
        if result.returncode:raise SystemExit(result.returncode)


if __name__=='__main__':main()
