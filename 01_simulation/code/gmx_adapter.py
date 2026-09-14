#!/usr/bin/env python3
"""Transparent GROMACS runtime adapter used by run.py; records every override."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys

args=sys.argv[1:]
original=list(args)
mode=os.environ['PETASE_RUN_MODE']
if args and args[0]=='mdrun':
    phase=args[args.index('-deffnm')+1] if '-deffnm' in args else ''
    if mode=='smoke' and Path(phase).name!='EM':
        if '-nsteps' in args:
            i=args.index('-nsteps');del args[i:i+2]
        args+=['-nsteps','200']
    if '-ntomp' not in args:args+=['-ntomp',os.environ['PETASE_THREADS']]
    args+=['-nb','gpu' if os.environ.get('PETASE_GPU') else 'cpu','-pme','cpu','-bonded','cpu','-update','cpu']
if args and args[0]=='grompp':
    if '-maxwarn' in args:args[args.index('-maxwarn')+1]='0'
    else:args+=['-maxwarn','0']

def record(**extra):
    with Path(os.environ['PETASE_COMMAND_LOG']).open('a') as f:
        f.write(json.dumps(dict(original=original,executed=args,mode=mode,**extra))+'\n')

record()
real=os.environ['PETASE_REAL_GMX']
if args and args[0]=='grompp':
    proc=subprocess.run([real,*args],capture_output=True,text=True)
    sys.stdout.write(proc.stdout);sys.stderr.write(proc.stderr)
    warnings=re.findall(r'WARNING \d+ \[.*?\]:\n(.*?)(?=\n\n)',proc.stderr,re.S)
    output=args[args.index('-o')+1] if '-o' in args else ''
    if proc.returncode and len(warnings)==1 and 'The Berendsen barostat' in warnings[0] and Path(output).name=='NPT.tpr':
        args[args.index('-maxwarn')+1]='1'
        record(reviewed_warning=warnings[0],scope='Historical equilibration warning only; all other warnings block')
        raise SystemExit(subprocess.call([real,*args]))
    raise SystemExit(proc.returncode)
raise SystemExit(subprocess.call([real,*args]))
