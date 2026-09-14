#!/usr/bin/env python3
"""List, inspect or replay the 126 delivered production MD inputs.

Examples from the repository root:
  python -B 01_simulation/run.py list
  python -B 01_simulation/run.py check --output-dir ../tpr_check
  python -B 01_simulation/run.py replay --name TfCut1_L10_bi/rep01 --steps 100 --output-dir ../md_smoke

Only production MD TPRs are included. Full upstream preparation is retained as
documented source code for inspection, not as this review entry's execution mode.
"""
from pathlib import Path
import subprocess
import sys


if __name__ == '__main__':
    entry = Path(__file__).resolve().parent / 'code/run_simulation.py'
    raise SystemExit(subprocess.run([sys.executable, '-B', str(entry), *sys.argv[1:]]).returncode)
