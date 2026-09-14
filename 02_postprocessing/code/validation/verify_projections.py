"""Verify three published projections using only this working copy's inputs.

Run with python -B 02_postprocessing/code/validation/verify_projections.py from the package root.
Only a temporary S1 CSV is written; expected data files remain read-only.
"""
import hashlib
import importlib.util
import json
import platform
import tempfile
from pathlib import Path

import pandas

ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / '02_postprocessing/results/aggregate_statistics'
PUBLISHED = ROOT / '02_postprocessing/results/tables'


def module(name):
    path = ROOT / '02_postprocessing/code/tables' / f'build_supplementary_{name}.py'
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in ROOT.rglob('*') if p.is_file()}
    s1 = module('data_s1')
    with tempfile.TemporaryDirectory(prefix='petase_projection_check_') as temp:
        s1.OUTPUT = Path(temp) / 'Data_S1.csv'
        s1.main()
        assert s1.OUTPUT.read_bytes() == (PUBLISHED / 'Supplementary_Data_S1_residue_class_assignments.csv').read_bytes()
    s2 = module('data_s2')
    content, counts = s2.build(s2.SOURCE, s2.SUMMARY, s2.TABLE)
    assert content == s2.OUTPUT.read_bytes()
    s5 = module('table_s5')
    projected = s5.project_source(INPUTS / 'interaction_shift_statistics.csv')
    assert projected == (PUBLISHED / f'{s5.STEM}.csv').read_bytes()
    assert s5.render(projected).encode() == (PUBLISHED / f'{s5.STEM}.tex').read_bytes()
    after = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in ROOT.rglob('*') if p.is_file()}
    assert before == after, 'Verification changed working-copy files'
    print(json.dumps(dict(status='PASS', python=platform.python_version(), pandas=pandas.__version__,
                          Data_S1='1532 rows; byte-exact', Data_S2=counts,
                          Table_S5='CSV and TeX byte-exact', working_copy_unchanged=True), indent=2))


if __name__ == '__main__':
    main()
