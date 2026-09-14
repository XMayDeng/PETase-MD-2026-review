"""Read-only validation of selected statistics from bundled trajectory summaries.

Run: python -B 02_postprocessing/code/validation/verify_statistics.py
This does not read trajectories, select candidates, or overwrite provided data.
"""
import hashlib
import importlib.util
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[3]
INPUTS=ROOT/'02_postprocessing/results/aggregate_statistics'
DERIVED=ROOT/'02_postprocessing/results/trajectory_summaries'


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/f'02_postprocessing/code/statistics/{name}.py')
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


def compare(observed,expected):
    if isinstance(observed,(str,bool,np.bool_)):
        assert observed==expected,(observed,expected)
        return 0.
    np.testing.assert_allclose(float(observed),float(expected),rtol=0,atol=1e-12,equal_nan=True)
    return abs(float(observed)-float(expected)) if np.isfinite(observed) else 0.


def main():
    before={str(p.relative_to(ROOT)):digest(p) for p in ROOT.rglob('*') if p.is_file()}
    for r in json.loads((ROOT/'02_postprocessing/code/statistics/source_manifest.json').read_text()):
        assert before[r['release_path']]==r['copy_sha256']
    con,loc=module('con03'),module('loc01')
    expected=pd.read_csv(INPUTS/'trajectory_level_contact_shift_candidates.csv')
    contacts=pd.read_csv(DERIVED/'candidate_case_contact_values.csv')
    assert len(expected)==95 and not contacts.duplicated(['case','residue']).any()
    con_cells=0;largest=0.
    for row in expected.to_dict('records'):
        group=contacts[(contacts.protein_id==row['protein_id'])&(contacts.residue==row['residue'])]
        x=group[group.pet_kind=='PET_L10'].sort_values('case').contact_fraction.to_numpy(float)
        y=group[group.pet_kind=='PET_L20'].sort_values('case').contact_fraction.to_numpy(float)
        rng=np.random.default_rng(con.stable_seed(row['protein_id'],int(row['residue'])))
        bx=con.bootstrap_means(x,rng);by=con.bootstrap_means(y,rng)
        dx=con.ci95(bx);dy=con.ci95(by);delta=float(y.mean()-x.mean())
        interval=con.ci95(by-bx);positive=float(np.mean(by-bx>0))
        loo=con.leave_one_out_deltas(x,y)
        checked=dict(n_l10=len(x),n_l20=len(y),underpowered=min(len(x),len(y))<3,
                     mean_contact_l10=x.mean(),mean_contact_l20=y.mean(),
                     median_contact_l10=np.median(x),median_contact_l20=np.median(y),
                     sd_contact_l10=np.std(x,ddof=1) if len(x)>1 else float('nan'),
                     sd_contact_l20=np.std(y,ddof=1) if len(y)>1 else float('nan'),
                     mean_contact_l10_ci95_low=dx[0],mean_contact_l10_ci95_high=dx[1],
                     mean_contact_l20_ci95_low=dy[0],mean_contact_l20_ci95_high=dy[1],
                     mean_delta_l20_minus_l10=delta,delta_ci95_low=interval[0],delta_ci95_high=interval[1],
                     bootstrap_probability_delta_positive=positive,
                     bootstrap_direction_probability=positive if delta>=0 else 1-positive,
                     cliffs_delta_l20_vs_l10=con.cliffs_delta(y,x),
                     loo_sign_agreement=con.sign_agreement(loo,delta),
                     loo_min_delta=min(loo) if loo else float('nan'),loo_max_delta=max(loo) if loo else float('nan'))
        for label,values in [('l10',x),('l20',y)]:
            checked[f'n_cases_ge_0_2_{label}']=int(np.sum(values>=.2))
            checked[f'n_cases_ge_0_5_{label}']=int(np.sum(values>=.5))
            checked[f'support_ge_0_5_fraction_{label}']=float(np.mean(values>=.5))
        direction_names=['n_common_directions','common_directions','direction_delta_values',
                         'direction_adjusted_delta','direction_sign_agreement','direction_delta_range',
                         'n_direction_strata_effect_consistent','direction_effect_consistent_fraction']
        checked.update(zip(direction_names,con.direction_metrics(group,delta)))
        for key,value in checked.items():
            largest=max(largest,compare(value,row[key]));con_cells+=1
        # The all-atom sensitivity sign remains an input from the provided summary.
        shifted,status,core=con.classify_shift({**row,**checked})
        assert (shifted,status,core)==(row['shift_class'],row['manuscript_status'],row['stable_contact_core'])
    s5=pd.read_csv(ROOT/'02_postprocessing/results/tables/Supplementary_Table_S5_per_ring_contact_geometry.csv').set_index('metric_id')
    metrics=pd.read_csv(DERIVED/'trajectory_metric_values.csv')
    selected=metrics[metrics.metric_id.isin(s5.index)]
    rebuilt=loc.build_shift_statistics(selected).set_index('metric_id')
    loc_cells=0
    for metric_id,row in s5.iterrows():
        for column,expected_value in row.items():
            largest=max(largest,compare(rebuilt.loc[metric_id,column],expected_value));loc_cells+=1
    after={str(p.relative_to(ROOT)):digest(p) for p in ROOT.rglob('*') if p.is_file()}
    assert before==after,'Package files changed during read-only verification'
    print(json.dumps(dict(status='PASS',python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
                          CON03_candidates=95,CON03_checked_cells=con_cells,
                          S5_rows=5,S5_trajectory_metric_rows=len(selected),S5_checked_cells=loc_cells,
                          max_absolute_numeric_difference=largest,tolerance=1e-12,
                          package_unchanged=True,
                          exclusions=['raw trajectory extraction','initial candidate ranking/selection','all-atom sensitivity calculation','structural QC generation']),indent=2))


if __name__=='__main__':main()
