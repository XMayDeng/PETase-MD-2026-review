"""Recompute STR-SENS-01 aggregate statistics from bundled trajectory summaries.

The source RMSDs and frame-state classifications are inputs, not re-extracted.
The original whole-trajectory bootstrap seed/order and all valid zeros persist.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import str_sens

ROOT = Path(__file__).resolve().parents[3]
AGG = ROOT/'02_postprocessing/results/aggregate_statistics'
TRAJ = ROOT/'02_postprocessing/results/trajectory_summaries'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(row):
    return (row.protein_id, row.pet_kind, row.direction, int(row.replica_id))


def validate_inputs(trajectory, interaction):
    retention = pd.read_csv(TRAJ/'retention_manifest.csv')
    if retention.duplicated(['protein_id','pet_kind','direction','replica_id']).any():
        raise ValueError('Duplicate retention identity')
    master = {key(row):row for row in retention.itertuples()}
    expected_targets = set(str_sens.TARGET_ORDER)
    if len(trajectory) != 252 or trajectory.duplicated(['trajectory_id','target_resid','window']).any():
        raise ValueError('Expected 252 distinct target/trajectory/window records')
    if set(zip(trajectory.protein_id,trajectory.target_resid)) != expected_targets:
        raise ValueError('Incorrect target identities')
    if not pd.api.types.is_bool_dtype(trajectory.retained) or not pd.api.types.is_bool_dtype(trajectory.strict):
        raise ValueError('Retention/strict flags must be explicit booleans')
    for (protein,target),group in trajectory.groupby(['protein_id','target_resid'],sort=False):
        expected = {k for k in master if k[0] == protein}
        for window,frames in [('start',1),('early',101),('primary',801)]:
            subset = group[group.window == window]
            if len(subset) != 21 or {key(row) for row in subset.itertuples()} != expected:
                raise ValueError('Missing trajectory/window for '+repr((protein,target,window)))
            if not (subset.n_sampled_frames == frames).all():
                raise ValueError('Window sampling count changed')
    for row in trajectory.itertuples():
        reference = master[key(row)]
        if row.retained != (reference.retained_fingerprint_eligible == 'yes') or row.strict != (reference.strict_fingerprint_eligible == 'yes'):
            raise ValueError('Retention flag differs from the canonical manifest')
        if row.stability_verdict != reference.stability_verdict:
            raise ValueError('Stability label differs from the canonical manifest')
    required = ['mean_primary_production_rmsd_angstrom','mean_primary_experimental_rmsd_angstrom',
                'mean_primary_reference_closeness_delta_angstrom','fraction_primary_experimental_overlap',
                'fraction_primary_experimental_closer']
    if not np.isfinite(trajectory[required].to_numpy(float)).all():
        raise ValueError('Missing required trajectory measurements')
    if not trajectory[['fraction_primary_experimental_overlap','fraction_primary_experimental_closer']].apply(lambda x:x.between(0,1).all()).all():
        raise ValueError('Fraction outside [0,1]')
    primary = trajectory[trajectory.window == 'primary']
    allowed = {(row.trajectory_id,int(row.target_resid),'pet_heavy_atom_contact')
               for row in primary.itertuples() if row.retained}
    allowed |= {(row.trajectory_id,166,'t166_native_sidechain_hbond')
                for row in primary.itertuples() if row.retained and row.protein_id == 'Pro00121' and row.target_resid == 166}
    actual = set(zip(interaction.trajectory_id,interaction.target_resid,interaction.interaction))
    if len(interaction) != 65 or len(actual) != len(interaction) or actual != allowed:
        raise ValueError('Incomplete/duplicate interaction-trajectory coverage')
    if not (interaction.n_sampled_frames == 801).all():
        raise ValueError('Interaction window must contain 801 sampled frames')
    for row in interaction.itertuples():
        record = master[key(row)]
        if record.retained_fingerprint_eligible != 'yes':
            raise ValueError('Interaction summary includes a nonretained trajectory')
        n,e = int(row.n_sampled_frames),int(row.n_event_frames)
        if row.n_event_frames != e or not np.isfinite(row.event_fraction) or not 0 <= e <= n or abs(row.event_fraction-e/n) > 1e-12:
            raise ValueError('Event count/fraction mismatch')
        state = count_from_fraction(row.experimental_closer_fraction,n)
        exp,prod,delta = (row.event_fraction_in_experimental_closer_state,
                          row.event_fraction_in_production_closer_state,row.event_fraction_exp_minus_prod_state)
        if not 0 <= state <= n:
            raise ValueError('Invalid reference-state count')
        if state == 0:
            if not np.isnan(exp) or not np.isnan(delta) or not np.isfinite(prod) or abs(prod-e/n)>1e-12:
                raise ValueError('Absent experimental state must remain undefined')
        elif state == n:
            if not np.isnan(prod) or not np.isnan(delta) or not np.isfinite(exp) or abs(exp-e/n)>1e-12:
                raise ValueError('Absent production state must remain undefined')
        elif not np.isfinite([exp,prod,delta]).all() or abs(delta-(exp-prod))>1e-12:
            raise ValueError('Missing/inconsistent within-trajectory state contrast')
        if np.isfinite(exp) and not 0 <= exp <= 1 or np.isfinite(prod) and not 0 <= prod <= 1:
            raise ValueError('State-conditioned event fraction outside [0,1]')
        if e == 0:
            if not np.isnan(row.experimental_closer_fraction_among_event_frames):
                raise ValueError('No-event conditional fraction must remain undefined')
        else:
            count_from_fraction(row.experimental_closer_fraction_among_event_frames,e)


def count_from_fraction(fraction, count):
    if not np.isfinite(fraction) or not 0 <= fraction <= 1:
        raise ValueError('Undefined or invalid required count fraction')
    value = float(fraction)*int(count)
    rounded = round(value)
    if abs(value-rounded)>1e-8:
        raise ValueError('Fraction does not encode an integer frame count')
    return rounded


def coupling_summary(trajectory):
    """Same grouped bootstrap; reconstruct pooled counts from exact count fractions."""
    records = []
    seed_index = 500
    for (protein,target,interaction),group in trajectory.groupby(['protein_id','target_resid','interaction'],sort=False):
        values = group.event_fraction_exp_minus_prod_state.dropna().to_numpy()
        if len(values):
            mean,low,high = str_sens.bootstrap_mean_ci(values,str_sens.BOOTSTRAP_SEED+seed_index)
            seed_index += 1
        else:
            mean=low=high=float('nan')
        n = int(group.n_sampled_frames.sum())
        event_n = int(group.n_event_frames.sum())
        state_n = sum(count_from_fraction(row.experimental_closer_fraction,row.n_sampled_frames) for row in group.itertuples())
        state_event_n = sum(count_from_fraction(row.experimental_closer_fraction_among_event_frames,row.n_event_frames)
                            for row in group.itertuples() if row.n_event_frames>0)
        first = group.iloc[0]
        records.append(dict(protein_id=protein,protein=first.protein,target_resid=int(target),target_label=first.target_label,
            interaction=interaction,n_trajectories=len(group),n_trajectories_with_event=int((group.n_event_frames>0).sum()),
            n_trajectories_with_both_reference_states=int(group.event_fraction_exp_minus_prod_state.notna().sum()),
            total_sampled_frames=n,total_event_frames=event_n,pooled_event_fraction=event_n/n,
            pooled_experimental_closer_fraction=state_n/n,
            pooled_experimental_closer_fraction_among_event_frames=state_event_n/event_n if event_n else float('nan'),
            mean_trajectory_event_fraction_exp_minus_prod_state=mean,bootstrap_low=low,bootstrap_high=high))
    return pd.DataFrame(records)


def compare(actual, expected, keys):
    if list(actual.columns) != list(expected.columns):
        raise ValueError('Unexpected output schema')
    if actual.duplicated(keys).any() or expected.duplicated(keys).any():
        raise ValueError('Duplicate output identities')
    a = actual.set_index(keys).sort_index()
    e = expected.set_index(keys).sort_index()
    pd.testing.assert_frame_equal(a,e,check_dtype=False,check_exact=False,rtol=0,atol=1e-12)
    largest = 0.0
    numeric_cells = 0
    for column in a.columns:
        if pd.api.types.is_numeric_dtype(a[column]):
            delta = np.abs(a[column].to_numpy(float)-e[column].to_numpy(float))
            valid = np.isfinite(delta)
            if valid.any():
                largest = max(largest,float(delta[valid].max()))
            numeric_cells += len(a)
    return numeric_cells,largest


def calculate():
    record = json.loads((ROOT/'metadata/analyses.json').read_text())['STR-SENS-01']
    for path,sha in record['inputs'].items():
        if digest(ROOT/path) != sha:
            raise ValueError('Registered analysis input changed: '+path)
    for path,sha in record['references'].items():
        if digest(ROOT/path) != sha:
            raise ValueError('Registered numerical reference changed: '+path)
    trajectory = pd.read_csv(TRAJ/'trajectory_target_summary.csv')
    interaction = pd.read_csv(TRAJ/'trajectory_interaction_state_coupling.csv')
    validate_inputs(trajectory,interaction)
    spec = json.loads((ROOT/'02_postprocessing/inputs/str_sens_statistics.json').read_text())
    outputs = {'target_subset_summary.csv':str_sens.target_summaries(trajectory,spec),
               'interaction_state_coupling_summary.csv':coupling_summary(interaction)}
    counts=0
    largest=0.0
    keys = {'target_subset_summary.csv':['protein_id','target_resid','subset'],
            'interaction_state_coupling_summary.csv':['protein_id','target_resid','interaction']}
    for name,data in outputs.items():
        cells,difference = compare(data,pd.read_csv(AGG/name),keys[name])
        counts += cells
        largest = max(largest,difference)
    base=pd.read_csv(AGG/'static_reference_reconciliation.csv').set_index(['protein_id','target_resid'])
    summary=outputs['target_subset_summary.csv'].set_index(['protein_id','target_resid','subset'])
    coupling=outputs['interaction_state_coupling_summary.csv'].set_index(['protein_id','target_resid','interaction'])
    projected=[]
    for protein,target in str_sens.TARGET_ORDER:
        b=base.loc[(protein,target)]
        s=summary.loc[(protein,target,'all')]
        interaction_name='t166_native_sidechain_hbond' if (protein,target)==('Pro00121',166) else 'pet_heavy_atom_contact'
        c=coupling.loc[(protein,target,interaction_name)]
        if b['pass'] != True or not np.isfinite([b.recomputed_production_to_experimental_rmsd_angstrom,
                                                b.absolute_reconciliation_difference_angstrom]).all() or float(b.absolute_reconciliation_difference_angstrom)>1e-12:
            raise ValueError('Static-reference reconciliation did not pass')
        projected.append(dict(target=s.target_label,
            metric='local backbone after pocket fit' if s.primary_metric=='target_backbone_after_pocket_fit' else 'side-chain conformation',
            initial_model_reference_rmsd_A=b.recomputed_production_to_experimental_rmsd_angstrom,
            experimental_overlap_trajectories=s.n_trajectories_exp_overlap_ge_0p05,total_trajectories=s.n_trajectories,
            median_overlap_fraction=s.median_primary_experimental_overlap_fraction,
            median_experimental_closer_fraction=s.median_primary_experimental_closer_fraction,
            mean_analysis_minus_start_experimental_rmsd_A=s.mean_analysis_minus_start_experimental_rmsd_angstrom,
            interaction_indicator='not estimable' if (protein,target)==('Pro00057',185) else
                'native side-chain hydrogen bond' if interaction_name=='t166_native_sidechain_hbond' else 'heavy-atom contact',
            interaction_state_contrast=c.mean_trajectory_event_fraction_exp_minus_prod_state,
            interaction_state_ci95_low=c.bootstrap_low,interaction_state_ci95_high=c.bootstrap_high))
    table=pd.read_csv(ROOT/'02_postprocessing/results/tables/Supplementary_Table_S4_local_starting_coordinate_sensitivity.csv')
    cells,difference=compare(pd.DataFrame(projected),table,['target'])
    counts+=cells
    largest=max(largest,difference)
    return outputs,dict(status='PASS',analysis='STR-SENS-01',unique_trajectories=trajectory.trajectory_id.nunique(),
        target_trajectory_windows=len(trajectory),interaction_trajectory_rows=len(interaction),
        target_subset_rows=12,interaction_summary_rows=5,table_s4_rows=4,numeric_cells=counts,
        max_absolute_numeric_difference=largest,tolerance=1e-12,bootstrap_samples=str_sens.BOOTSTRAP_SAMPLES,
        undefined_H185_contrast_preserved=True,scope='trajectory_summaries_to_aggregate_statistics',
        not_rerun=['coordinate extraction','RMSD calculation','per-frame state/contact/hydrogen-bond assignment'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='Read-only recomputation (also the default)')
    parser.add_argument('--output-dir',type=Path,help='New directory for actual recomputed CSVs; numeric equivalence, not byte identity')
    args=parser.parse_args()
    if args.check and args.output_dir is not None:
        parser.error('--check and --output-dir are mutually exclusive')
    outputs,report=calculate()
    if args.output_dir is not None:
        directory=args.output_dir.resolve()
        if directory.exists():
            raise FileExistsError('Use a new output directory; no overwrite: '+str(directory))
        encoded={name:data.to_csv(index=False,lineterminator='\n').encode() for name,data in outputs.items()}
        directory.mkdir(parents=True,exist_ok=False)
        for name,data in encoded.items():
            with (directory/name).open('xb') as handle:
                handle.write(data)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
