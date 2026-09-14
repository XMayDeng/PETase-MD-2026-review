#!/usr/bin/env python3
"""Recompute declared statistics from delivered observations, never reference outputs."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent/'statistics'))
import build_pooled_top8_unfiltered as pool8
import build_pooled_top20_composition_unfiltered as pool20
import class_uncertainty
import pet_rg
import wloop
import phl7_rotamers
import candidate_hbonds
import hbond_angles
import contact_composition
import loc01


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--groups',nargs='+',choices=['contacts','class','rg','wloop','rotamers','hbonds','local','all'],default=['all'])
    args=parser.parse_args()
    out=args.output_dir.resolve()
    if out.exists() or ROOT==out or ROOT in out.parents: parser.error('Use a new output directory outside the package')
    out.mkdir(parents=True)
    base=ROOT/'02_postprocessing'
    derived=base/'inputs/derived_data'
    results=base/'results'
    checks=[]
    groups={'contacts','class','rg','wloop','rotamers','hbonds','local'} if 'all' in args.groups else set(args.groups)
    def read(folder,name): return pd.read_csv(folder/(name+'.csv'),float_precision='round_trip')
    def compare(frame,name,folder='aggregate_statistics',keys=None,atol=1e-12):
        frame.to_csv(out/(name+'.csv'),index=False)
        expected=read(results/folder,name)
        if set(frame.columns)!=set(expected.columns):
            raise ValueError(f'{name}: column mismatch {set(frame.columns)^set(expected.columns)}')
        actual=frame[expected.columns]
        for col in actual:
            if isinstance(actual[col].dtype,pd.CategoricalDtype):actual=actual.assign(**{col:actual[col].astype(object)})
        if keys:
            assert not actual.duplicated(keys).any() and not expected.duplicated(keys).any()
            actual=actual.sort_values(keys).reset_index(drop=True)
            expected=expected.sort_values(keys).reset_index(drop=True)
        assert actual.shape==expected.shape,(name,actual.shape,expected.shape)
        max_error=0.;failures=[]
        for col in expected:
            a=actual[col];b=expected[col]
            if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
                av=a.to_numpy(float);bv=b.to_numpy(float)
                valid=np.isfinite(av)&np.isfinite(bv)
                error=float(np.max(np.abs(av[valid]-bv[valid]))) if valid.any() else 0.
                max_error=max(max_error,error)
                if not np.allclose(av,bv,atol=atol,rtol=0,equal_nan=True): failures.append(dict(column=col,max_error=error))
            elif not a.astype(object).fillna('').astype(str).equals(b.astype(object).fillna('').astype(str)):
                failures.append(dict(column=col,kind='categorical_or_text'))
        checks.append(dict(output=name,status='FAIL' if failures else 'PASS',rows=len(actual),
                           numeric_tolerance=atol,maximum_absolute_numeric_error=max_error,failures=failures))
        print(json.dumps(checks[-1]),flush=True)
    if groups & {'contacts','class'}:
        contacts=read(derived,'residue_contact_fractions')
        manifest=read(derived,'retained_case_metadata')
        annotation=read(base/'inputs/annotations','full_panel_contact_landmark_annotation')
    if 'contacts' in groups:
        joined,retained,population,qa=pool8.validate_and_join(contacts,manifest,annotation)
        assert all(r['status']=='PASS' for r in qa),qa
        pooled,top8=pool8.build_pooled_summary(joined)
        top20,composition,qa=pool20.build_composition(pooled,population)
        assert all(r['status']=='PASS' for r in qa),qa
        compare(top8,'pooled_top8_contacts_unfiltered',keys=['protein_id','residue'])
        compare(composition,'pooled_top20_contact_class_composition_unfiltered',keys=['protein_id','landmark_class'])
        shares=contact_composition.build_share_table(contacts[['case','residue','contact_fraction','contact_count','total_frames']],
                retained,annotation,min_contact_fraction=.20,analysis_label='contact_fraction_ge_0_2')
        compare(shares,'pet_length_landmark_share_filtered_contact_ge_0_2',keys=['protein_id','pet_kind','landmark_class'])
    if 'class' in groups:
        retained=manifest[pool8.bool_series(manifest.retained_bool)].copy()
        mass=class_uncertainty.build_case_class_mass(contacts[['case','residue','contact_fraction','contact_count','total_frames']],
                    retained[['case','protein_id','protein_label','pet_kind','direction','replica_id','stability_pass_bool']],
                    annotation[['protein_id','residue','landmark_class']].drop_duplicates())
        # This is an independently derived point check, not a bootstrap reference.
        shares=contact_composition.build_share_table(
            contacts[['case','residue','contact_fraction','contact_count','total_frames']],
            retained,annotation,min_contact_fraction=.20,analysis_label='contact_fraction_ge_0_2')
        shares['protein_short']=shares.protein_short.replace({'CUT1':'TfCut1'})
        wide=shares.pivot(index=['protein_short','landmark_class'],columns='pet_kind',values='contact_mass_share')
        wide=wide.reindex(pd.MultiIndex.from_product([class_uncertainty.PROTEIN_NAMES.values(),class_uncertainty.CLASS_ORDER],
                            names=['protein_short','landmark_class'])).fillna(0).reset_index()
        points=wide[['protein_short','landmark_class']].copy()
        points['comparison']='PET_L20_minus_PET_L10'
        points['delta_share_pct_points']=100*(wide.PET_L20-wide.PET_L10)
        compare(class_uncertainty.analyze(mass,points),'class_contact_shift_bootstrap',keys=['protein_id','landmark_class'])
    if 'rg' in groups:
        group,comparison=pet_rg.summarize_groups(read(results/'trajectory_summaries','pet_rg_trajectory_summary'))
        group.to_csv(out/'pet_rg_group_summary.csv',index=False)
        compare(comparison,'pet_rg_group_comparison')
    if 'wloop' in groups:
        compare(wloop.summarize_profile(read(derived,'wloop_residue_rmsf')),'protein_wloop_profile_summary',
                keys=['protein_id','relative_position'])
        trajectory=read(results/'trajectory_summaries','trajectory_wloop_summary')
        summary=pd.concat([wloop.summarize_metric_by_protein(trajectory,m,m)
                    for m in ['wloop_mean_rmsf_nm','marker_rmsf_nm']],ignore_index=True)
        compare(summary,'protein_wloop_metric_summary',keys=['protein_id','metric'])
        compare(wloop.summarize_pair_sites(read(derived,'pair_site_rmsf')),'mapped_pair_site_summary',
                keys=['protein_id','residue'])
    if 'rotamers' in groups:
        frames=pd.read_csv(results/'frame_observables/frame_level_w156_rotamer_contacts.csv.gz')
        trajectory=phl7_rotamers.trajectory_summary(frames)
        # Frame angles are distributed to six decimal places; output references
        # were originally calculated before that export rounding. Report errors.
        compare(trajectory,'trajectory_rotamer_summary','trajectory_summaries',keys=['case'],atol=1e-6)
        group,difference=phl7_rotamers.bootstrap_group_and_difference(trajectory,20000,20260730)
        compare(group,'length_group_summary',keys=['metric','pet_kind'],atol=1e-6)
        difference.to_csv(out/'length_difference_bootstrap.csv',index=False)
        compare(phl7_rotamers.mixture_model_sensitivity(frames,20260730),'mixture_model_sensitivity',atol=1e-6)
    if 'hbonds' in groups:
        values=read(results/'trajectory_summaries','candidate_case_hbond_values_zero_filled')
        # Original occupancy = integer frame count / total_frames. Eight-place
        # CSV rounding uniquely identifies that integer (8001 frames); recover
        # it before strict-sign bootstrap probabilities, without changing data.
        for col in ['any_hbond_occupancy','sidechain_hbond_occupancy','backbone_hbond_occupancy',
                    'protein_donor_hbond_occupancy','protein_acceptor_hbond_occupancy']:
            count=np.rint(values[col]*values.total_frames).astype(int)
            exact=count/values.total_frames
            assert np.allclose(exact,values[col],atol=5.01e-9,rtol=0)
            if col=='any_hbond_occupancy': assert np.array_equal(count,values.frames_with_any_candidate_hbond)
            values[col]=exact
        values['mean_triplet_count_per_frame']=values.hbond_triplet_events/values.total_frames
        compare(candidate_hbonds.summarize_candidate_pet(values),'candidate_pet_hbond_summary',
                keys=['protein_id','residue_number','pet_kind'],atol=1e-8)
        shifts,directions=candidate_hbonds.analyze_shifts(values)
        compare(shifts,'candidate_hbond_shift_uncertainty',keys=['protein_id','residue_number'],atol=1e-8)
        directions.to_csv(out/'candidate_hbond_direction_sensitivity.csv',index=False)
        angles=read(derived,'criterion_case_hbond_occupancy')
        compare(hbond_angles.build_shift_statistics(angles),'criterion_residue_shift_statistics',
                keys=['criterion','protein_id','residue_number','metric'])
    if 'local' in groups:
        metrics=read(results/'trajectory_summaries','trajectory_metric_values')
        compare(loc01.build_shift_statistics(metrics),'interaction_shift_statistics',keys=['analysis_target','metric_id'])
    report=dict(status='PASS' if all(c['status']=='PASS' for c in checks) else 'FAIL',checks=checks,
                scope='Delivered derived observations to numerical summaries; no raw trajectory extraction')
    (out/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
    if report['status']!='PASS': raise SystemExit(1)


if __name__=='__main__':main()
