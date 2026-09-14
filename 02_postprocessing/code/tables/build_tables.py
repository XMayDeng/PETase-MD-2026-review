"""Assemble five registered tables without changing their approved contents.

Table 1/S1 are curated-cell assembly. S2/S3 project registered measurements.
S4 recomputes its aggregates and bootstraps from trajectory summaries first.
This command does not extract trajectory observations or calculate RMSDs.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import sys

import contact_shift_table
import structural_concordance

ROOT = Path(__file__).resolve().parents[3]
AGG = ROOT / '02_postprocessing/results/aggregate_statistics'
sys.path.insert(0,str(ROOT/'02_postprocessing/code/statistics'))
import recompute_structural_sensitivity


def digest(data):
    return hashlib.sha256(data).hexdigest()


def rows(path):
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        result = list(reader)
        if not result or not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Empty table or duplicate header: '+str(path))
        if any(None in row or any(value is None for value in row.values()) for row in result):
            raise ValueError('Malformed CSV row: '+str(path))
        return result


def unique(records, keys):
    result = {}
    for row in records:
        key = tuple(row[name] for name in keys)
        if key in result:
            raise ValueError('Duplicate identity: '+repr(key))
        result[key] = row
    return result


def finite(row, name):
    number = float(row[name])
    if not math.isfinite(number):
        raise ValueError('Missing/nonfinite required value: '+name)
    return number


def s4_rows():
    targets = [('Pro00057','185'), ('Pro00075','184'), ('Pro00121','164'), ('Pro00121','166')]
    base = unique(rows(AGG/'static_reference_reconciliation.csv'), ['protein_id','target_resid'])
    recomputed, _ = recompute_structural_sensitivity.calculate()
    summary = unique(list(csv.DictReader(io.StringIO(recomputed['target_subset_summary.csv'].to_csv(index=False)))),
                     ['protein_id','target_resid','subset'])
    coupling = unique(list(csv.DictReader(io.StringIO(recomputed['interaction_state_coupling_summary.csv'].to_csv(index=False)))),
                      ['protein_id','target_resid','interaction'])
    if set(base) != set(targets) or set(summary) != {(*key,subset) for key in targets for subset in ['all','retained','retained_l10_l20']}:
        raise ValueError('S4 target/subset coverage mismatch')
    expected_interactions = {(*key,'pet_heavy_atom_contact') for key in targets} | {('Pro00121','166','t166_native_sidechain_hbond')}
    if set(coupling) != expected_interactions:
        raise ValueError('S4 interaction coverage mismatch')
    result = []
    for key in targets:
        b, s = base[key], summary[(*key,'all')]
        interaction = 't166_native_sidechain_hbond' if key == ('Pro00121','166') else 'pet_heavy_atom_contact'
        c = coupling[(*key,interaction)]
        if b['pass'] != 'True' or s['n_trajectories'] != '21' or b['primary_metric'] != s['primary_metric']:
            raise ValueError('S4 reference/population mismatch: '+repr(key))
        contrast_names = ['mean_trajectory_event_fraction_exp_minus_prod_state','bootstrap_low','bootstrap_high']
        if key == ('Pro00057','185'):
            if c['n_trajectories_with_both_reference_states'] != '0' or any(c[name] != '' for name in contrast_names):
                raise ValueError('H185 must retain its undefined within-trajectory contrast')
            event = 'Not estimable ($n=0$)'
        else:
            values = [finite(c,name) for name in contrast_names]
            label = 'Side-chain H bond' if interaction == 't166_native_sidechain_hbond' else 'Heavy-atom contact'
            eligible_n = int(c['n_trajectories_with_both_reference_states'])
            if eligible_n <= 0:
                raise ValueError('An estimable event contrast requires eligible trajectories')
            label += f' ($n={eligible_n}$)'
            event = (r'\shortstack[l]{'+label+r':\\'+f'${values[0]:+.3f}$ (${values[1]:+.3f}$, ${values[2]:+.3f}$)'+'}')
        result.append([
            s['target_label'], 'Local backbone' if s['primary_metric'] == 'target_backbone_after_pocket_fit' else 'Side-chain conformation',
            f"{finite(b,'recomputed_production_to_experimental_rmsd_angstrom'):.3f}",
            s['n_trajectories_exp_overlap_ge_0p05']+'/21',
            f"{finite(s,'median_primary_experimental_overlap_fraction'):.3f}",
            f"{finite(s,'median_primary_experimental_closer_fraction'):.3f}",
            f"{finite(s,'mean_analysis_minus_start_experimental_rmsd_angstrom'):+.3f}", event])
    return result


def s3_rows():
    raw = rows(structural_concordance.SUMMARY_CSV)
    if len(raw) != 6 or len(unique(raw,['protein_id'])) != 6:
        raise ValueError('Expected exactly six unique structural summaries')
    output = []
    for row in structural_concordance.build_rows():
        output.append([
            row['enzyme'], row['experimental_reference']+f" ({row['experimental_resolution_angstrom']} \\AA)",
            f"{row['aligned_residues_n']} ({row['aligned_sequence_identity_percent']})",
            row['production_model_coverage_percent'], row['global_calpha_rmsd_angstrom'],
            row['pocket_backbone_rmsd_angstrom'],
            row['pocket_sidechain_conformation_rmsd_median_angstrom']+'/'+row['pocket_sidechain_conformation_rmsd_p90_angstrom'],
            row['maximum_triad_distance_delta_angstrom'], structural_concordance.latex_flag(row['prespecified_local_flag'])])
    return output


def build(key, record):
    for name, expected in record['inputs'].items():
        if digest((ROOT/name).read_bytes()) != expected:
            raise ValueError('Registered input changed: '+name)
    if key == 'Table_S02':
        data = rows(AGG/'manuscript_candidate_residues.csv')
        if len(data) != 29 or len(unique(data,['protein_id','residue'])) != 29:
            raise ValueError('Expected 29 unique approved S2 candidates')
        result = contact_shift_table.build_table(data).encode()
    else:
        if key in ['Table_01','Table_S01']:
            path = ROOT/f'02_postprocessing/inputs/table_assembly/{key}_approved_cells.csv'
            data = rows(path)
            if len(data) != 6:
                raise ValueError('Expected six curated rows')
            cells = [list(row.values()) for row in data]
        elif key == 'Table_S03':
            cells = s3_rows()
        elif key == 'Table_S04':
            cells = s4_rows()
        else:
            raise ValueError('Unregistered table: '+key)
        body = '\n'.join('    '+' & '.join(row)+r' \\' for row in cells)
        template = (ROOT/record['template']).read_text()
        if template.count('{{ROWS}}') != 1:
            raise ValueError('Expected exactly one table-body placeholder')
        result = template.replace('{{ROWS}}',body).encode()
    if digest(result) != record['reference_sha256']:
        raise ValueError('Generated table differs from the approved endpoint: '+key)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables',nargs='+',default=['all'])
    parser.add_argument('--check',action='store_true',help='Compare existing outputs without writing')
    parser.add_argument('--output-dir',type=Path,help='Destination for generated table environments')
    args = parser.parse_args()
    registry = json.loads((ROOT/'metadata/tables.json').read_text())
    keys = sorted(registry) if args.tables == ['all'] else args.tables
    if len(set(keys)) != len(keys) or any(key not in registry for key in keys):
        parser.error('Unknown or duplicate table identifiers')
    generated = []
    for key in keys:
        record = registry[key]
        target = ROOT/record['output'] if args.output_dir is None else args.output_dir.resolve()/Path(record['output']).name
        data = build(key,record)
        if target.exists():
            if target.read_bytes() != data:
                raise ValueError('Existing output differs; no overwrite: '+str(target))
        elif args.check:
            raise FileNotFoundError('Missing output in check mode: '+str(target))
        generated.append((key,target,data))
    for key,target,data in generated:
        if not args.check and not target.exists():
            target.parent.mkdir(parents=True,exist_ok=True)
            with target.open('xb') as handle:
                handle.write(data)
    print(json.dumps(dict(status='PASS',tables=keys,byte_exact=True,
                         scopes={key:registry[key]['scope'] for key in keys}),indent=2))


if __name__ == '__main__':
    main()
