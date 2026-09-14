"""Check packaged preparation metadata, not MD execution or chemical validity."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / '01_simulation/metadata'


def rows(path):
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        result = list(reader)
    if not result or any(None in row or None in row.values() for row in result):
        raise ValueError('Empty or malformed metadata: ' + str(path))
    return result


def indexed(records, key):
    result = {row[key]: row for row in records}
    if len(result) != len(records):
        raise ValueError('Duplicate identity: ' + key)
    return result


def validate_mdps(trajectories):
    records = rows(ROOT / '01_simulation/configurations/mdp_index.csv')
    mdps = {(r['trajectory_id'], r['stage']): r for r in records}
    expected = {(key, stage) for key in trajectories for stage in ('NVT', 'NPT', 'MD')}
    if len(mdps) != len(records) or set(mdps) != expected:
        raise ValueError('Expected exactly one NVT/NPT/MD record per trajectory')
    direction_names = {'single': 'single', 'head-side': 'head',
                       'tail-side': 'tail', 'bidirectional': 'bi'}
    paths = set()
    for (key, stage), record in mdps.items():
        t = trajectories[key]
        relative = (f"01_simulation/configurations/mdp/{t['enzyme']}/L{t['pet_length']}/"
                    f"{direction_names[t['direction']]}/rep{int(t['replica_id']):02d}/{stage}.mdp")
        if record['mdp_path'] != relative:
            raise ValueError('MDP path does not match trajectory identity: ' + key)
        path = ROOT / relative
        if (path.is_symlink() or not path.resolve().is_relative_to(ROOT.resolve())
                or not path.is_file()):
            raise ValueError('MDP must be a real file inside the package: ' + relative)
        if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
            raise ValueError('MDP content hash mismatch: ' + relative)
        paths.add(path)
    actual = set((ROOT / '01_simulation/configurations/mdp').rglob('*.mdp'))
    if len(paths) != len(records) or paths != actual:
        raise ValueError('MDP files are missing, duplicated or not indexed')
    return mdps


def main():
    trajectories = indexed(rows(ROOT / 'metadata/trajectory_index.csv'), 'trajectory_id')
    seeds = indexed(rows(DATA / 'simulation_seeds.csv'), 'trajectory_id')
    pme = indexed(rows(DATA / 'production_pme_tuning.csv'), 'trajectory_id')
    systems = indexed(rows(DATA / 'simulation_systems.csv'), 'condition_id')
    models = indexed(rows(DATA / 'protein_sources.csv'), 'protein_id')
    preparation = indexed(rows(DATA / 'protein_preparation.csv'), 'enzyme')
    if len(trajectories) != 126 or set(seeds) != set(trajectories) or set(pme) != set(trajectories):
        raise ValueError('Seed/PME coverage must equal the 126 registered trajectories')
    conditions = {key.rsplit('__r', 1)[0] for key in trajectories}
    if len(conditions) != 42 or set(systems) != conditions:
        raise ValueError('System coverage must equal the 42 registered conditions')
    enzymes = {row['enzyme'] for row in trajectories.values()}
    if len(models) != 6 or {row['enzyme'] for row in models.values()} != enzymes or set(preparation) != enzymes:
        raise ValueError('Model/preparation coverage must equal the six registered enzymes')
    for model in models.values():
        path = (ROOT / model['model_path']).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('Model must be a real package file')
        if hashlib.sha256(path.read_bytes()).hexdigest() != model['sha256']:
            raise ValueError('Source-model hash mismatch: ' + model['protein_id'])
        if model['sha256'] != preparation[model['enzyme']]['source_model_sha256']:
            raise ValueError('Preparation/source-model identity mismatch')
    mdps = validate_mdps(trajectories)
    for key, seed in seeds.items():
        text = (ROOT / mdps[key, 'NVT']['mdp_path']).read_text()
        parameters = {}
        for line in text.splitlines():
            line = line.split(';', 1)[0]
            if '=' in line:
                name, value = line.split('=', 1)
                parameters[name.strip().replace('_', '-').lower()] = value.strip()
        if int(parameters['gen-seed']) != int(seed['nvt_gen_seed']):
            raise ValueError('NVT seed/MDP mismatch: ' + key)
        if int(trajectories[key]['nvt_gen_seed']) != int(seed['nvt_gen_seed']):
            raise ValueError('NVT seed/trajectory index mismatch: ' + key)
        int(seed['production_ld_seed'])
        if len(seed['production_tpr_sha256']) != 64:
            raise ValueError('Missing production TPR provenance: ' + key)
    for condition, system in systems.items():
        members = [key for key in trajectories if key.rsplit('__r', 1)[0] == condition]
        if len(members) != 3 or len({seeds[key]['nvt_gen_seed'] for key in members}) != 3:
            raise ValueError('Expected three replica seeds per condition: ' + condition)
        if int(system['replicas']) != 3 or int(system['system_atoms']) <= int(system['protein_atoms']):
            raise ValueError('Invalid system composition: ' + condition)
    for path in DATA.glob('*.csv'):
        if any(token in path.read_text() for token in ('/workspace/', '/scratch/', '/home/')):
            raise ValueError('Internal absolute path in preparation metadata: ' + path.name)
    print(json.dumps(dict(status='PASS', models=6, systems=42, trajectories=126,
                         mdp_files=len(mdps), readable_mdp_identity_and_hashes_checked=True,
                         unique_nvt_seeds=len({r['nvt_gen_seed'] for r in seeds.values()}),
                         scope='source identity, metadata coverage and packaged NVT seed consistency',
                         tpr_reparsed=False, md_rerun=False), indent=2))


if __name__ == '__main__':
    main()
