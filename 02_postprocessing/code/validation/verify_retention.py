"""Verify RET-01 derived tables without reading trajectories or writing files."""
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
ORDER = ['IsPETase', 'TfCut1', 'LCC', 'FoCut5a', 'HiC', 'PHL7']
NAMES = dict(zip(['Pro00083', 'Pro00057', 'Pro00062', 'Pro00075', 'Pro00121', 'Pro00137'], ORDER))
DIRECTIONS = {'bidirectional': 'bi', 'head-side': 'head', 'tail-side': 'tail', 'single': 'single'}


def tables(frame):
    frame = frame.copy()
    frame['protein_short'] = frame.protein_id.map(NAMES)
    assert frame.protein_short.notna().all()
    frame['condition'] = ['L%d/%s' % (n, DIRECTIONS[d]) for n, d in zip(frame.pet_length, frame.direction)]
    condition = frame.groupby(['protein_id', 'protein_short', 'condition'], as_index=False).agg(
        mean_retention_fraction=('retention_fraction', 'mean'),
        retained_ready_n=('retained_fingerprint_eligible', lambda x: int(x.eq('yes').sum())),
        trajectory_n=('retention_fraction', 'size'))
    counts = []
    for name in ORDER:
        group = frame[frame.protein_short.eq(name)]
        counts.append(dict(protein_short=name,
            detached=int(group.retention_class.eq('detached').sum()),
            partial=int(group.retention_class.eq('partial').sum()),
            retained_warn=int((group.retention_class.eq('retained') & group.stability_verdict.eq('WARN')).sum()),
            retained_pass=int((group.retention_class.eq('retained') & group.stability_verdict.eq('PASS')).sum())))
    return condition, pd.DataFrame(counts)


def main():
    before = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in ROOT.rglob('*') if p.is_file()}
    frame = pd.read_csv(ROOT/'02_postprocessing/results/trajectory_summaries/retention_manifest.csv')
    master = pd.read_csv(ROOT/'metadata/trajectory_index.csv')
    keys = ['protein_id', 'pet_length', 'direction', 'replica_id']
    assert len(frame) == len(master) == 126
    assert not frame.duplicated(keys).any() and not master.duplicated(keys).any()
    merged = frame.merge(master, on=keys, suffixes=('', '_index'), validate='one_to_one')
    assert len(merged) == 126
    for col in ['run_id', 'retention_fraction', 'retention_class', 'stability_verdict', 'pbc_verdict', 'clash_verdict']:
        assert merged[col].equals(merged[col+'_index']), col
    assert frame.retention_fraction.between(0, 1).all() and frame.burn_in_ns.eq(20).all()
    expected = frame.retention_fraction.map(lambda x: 'retained' if x >= .8 else 'partial' if x >= .5 else 'detached')
    assert expected.equals(frame.retention_class)
    assert frame.retention_class.value_counts().to_dict() == {'retained': 71, 'detached': 39, 'partial': 16}
    assert frame.retained_fingerprint_eligible.eq('yes').equals(expected.eq('retained'))
    conditions, counts = tables(frame)
    assert len(conditions) == 42 and conditions.trajectory_n.eq(3).all()
    assert counts.iloc[:, 1:].sum(axis=1).eq(21).all()
    for table, name in [(conditions, 'retention_conditions.csv'), (counts, 'retention_qc_counts.csv')]:
        target = pd.read_csv(ROOT/'02_postprocessing/results/aggregate_statistics'/name)
        pd.testing.assert_frame_equal(table, target, check_dtype=False, atol=1e-15, rtol=0)
    after = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in ROOT.rglob('*') if p.is_file()}
    assert before == after
    print(json.dumps(dict(status='PASS', trajectories=126, conditions=42, retained=71,
                         partial=16, detached=39, package_unchanged=True,
                         scope='Derived-data aggregation only; no distance-series or trajectory recomputation'), indent=2))


if __name__ == '__main__':
    main()
