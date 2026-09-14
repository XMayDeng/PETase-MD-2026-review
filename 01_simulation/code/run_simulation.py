#!/usr/bin/env python3
"""Inspect or replay an indexed production-MD input in the review package.

This starts a new execution from the stored production starting state. Earlier
preparation stages and the original trajectories are outside this package.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT.parent


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def catalog():
    with (ROOT / 'configurations/stages.csv').open(newline='') as f:
        records = list(csv.DictReader(f))
    seen = set()
    for r in records:
        identity = (r['unit'], r['role'], r['stage'])
        if identity in seen:
            raise ValueError('Duplicate stage identity')
        seen.add(identity)
        p = PACKAGE / r['tpr_path']
        if not p.resolve().is_relative_to(PACKAGE.resolve()) or not p.is_file():
            raise ValueError('Missing/unsafe package input')
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['list', 'check', 'replay'])
    p.add_argument('--name', help='Readable name, e.g. TfCut1_L10_bi/rep01')
    p.add_argument('--stage', choices=['MD'], default='MD')
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--gmx', default='gmx')
    p.add_argument('--steps', type=int, help='Bounded smoke test; not scientific production')
    p.add_argument('--full-stage', action='store_true', help='Execute the full recorded stage')
    p.add_argument('--gpu', type=int, help='Use this currently available GPU; otherwise CPU')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--workers', type=int, default=2, help='Concurrent read-only input checks')
    args = p.parse_args()
    records = [r for r in catalog() if (not args.name or args.name == r['name'])
               and (not args.stage or args.stage == r['stage'])]
    if not records:
        p.error('No matching stage')
    if args.action == 'list':
        for r in records:
            print(r['name'], r['stage'], r['role'], r['tpr_path'], sep='\t')
        return
    if args.output_dir is None:
        p.error('--output-dir is required')
    out = args.output_dir.resolve()
    if out.is_relative_to(PACKAGE.resolve()) or out.exists():
        p.error('Use a new output directory outside the package')
    if args.action == 'replay':
        if len(records) != 1:
            p.error('Replay requires exactly one --name and --stage')
        if (args.steps is None) == (not args.full_stage):
            p.error('Specify exactly one of --steps or --full-stage')
        if args.steps is not None and args.steps <= 0:
            p.error('--steps must be positive')
    if not 1 <= args.workers <= 4 or args.threads <= 0:
        p.error('Use 1--4 check workers and positive threads')
    out.mkdir(parents=True)
    env = dict(os.environ, GMX_MAXBACKUP='-1', OMP_NUM_THREADS=str(args.threads),
               PYTHONDONTWRITEBYTECODE='1')
    version = subprocess.run([args.gmx, '--version'], capture_output=True, text=True, check=True)
    (out / 'gromacs_version.txt').write_text(version.stdout + version.stderr)

    def check(r):
        inp = PACKAGE / r['tpr_path']
        if sha(inp) != r['sha256']:
            raise ValueError('Input checksum changed: ' + str(inp))
        work = out / r['name'] / r['role'] / r['stage']
        work.mkdir(parents=True)
        cmd = [args.gmx, 'dump', '-s', str(inp)]
        with (work / 'read_input.log').open('w') as f:
            proc = subprocess.run(cmd, cwd=work, env=env, stdout=subprocess.DEVNULL,
                                  stderr=f, timeout=180)
        if proc.returncode:
            raise RuntimeError('TPR read failed: ' + r['name'] + '/' + r['role'])
        return dict(name=r['name'], role=r['role'], stage=r['stage'], status='PASS', sha256=r['sha256'])

    if args.action == 'check':
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = []
            for result in pool.map(check, records):
                results.append(result)
                if len(results) % 50 == 0:
                    print(f'GROMACS read {len(results)}/{len(records)} TPR files', flush=True)
        (out / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        print(json.dumps(dict(status='PASS', checked_tprs=len(results), scope='GROMACS readability, not execution')))
        return
    r = records[0]
    check(r)
    inp = PACKAGE / r['tpr_path']
    cmd = [args.gmx, 'mdrun', '-s', str(inp), '-deffnm', 'test' if args.steps else r['stage'],
           '-ntomp', str(args.threads), '-pme', 'cpu', '-bonded', 'cpu', '-update', 'cpu']
    if args.gpu is None:
        cmd += ['-nb', 'cpu']
    else:
        cmd += ['-nb', 'gpu', '-gpu_id', str(args.gpu)]
    if args.steps:
        cmd += ['-nsteps', str(args.steps)]
    with (out / 'execution.log').open('w') as f:
        proc = subprocess.run(cmd, cwd=out, env=env, stdout=f, stderr=subprocess.STDOUT)
    log = out / ('test.log' if args.steps else r['stage'] + '.log')
    text = log.read_text(errors='replace') if log.exists() else ''
    result = dict(status='PASS' if proc.returncode == 0 and 'Finished mdrun' in text else 'FAIL',
                  name=r['name'], role=r['role'], source_tpr_sha256=r['sha256'],
                  scope='bounded_execution_smoke_test' if args.steps else 'full_recorded_stage_replay',
                  command=cmd, requested_steps=args.steps, gpu_requested=args.gpu,
                  input_preserved=sha(inp) == r['sha256'])
    (out / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    if result['status'] != 'PASS' or not result['input_preserved']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
