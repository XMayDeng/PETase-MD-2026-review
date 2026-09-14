# Protein sources and executed-system metadata

This local preparation supplement contains six source-model identities,
42 system-composition records, and 126 replica seed/PME records. The source
models are the PlasticDB inputs, not equilibrated solvated starting systems.
Source-model paths are explicitly indexed. Additional structure copies used by
the original protocol or figure code retain their recorded paths and hashes.

## Files

- `metadata/protein_sources.csv`: source-model identity and relative package path.
  The original files are unchanged. B-column summaries describe the supplied
  files; an all-zero B column is not interpreted as pLDDT zero or as confidence
  metadata recovered from the prediction pipeline.
- `metadata/protein_preparation.csv`: actual termini, histidine states, charges,
  disulfide connectivity and sequence extensions recorded in the prepared
  systems and checked against the production TPRs in the source audit.
- `metadata/simulation_systems.csv`: one row per enzyme/length/presentation
  condition, including molecule/atom counts, initial box size and ion seed.
- `metadata/simulation_seeds.csv`: one row per production trajectory. NVT
  velocity-generation seeds and production stochastic seeds are distinct fields.
- `metadata/production_pme_tuning.csv`: nominal input settings and final logged
  tuning choices. Source-log hashes are retained; local machine paths are omitted.

Histidine `delta`, `epsilon` and `both` identify N-delta1 protonation,
N-epsilon2 protonation and double protonation, respectively. Residue indices are
simulation indices. These records describe the simulations actually performed;
they do not certify the scientific appropriateness of the chemical states or
replace the manuscript's methods and interpretation.

Run from the package root:

```bash
python -B 01_simulation/verify_preparation.py
```

This checks identity/coverage and the recorded NVT seeds against the packaged
MDPs. All 126 production TPRs are included and indexed separately in
`configurations/stages.csv`. Source-log hashes remain provenance identifiers;
the original logs and trajectories are not included. Readable topology inputs
and force-field includes are supplied. Starting coordinates can also be
extracted from each production TPR, as shown in the 01 README. These records
support inspection and production-stage replay, not a claim that every
historical upstream stage or scientific validation has been repeated.
