# 01 — Simulation inputs and protocol

This review package contains 126 original production `MD.tpr` files covering
six enzymes, seven PET length/presentation conditions per enzyme and three
replicas per condition. The TPR files preserve the actual production starting
states, molecular topology and compiled run parameters.

## Inspect or replay a production input

Run from the repository root with GROMACS 2025.4 available:

```bash
python -B 01_simulation/run.py list
python -B 01_simulation/run.py check --output-dir ../tpr_check
python -B 01_simulation/run.py replay --name TfCut1_L10_bi/rep01 --steps 100 --output-dir ../md_smoke
```

Use `--gmx /path/to/gmx` for an explicit executable. Replay defaults to CPU
execution with four OpenMP threads; `--gpu 0` requests the named GPU.
`--steps 100` is a short execution check. To run every production step recorded
in one TPR, replace it with `--full-stage`. This is computationally expensive
and will generate a new trajectory. It is not a continuation from the original
final checkpoint and is not expected to recover the original frames exactly.

All output paths must be new and outside this repository. No bundled input is
overwritten. The checker verifies TPR hashes and GROMACS readability separately
from execution. Successful execution does not establish convergence or the
scientific suitability of a model.

To inspect a recorded input or extract its starting coordinates:

```bash
gmx dump -s 01_simulation/inputs/runs/TfCut1_L10_bi/rep01/MD.tpr
gmx editconf -f 01_simulation/inputs/runs/TfCut1_L10_bi/rep01/MD.tpr -o ../TfCut1_L10_bi_rep01_initial.gro
```

Choose a new output filename when extracting coordinates. These are production
starting coordinates, not the unprepared PlasticDB structures.

## Files

- `inputs/runs/<enzyme>_<length>_<presentation>/repNN/MD.tpr`: production
  input for each trajectory, indexed in [stages.csv](configurations/stages.csv).
- `configurations/mdp/`: 378 effective NVT/NPT/MD MDP files, with the
  human-readable mapping in [mdp_index.csv](configurations/mdp_index.csv).
- `inputs/preparation/`: the recorded EM and chain-relaxation MDP files.
  Preparation-stage TPRs are intentionally not included.
- `configurations/templates/` and `configurations/cases/`: reusable
  protocol settings. Fifteen case YAMLs are historical deployed copies; the
  three TfCut1 YAMLs retain their documented current-source fallback status.
- `inputs/proteins/`, `production_structures/`, and
  `reference_structures/`: source models and structures needed by the figures.
- `inputs/pets/`: PET inputs, including the frozen preequilibrated conformers
  used for docking. These are distinct from the omitted full-panel cluster archive.
- `inputs/systems/*/topology/`: readable protein/PET topologies and restraint
  definitions; their force-field include files are present in `inputs/forcefield/`.
- `inputs/forcefields/`: the source force-field templates for protocol inspection.

Protein asset keys map to enzymes as follows: 00083=IsPETase, 00057=TfCut1,
00062=LCC, 00075=FoCut5a, 00121=HiC and 00137=PHL7.

## Parameter filenames

For example, `configurations/mdp/TfCut1/L10/bi/rep01/` contains
`NVT.mdp`, `NPT.mdp`, and `MD.mdp`. The 126 replicas therefore have 378
stage-parameter files. L4 uses `single`; L10/L20 use `head`, `tail`, and `bi`
(the source labels are `head-side`, `tail-side`, and `bidirectional`).
`rep01`–`rep03` identify replicas. Checksums are in the indexes, not filenames.

## Preparation and system records

| File in `metadata/` | Contents |
|---|---|
| [protein_sources.csv](metadata/protein_sources.csv) | Six PlasticDB source-model identities and paths |
| [protein_preparation.csv](metadata/protein_preparation.csv) | Actual termini, histidine states, charges, disulfide connectivity, and sequence extensions, checked against production TPRs |
| [simulation_systems.csv](metadata/simulation_systems.csv) | 42 enzyme/length/presentation conditions, molecule and atom counts, initial boxes, and ion seeds |
| [simulation_seeds.csv](metadata/simulation_seeds.csv) | 126 replicas with separate NVT velocity-generation and production stochastic seeds |
| [production_pme_tuning.csv](metadata/production_pme_tuning.csv) | Nominal settings and final logged PME tuning choices |

Source models are not equilibrated solvated systems. An all-zero source-model
B column is not interpreted as pLDDT zero or recovered prediction confidence.
Histidine `delta`, `epsilon`, and `both` denote N-delta1, N-epsilon2, and double
protonation. Residue indices are simulation indices. Source-log hashes identify
the original records; the logs themselves are not included.

```bash
python -B 01_simulation/verify_preparation.py
```

This checks metadata coverage, model identity, and recorded NVT seeds against
the delivered MDPs. For TPR parsing or execution, use the commands above.

## Protocol source versus supported execution

`code/scripts/`, `code/pipelines/`, `code/profiles/` and `gmx_adapter.py`
retain the preparation/docking/MD protocol source for inspection. The supported
review entry above starts at an existing production TPR. This slim repository
does not deliver the complete historical docking/relaxation outputs, PET-only
trajectories, equilibration TPRs, solvated GRO archives or original checkpoints.
Those upstream stages are not a tested end-to-end reconstruction in this
repository. `build_system.py` requires an explicitly selected protein input
via `--protein` or its configuration file; there is no default protein.

For the reported statistical and figure reproduction, use 02 and 03. Those
entries do not need GROMACS or the full original trajectories.
