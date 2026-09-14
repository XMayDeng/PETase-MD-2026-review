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
- [Preparation metadata](PREPARATION_METADATA.md): model identity, chemical
  preparation, system composition, seeds and recorded PME tuning.

Protein asset keys map to enzymes as follows: 00083=IsPETase, 00057=TfCut1,
00062=LCC, 00075=FoCut5a, 00121=HiC and 00137=PHL7.

## Protocol source versus supported execution

`code/scripts/`, `code/pipelines/`, `code/profiles/` and `gmx_adapter.py`
retain the preparation/docking/MD protocol source for inspection. The supported
review entry above starts at an existing production TPR. This slim repository
does not deliver the complete historical docking/relaxation outputs, PET-only
trajectories, equilibration TPRs, solvated GRO archives or original checkpoints.
Do not interpret the presence of protocol code as a tested end-to-end
reconstruction of those omitted upstream stages.

For the reported statistical and figure reproduction, use 02 and 03. Those
entries do not need GROMACS or the full original trajectories.
