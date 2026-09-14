# PETase-MD-2026-review

## PET Chain Extension Remodels Protein–Polymer Interfaces across Natural PET Hydrolases

Xiaomei Deng and Christian D. Lorenz

![Graphical abstract showing distinct interfaces, enzyme-specific interface remodeling, and interaction geometry](03_figures/toc/TOC_graphic.png)

## Abstract

Enzymatic degradation of poly(ethylene terephthalate) (PET) requires PET
hydrolases to bind and accommodate polymer chains at their surfaces, but how
substrate engagement varies among enzymes and with polymer chain length
remains poorly understood. Here, we used all-atom molecular dynamics simulations
to compare six natural PET hydrolases interacting with oligomers containing 4,
10, or 20 repeat units. In total, we analyzed 126 independent 100-ns simulations,
of which 71 met a prespecified global retention threshold and entered
residue-level contact analysis. The retained ensembles displayed distinct
interface architectures: IsPETase, TfCut1, and LCC used combinations of residues
around the catalytic site, W-loop, and binding cleft; FoCut5a interactions were
concentrated around its flap and binding-loop regions; HiC predominantly used
cleft-lining residues; and PHL7 strongly engaged its substrate-binding subsites.
Increasing PET length from 10 to 20 repeat units did not simply increase
association. Instead, it redistributed contacts in an enzyme-specific manner,
with the clearest residue-level changes in TfCut1 and HiC. Longer PET chains
increased direct hydrogen-bond occupancies involving TfCut1 Thr64 and Ser67 and
HiC Thr166, and increased aromatic-ring proximity involving TfCut1 Phe210.
These results show that PET hydrolases differ both in their retention of PET
and in how their retained interfaces accommodate longer polymer chains. The
identified interaction patterns and residues provide testable hypotheses for
mutagenesis, reactive calculations, and experimental studies of PET recognition
and hydrolytic activity.

## Project structure

This repository provides the data and code for reviewing the reported analyses,
regenerating the figures and tables, and replaying
the 126 production MD inputs. All input data for these documented review
commands are stored here; no Hugging Face download or original project
directory is required. Software dependencies are installed separately.

The manuscript and Supporting Information are handled separately through the
journal submission system. This repository contains the supporting data and code,
not the manuscript source files or PDFs.

| Folder | Contents |
|---|---|
| [01_simulation](01_simulation/README.md) | Production MD TPRs, readable parameters and topologies, source structures, and protocol code |
| [02_postprocessing](02_postprocessing/README.md) | Delivered measurements, statistical calculations, and table-generation code |
| [03_figures](03_figures/README.md) | Figure-generation code, structural views, and the TOC graphic |

## Getting started

Install the [documented environment](environment/README.md), then run:

```bash
python -B verify.py --checksums-only
python -B verify.py
```

The first command checks every delivered file and registered dependency.
The second also recalculates the supported statistics and tables.
Each numbered folder has instructions for its own commands. Reproduction
outputs should be written to new directories outside this checkout.

## Documentation

- [Reproduction scope and data coverage](REPRODUCIBILITY.md)
- [File index and checksums](metadata/file_manifest.csv)
- [Trajectory and retained-population index](metadata/trajectory_index.csv)
- [Validation results](metadata/TEST_RESULTS.md)

The original full trajectories and the separate full-panel cluster archive
are not included. Numerical reproduction starts from the supplied measurements;
the included PDB snapshots support the displayed structural figures.
Production TPR replay is a new simulation, not recovery of the original frames.
See the scope guide before interpreting a reproduction check.

## Code and data

[GitHub — PETase-MD-2026-review](https://github.com/XMayDeng/PETase-MD-2026-review)

This repository is the single entry point for the review materials. If it is
private, reviewers require an access arrangement before they can use the link.
The included [code license](metadata/CODE_LICENSE.txt) is retained from the
source package; third-party structures, parameters and software retain their
own terms.

## Citation

If you use this project, please cite
[PETase-MD-2026-review](https://github.com/XMayDeng/PETase-MD-2026-review).
