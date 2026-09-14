# Review export validation

Checked on 2026-09-14 in the recorded Linux environment, using the independent
`PETase-MD-2026-review` export. Validation outputs were written outside the
repository. The original simulation directories and the concurrently maintained
HF package were not modified.

| Check | Result | Scope |
|---|---|---|
| File hashes and dependency registries | PASS | Every included file; analysis/figure/table inputs; 126 production TPR identities |
| Preparation metadata | PASS | Six enzymes, 42 conditions, 126 trajectories, 378 effective MDPs |
| GROMACS TPR parsing | PASS | All 126 production TPRs |
| Bounded production execution | PASS | Six runs, one L10/bi/rep01 input per enzyme, 100 steps each, CPU |
| Initial-coordinate extraction | PASS | One recorded TfCut1 production input, 214,972 atoms |
| Topology include closure | PASS | 42 system topologies and their 212-file include closure |
| Numerical recalculation | PASS | Eight analysis-entry steps and 15 registered numerical comparisons |
| Tables and supplementary data | PASS | Six TeX tables and Data S1/S2 match the delivered references byte for byte |
| Numbered figure reconstruction | PASS | All 15 PNGs match their reference hashes |
| TOC reconstruction | PASS | Independently redrawn PNG matches the delivered TOC exactly |
| CLI rejection tests | PASS | Eight invalid/unsafe requests rejected without overwriting package inputs |
| Source syntax and documentation links | PASS | 78 Python files parsed and 39 local links checked |
| Root verification entry | PASS | `python -B verify.py`, including numerical recalculation and exact regenerated-table checks |

Machine-readable details, numerical tolerances, output hashes and tested input
identities are in [validation_results.json](validation_results.json).
The final data/code verification repeated file-integrity checks, numerical
recalculation, table reconstruction, all 126 TPR reads, all 15 numbered figure
renders, the TOC render, Python syntax checks and documentation-link checks.
The bounded MD runs, coordinate extraction, topology closure and CLI rejection
tests refer to the recorded tests on the same unchanged inputs and code.
Manuscript and SI sources/PDFs are not included. Their compilation and page
checks are outside the scope of this data/code repository.

## Repeating the checks

From the repository root after installing the documented dependencies:

```bash
python -B verify.py --checksums-only
python -B verify.py
python -B 01_simulation/run.py check --output-dir ../tpr_check
python -B 03_figures/run.py --output-dir ../figures_and_tables --chimerax /path/to/ChimeraX
python -B 03_figures/toc/make_figure.py --output ../toc_result/TOC_graphic.png
```

## What these results do not claim

The six short MD runs are execution checks, not six new scientific trajectories
or 126 full 100 ns reruns. The numerical checks begin at the supplied
measurements; they do not repeat raw-trajectory extraction or initial candidate
selection. No full-panel cluster archive is needed by these analysis/figure
commands. The protocol source for upstream preparation is available for
inspection, but an end-to-end upstream workflow is not the supported review
entry. Exact graphical equality is tied to the recorded software/graphics stack.

Submission fields and reviewer access remain author-controlled. The repository
is kept private at the author's request. A private URL alone does not authorize
reviewer access, and repository verification does not constitute journal
submission or approval of the manuscript.
