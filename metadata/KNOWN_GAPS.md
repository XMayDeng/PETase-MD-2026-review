# Review scope and remaining decisions

The included-file manifest and [reproduction guide](../REPRODUCIBILITY.md)
define this review export. Its documented numerical, figure and production-TPR
commands do not require the separately prepared HF package.

- Full original trajectories, the full-panel cluster archive, earlier-stage
  TPRs, original final checkpoints and historical docking output archives are
  intentionally excluded. New atom-level extraction, altered selections or
  reclustering of original trajectories require additional original data.
- Numerical reproduction starts from delivered measurements. Curated mappings,
  initial screening and structural measurements remain source inputs where
  specified in `TABLES_AND_STRUCTURAL_SENSITIVITY.md`.
- Production replay starts from the stored TPR state. It is not an end-to-end
  rerun of preparation, docking and equilibration, nor recovery of original
  trajectory frames. Bounded execution tests are not 126 complete 100 ns reruns.
- Fifteen case YAMLs are historical deployed copies. The three TfCut1 cases
  retain their explicit current-source fallback provenance.
- The recorded chemical preparation states are preserved in the TPRs and
  metadata. Execution tests do not establish scientific appropriateness or
  sampling convergence.
- Exact PNG equality uses the recorded Python, fonts, ChimeraX and graphics
  environment; it is not asserted for arbitrary operating systems or drivers.
- The copied manuscript's acknowledgments and Data and Software Availability
  section require author-approved completion. Repository packaging does not
  modify the canonical manuscript.
- A private repository needs reviewer-access arrangements. Long-term archival
  identifiers, any DOI, and third-party redistribution permissions remain
  author decisions. The original code license is retained, not expanded to
  cover third-party material.
