# Human-readable file and folder names

Author-confirmed convention, 2026-09-13: names in the reader-facing package
should communicate what an item contains. Use scientific objects, conditions,
replicates and processing stages rather than opaque hashes or internal job IDs.
This convention applies to future packaging work as well as the MDP layout below.

## Current simulation-parameter layout

```text
01_simulation/configurations/mdp/
  TfCut1/
    L10/
      bi/
        rep01/
          NVT.mdp
          NPT.mdp
          MD.mdp
```

- Enzyme folders use the established labels: IsPETase, TfCut1, LCC, FoCut5a,
  HiC and PHL7.
- PET-length folders are L4, L10 and L20.
- L4 uses `single`. L10/L20 use `head`, `tail` and `bi`, mapped respectively
  from `head-side`, `tail-side` and `bidirectional` in the source records.
- `rep01`, `rep02` and `rep03` identify the three replicas within a condition.
- `NVT.mdp`, `NPT.mdp` and `MD.mdp` contain the respective equilibration or
  production-stage effective parameters. Active parameter lines are unchanged;
  machine-specific generated comment headers are omitted in the delivery copy.
  There are 126 replica folders and 378 MDP files.

The authoritative crosswalk is
[mdp_index.csv](../01_simulation/configurations/mdp_index.csv). Its original
`trajectory_id` and `stage` remain unchanged, alongside the readable relative
`mdp_path` and content `sha256`. Names are a navigation aid, not a substitute
for source identity or checksums. The executed-stage index is
`01_simulation/configurations/stages.csv`: in this review package it covers
the 126 production MD TPRs. Earlier-stage TPRs are not included.

## Rules for subsequent materials

1. Prefer short descriptive English names, established enzyme labels and clear
   directory hierarchy. Use zero-padded replica/cluster numbers where sorting
   matters. Explain study-specific abbreviations in the responsible README.
2. Give a file a functional name such as `topology.top`, `cluster_summary.csv`
   or `representatives.pdb` within a clearly identified system directory. Do not
   rename formats that require a standard filename without checking their tools.
3. Store full SHA-256 values in manifests and provenance tables, not filenames.
   Preserve original database IDs, simulation IDs and original names in source
   mappings. A database accession may remain a useful filename when it is the
   scientific identity, provided the associated enzyme is clearly indexed.
4. Avoid machine/user paths, opaque job timestamps, and ambiguous suffixes such
   as `final2`, `new` or `unpacked` in reader-facing names. Preserve actual dates
   or versions where they distinguish genuine releases.
5. Every rename must update active indexes, manifests and executable references,
   verify unchanged scientific content, and retain a reversible old/new mapping.
   Never rename original experimental files as a side effect of package cleanup.
6. Do not duplicate a data file merely to offer several names. Keep one maintained
   package copy and use the crosswalk for alternative scientific/source identities.

Prepared systems, stage TPRs and representative structures use readable condition
folders such as `TfCut1_L10_bi/rep01/`. Original database keys are retained only
where they identify source assets and are mapped to enzyme names in the 01 README.
All changes apply to independent package copies, never the source experiments.
