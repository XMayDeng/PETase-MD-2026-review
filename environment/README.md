# Reproduction environment

## Numerical analysis and figure generation

Use Python 3.12.3 and the exact versions in [requirements.txt](requirements.txt).
The reference environment is recorded in [versions.json](versions.json).
Install dependencies in a separately managed environment:

```bash
python -m pip install -r environment/requirements.txt
```

Run without `-O` or `-OO`; validation assertions must remain enabled.
The renderer checks the exact Matplotlib, NumPy, pandas, SciPy and Pillow
versions and bundled DejaVu Sans fonts. Do not change these versions simply
to bypass a failed reference comparison.

Figures 1, 5, 6 and S2 require UCSF ChimeraX
`1.13.dev202606262139`. Supply its executable with `--chimerax`.
It is installed separately and is not auto-downloaded. The reference uses the
Linux technical-preview build with a working offscreen OpenGL stack. Scenes
are rendered from local PDB files in fresh temporary sessions. Existing PNGs
are reference outputs, not substitutes for molecular rendering.

Exact PNG identity is tied to the recorded environment and graphics stack.
Cross-platform identity is not assumed. Numerical comparison and the cause of
any graphical mismatch should be inspected separately.

## Production TPR replay

GROMACS 2025.4 is the tested version. The original build report is in
[gromacs_build.txt](gromacs_build.txt): external MPI with OpenMP, mixed precision,
Open MPI 4.1.6, and CUDA 12.8 support. CPU replay is available; CUDA hardware is
only needed when explicitly requesting GPU execution. A software installation
and drivers are not part of this repository.

GROMACS is not needed for the supplied-measurement statistical calculations or
for the pure-Python TOC schematic. LaTeX is not required by the documented
review commands.

## Archived upstream environment descriptions

[docking-linux-64.lock](docking-linux-64.lock), [docking_packages.csv](docking_packages.csv),
[simulation_versions.json](simulation_versions.json), and
[trajectory_requirements.txt](trajectory_requirements.txt) document the original
protocol/extraction environments. They are retained for provenance, not
additional installation requirements for the supported review commands.

The docking package is Vina 1.2.6 (build `py310hea6c23e_2`); its executable
reports `AutoDock Vina 7ac2999-mod`. The historical extraction environment's
NumPy 2.4.4 pin is distinct from the renderer's NumPy 2.4.3 pin. Do not upgrade
the figure environment to match that archived extraction environment.
