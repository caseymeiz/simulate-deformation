# Liver deformation batches

Generate neo-Hookean liver deformations with Elmer from tetrahedral meshes.

## Setup

The CPU environment includes Python 3.12, NumPy, VTK, and Elmer 26.2.1.
Supply your own mesh dataset.

Install [Pixi](https://pixi.sh), then from this repository:

```bash
pixi install --locked
pixi run install-elmer
pixi run test
```

The installer builds the pinned Elmer commit
`a19504ac53ec222e3355e182b08f2ff280c2203a` into `.elmer/`, using eight build jobs
by default (`ELMER_BUILD_JOBS` overrides this), with the Pixi C/C++/Fortran
compilers and the `ElasticSolve` and `ResultOutputSolve` modules used by this project. The first installation downloads
Conda packages and Elmer source. `pixi.lock` fixes the environment versions.

## Batch presets

Run from the repository root:

```bash
# Liver 18: 1,200 attempts.
pixi run liver18-1200

# All 131 livers: 20 attempts each, 2,620 total.
pixi run all131-20
```

The wrapper defaults to `all131_20`. Each invocation creates a new UTC-timestamped
`batch_...` directory under
`/shared/users/cjm6121/hybrid/MSD-RIT/Task03_Liver_NeoHookean/`.
Inputs come from `/shared/users/cjm6121/hybrid/MSD-RIT/Task03_Liver`.

| Setting | `liver18_1200` | `all131_20` |
|---|---|---|
| Mesh pattern | `liver_18_hcp_tet4.vtk` | `*_hcp_tet4.vtk` |
| Attempts per liver | 1,200 | 20 |
| Batch seed | 20261004 | 20261015 |
| Relaxation factor | 0.1 | 0.05 |
| Maximum nonlinear iterations | 100 | 200 |

Both presets use CG/Cholesky2 with fixed nonlinear relaxation.
Failed attempts count toward the requested total and are not retried or replaced.

Both presets use 8 workers, input scale 0.001 (mm to m), force factor 3,
Poisson ratio 0.35, 2–3 force patches, fixed-area scale 2, fixed margin 0.02 m,
and maximum force-patch radius 0.07 m. Young's modulus is sampled from 500–5000 Pa.
The all-liver preset processes every matching mesh; reproducing 131 livers requires
the original input set. Meshes are sorted by path, and the batch seed determines
each attempt's seed. Explicit CLI options override preset values.

To choose input and output directories:

```bash
pixi run python \
  batch.py /path/to/meshes /path/to/new-batch \
  --preset all131_20 --workers 8
```

Historical launch records are available in `docs/launches/`; their paths refer
to the original repository.

## Files

| File | Purpose |
|---|---|
| `generate_msd.sh` | Shared-dataset launcher with a fresh output directory |
| `batch.py` | Presets, job seeds, worker pool, rigid TRE and reports |
| `main.py` | Prepare and solve one deformation; validate and export results |
| `boundaries.py` | Connected fixed/loading patches and clearance |
| `elmer_mesh.py` | Native Elmer mesh writer |
| `mesh_utils.py` | Surface extraction, geometry checks and VTK output |
| `case.sif.template` | CG/Cholesky2 neo-Hookean solver settings |
| `boundary_to_vtk.py` | Optional labeled surface export for inspection |
| `test_*.py` | Sampling, batch and solver checks |
| `pixi.toml`, `pixi.lock` | Environment and task definitions |
| `scripts/install_elmer.sh` | Pinned local Elmer build |
| `docs/` | Historical method write-up and PDF; solver settings and commands differ from this project |

The target-area and force distributions derive from Volume2SurfaceCNN;
`VOLUME2SURFACECNN_LICENSE.txt` retains its license.

## Outputs

`batch_manifest.json` records settings, source meshes and seeds. Each attempt saves
`reference.vtk` in metres, a `simulation/` directory containing the native Elmer
mesh, `case.sif`, `solver.log`, and `manifest.json`, plus `metrics.json`.
Successful cases also save `simulation/deformed.vtk` and `rigid_registration.json`.
`results.csv` and `summary.json` update as attempts finish, including failures,
solve times and rigid TRE in millimetres. TRE uses all corresponding volume nodes
after a least-squares rotation and translation, without scaling or reflection.

Exported coordinates are already deformed. The exporter preserves input mesh
arrays and adds a displacement array. Success requires solver convergence,
finite geometry and no inverted tetrahedra.

## Single case and inspection

`main.py` accepts a linear tetrahedral VTK/VTU mesh already in metres and a new
output directory. For example, using the single-liver batch's solver settings:

```bash
pixi run python \
  main.py /path/to/liver-in-metres.vtk /tmp/liver-case \
  --seed 749935133 --force-factor 3 --min-force-patches 2 \
  --fixed-area-scale 2 --fixed-margin 0.02 --force-radius-max 0.07 \
  --relaxation-factor 0.1 --nonlinear-max-iterations 100 --prepare-only
```

Omit `--prepare-only` to solve. Preparation preserves the input tetrahedra and
vertices. It samples connected surface patches without changing their target
areas when growth fails. Fixed and force patches have the requested clearance;
force patches cannot share triangles. Forces become uniform Cartesian tractions
by dividing each force vector by its patch area.

For a labeled surface, run `boundary_to_vtk.py CASE/mesh.boundary surface.vtk`
in the same Pixi environment. Labels are 1 = free, 2 = fixed, 3+ = force patches.

## Checks

```bash
pixi run test
```
