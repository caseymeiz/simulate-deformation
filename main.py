"""Generate one neo-Hookean deformation of a tetrahedral VTK/VTU mesh using Elmer."""

# Assumption: input coordinates and output coordinates/displacements are in metres.
# No unit conversion is performed.

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import random
import re
from string import Template
import subprocess
import sys
import time

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from boundaries import FIXED_AREA_SCALE, FIXED_MARGIN, FORCE_RADIUS_MAX, PatchSizeError, sample_boundaries
from elmer_mesh import write_boundary, write_volume
from mesh_utils import boundary_faces, check_result_geometry, write_legacy


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path, help="Linear tetrahedral .vtk or .vtu mesh in metres")
    parser.add_argument("output", type=Path, help="New directory for this simulation")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--force-factor", type=float, default=1.5)
    parser.add_argument("--min-force-patches", type=int, choices=[1, 2, 3], default=1)
    parser.add_argument("--youngs-modulus", type=float, help="Pa; otherwise sampled from 500–5000 Pa")
    parser.add_argument("--poisson-ratio", type=float, default=0.35)
    parser.add_argument("--fixed-area-scale", type=float, default=FIXED_AREA_SCALE)
    parser.add_argument("--force-radius-max", type=float, default=FORCE_RADIUS_MAX,
                        help="Maximum sampled force-patch radius in metres; minimum is 0.01")
    parser.add_argument("--fixed-margin", type=float, default=FIXED_MARGIN,
                        help="Fixed-to-force patch clearance in metres; 0 disables clearance")
    parser.add_argument("--relaxation-factor", type=float, default=.1,
                        help="Multiplier for each nonlinear displacement update")
    parser.add_argument("--nonlinear-max-iterations", type=int, default=100)
    parser.add_argument("--prepare-only", action="store_true", help="Write inputs without running Elmer")
    parser.add_argument("--solver", default="ElmerSolver", help="Executable name or absolute path")
    return parser.parse_args(argv)


def read_mesh(path):
    readers = {".vtk": vtk.vtkUnstructuredGridReader, ".vtu": vtk.vtkXMLUnstructuredGridReader}
    reader = readers[path.suffix.lower()]()
    reader.SetFileName(str(path))
    reader.Update()
    mesh = reader.GetOutput()
    types = vtk_to_numpy(mesh.GetCellTypes())
    offsets = vtk_to_numpy(mesh.GetCells().GetOffsetsArray())[:-1][types == vtk.VTK_TETRA]
    connectivity = vtk_to_numpy(mesh.GetCells().GetConnectivityArray())
    tetra = connectivity[offsets[:, None] + np.arange(4)]
    return mesh, tetra


def prepare(args, source, tetra):
    points = vtk_to_numpy(source.GetPoints().GetData()).astype(float)
    faces, parents = boundary_faces(tetra)
    write_volume(args.output, points, tetra, len(faces))
    rng = random.Random(args.seed)
    labels, fixed, tractions = sample_boundaries(points, faces, rng, args.force_factor,
                                                args.min_force_patches, args.fixed_area_scale,
                                                args.fixed_margin, args.force_radius_max)
    write_boundary(args.output, faces, parents, labels)
    modulus = args.youngs_modulus if args.youngs_modulus is not None else rng.uniform(500, 5000)
    settings = {key: value for key, value in vars(args).items() if not isinstance(value, Path)}
    settings["youngs_modulus"] = round(modulus, 6)
    blocks = []
    for load in tractions:
        index, traction = load["boundary_id"], load["traction_pa"]
        blocks.append(f'''Boundary Condition {index}
  Target Boundaries(1) = {index}
  Name = "Load {index - 2}"
  ! Cartesian traction components in Pa (force per unit area).
  Force 1 = {traction[0]:.6f}
  Force 2 = {traction[1]:.6f}
  Force 3 = {traction[2]:.6f}
End''')
    template = Template(Path(__file__).with_name("case.sif.template").read_text())
    sif = template.substitute(settings, force_boundaries="\n\n".join(blocks))
    (args.output / "case.sif").write_text(sif)
    return settings, fixed, tractions


def save_result(args, source, tetra):
    result_path = args.output / "result/case_t0001.vtu"
    geometry = check_result_geometry(result_path, tetra)
    if not geometry["finite_geometry"] or geometry["inverted_tetrahedra"]:
        return {"status": "invalid_geometry", **geometry}
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(result_path))
    reader.Update()
    result = reader.GetOutput()
    target = vtk_to_numpy(result.GetPoints().GetData()).astype(float)
    displacement = vtk_to_numpy(result.GetPointData().GetArray("displacement")).astype(float)
    output = vtk.vtkUnstructuredGrid()
    output.DeepCopy(source)
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(target, deep=True))
    output.SetPoints(points)
    vectors = numpy_to_vtk(displacement, deep=True)
    vectors.SetName("displacement")
    output.GetPointData().AddArray(vectors)
    write_legacy(args.output / "deformed.vtk", output)
    return {"status": "success", **geometry}


def solve(args, source, tetra, manifest):
    case_sif = manifest["case_sif"]
    nonlinear_tolerance = float(re.search(
        r"^\s*Nonlinear System Convergence Tolerance\s*=\s*(\S+)", case_sif, re.M)[1])
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    print(f"Elmer log: {args.output / 'solver.log'}", flush=True)
    start = time.monotonic()
    with (args.output / "solver.log").open("w") as log:
        process = subprocess.run([args.solver, "case.sif"], cwd=args.output, env=env,
                                 stdout=log, stderr=subprocess.STDOUT)
    elapsed = time.monotonic() - start
    log = (args.output / "solver.log").read_text()
    changes = re.findall(r"NS \(ITER=\d+\) \(NRM,RELC\): \(\s*\S+\s+(\S+)", log)
    change = float(changes[-1]) if changes else None
    errors = [line for line in log.splitlines() if "ERROR::" in line]
    outcome = dict(status="solver_failed", solver_returncode=process.returncode,
                   nonlinear_relative_change=change, errors=errors)
    if process.returncode == 0 and not errors and (args.output / "result/case_t0001.vtu").exists():
        outcome["status"] = "not_converged"
        if change is not None and np.isfinite(change) and change < nonlinear_tolerance:
            outcome.update(save_result(args, source, tetra))
    manifest.update(outcome, elapsed_seconds=elapsed)
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main():
    args = arguments()
    source, tetra = read_mesh(args.mesh)
    args.output.mkdir(parents=True)
    (args.output / "result").mkdir()
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(), "argv": sys.argv[1:],
        "input_mesh": str(args.mesh.resolve()), "output": str(args.output.resolve()),
        "settings": {key: value for key, value in vars(args).items() if not isinstance(value, Path)},
        "input_mesh_units": "m", "simulation_units": "m", "deformed_mesh_units": "m",
        "status": "prepared",
    }
    manifest_path = args.output / "manifest.json"
    try:
        settings, fixed, tractions = prepare(args, source, tetra)
    except PatchSizeError as error:
        manifest.update(status="preparation_failed", error=str(error))
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        raise SystemExit(str(error)) from error
    case_sif = (args.output / "case.sif").read_text()
    manifest.update(settings=settings, fixed_boundary=fixed, tractions=tractions, case_sif=case_sif)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    if args.prepare_only:
        print(f"Prepared {args.output / 'case.sif'}")
        return
    solve(args, source, tetra, manifest)
    print(json.dumps(manifest, indent=2))
    if manifest["status"] != "success":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
