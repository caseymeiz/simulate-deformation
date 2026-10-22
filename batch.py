"""Run independent random liver deformations and measure TRE after rigid fitting."""

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import multiprocessing
import os
from pathlib import Path
import random
import subprocess
import sys
import time

import numpy as np
import vtk
from tqdm import tqdm
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from main import FIXED_AREA_SCALE, FIXED_MARGIN, FORCE_RADIUS_MAX, read_mesh
from mesh_utils import write_legacy


PRESETS = {
    "liver18_1200": dict(pattern="liver_18_hcp_tet4.vtk", deformations_per_liver=1200,
                         seed=20261004, relaxation_factor=.1, nonlinear_max_iterations=100),
    "all131_20": dict(pattern="*_hcp_tet4.vtk", deformations_per_liver=20,
                      seed=20261015, relaxation_factor=.05, nonlinear_max_iterations=200),
}


def rigid_metrics(reference, deformed):
    """Least-squares proper rigid fit over corresponding volume nodes; input metres."""
    x, y = reference * 1000, deformed * 1000
    xc, yc = x.mean(axis=0), y.mean(axis=0)
    u, _, vt = np.linalg.svd((x - xc).T @ (y - yc))
    correction = np.diag([1., 1., np.linalg.det(vt.T @ u.T)])
    rotation = vt.T @ correction @ u.T
    translation = yc - xc @ rotation.T
    distances = np.linalg.norm(y - (x @ rotation.T + translation), axis=1)
    return {
        "mean_tre_mm": float(distances.mean()),
        "rms_tre_mm": float(np.sqrt(np.mean(distances ** 2))),
        "max_tre_mm": float(distances.max()),
        "rotation": rotation.tolist(), "translation_mm": translation.tolist(),
        "number_of_nodes": len(x),
    }


def run_sample(job, args):
    start = time.monotonic()
    folder = args.output / job["sample"]
    folder.mkdir()
    source, _ = read_mesh(Path(job["input_mesh"]))
    reference = vtk_to_numpy(source.GetPoints().GetData()).astype(float) * args.input_scale
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(reference, deep=True))
    source.SetPoints(points)
    mesh_path = folder / "reference.vtk"
    write_legacy(mesh_path, source)
    case = folder / "simulation"
    command = [sys.executable, str(Path(__file__).with_name("main.py")),
               str(mesh_path), str(case), "--seed", str(job["seed"]),
               "--force-factor", str(args.force_factor),
               "--poisson-ratio", str(args.poisson_ratio),
               "--relaxation-factor", str(args.relaxation_factor),
               "--nonlinear-max-iterations", str(args.nonlinear_max_iterations),
               "--fixed-area-scale", str(args.fixed_area_scale),
               "--force-radius-max", str(args.force_radius_max),
               "--fixed-margin", str(args.fixed_margin),
               "--min-force-patches", str(args.min_force_patches)]
    with (folder / "run.log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    row = {**job, "status": "invocation_failed", "returncode": result.returncode,
           "failure_reason": "invocation_failed",
           "force_factor": args.force_factor, "poisson_ratio": args.poisson_ratio,
           "relaxation_factor": args.relaxation_factor,
           "nonlinear_max_iterations": args.nonlinear_max_iterations,
           "fixed_area_scale": args.fixed_area_scale,
           "fixed_margin": args.fixed_margin,
           "force_radius_max": args.force_radius_max,
           "youngs_modulus_pa": None,
           "force_patches": None, "solve_seconds": None, "sample_seconds": None,
           "mean_tre_mm": None, "rms_tre_mm": None, "max_tre_mm": None}
    manifest_path = case / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        row.update(status=manifest["status"],
                   youngs_modulus_pa=manifest["settings"].get("youngs_modulus"),
                   force_patches=len(manifest["tractions"]) if "tractions" in manifest else None,
                   solve_seconds=manifest.get("elapsed_seconds"))
        row["failure_reason"] = None if row["status"] == "success" else (
            manifest["errors"][0].strip() if manifest.get("errors") else
            manifest.get("error", row["status"]))
    if row["status"] == "success":
        deformed, _ = read_mesh(case / "deformed.vtk")
        metrics = rigid_metrics(reference, vtk_to_numpy(deformed.GetPoints().GetData()))
        (folder / "rigid_registration.json").write_text(json.dumps(metrics, indent=2) + "\n")
        row.update({key: metrics[key] for key in ("mean_tre_mm", "rms_tre_mm", "max_tre_mm")})
    row["sample_seconds"] = time.monotonic() - start
    (folder / "metrics.json").write_text(json.dumps(row, indent=2) + "\n")
    return row


def save_report(output, rows, elapsed):
    rows = sorted(rows, key=lambda row: row["sample"])
    successful = [row for row in rows if row["status"] == "success"]
    solved = [row for row in rows if row["solve_seconds"] is not None]
    def mean(group, key):
        return float(np.mean([row[key] for row in group])) if group else None
    summary = {
        "completed_samples": len(rows), "successful_samples": len(successful),
        "total_solves": len(solved),
        "failed_samples": len(rows) - len(successful),
        "status_counts": dict(Counter(row["status"] for row in rows)),
        "failure_reason_counts": dict(Counter(
            row["failure_reason"] for row in rows if row["status"] != "success")),
        "mean_successful_solve_seconds": mean(successful, "solve_seconds"),
        "mean_successful_sample_seconds": mean(successful, "sample_seconds"),
        "mean_all_solve_seconds": mean(solved, "solve_seconds"),
        "mean_tre_mm": mean(successful, "mean_tre_mm"),
        "batch_wall_seconds": elapsed,
        "tre_definition": "Equal-weight mean of successful samples' mean Euclidean node errors after least-squares rigid fitting; all corresponding volume nodes; no scaling or reflection.",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output / "results.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return summary


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Directory of tetrahedral liver meshes")
    parser.add_argument("output", type=Path, help="New batch directory")
    parser.add_argument("--preset", choices=PRESETS, help="Settings and seed from a retained batch")
    parser.add_argument("--pattern", default="*_hcp_tet4.vtk")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--samples", type=int,
                           help="Number of randomly selected livers (default: 20)")
    selection.add_argument("--deformations-per-liver", type=int,
                           help="Generate this many deformations for every matching liver")
    parser.add_argument("--seed", type=int, default=20260916)
    parser.add_argument("--force-factor", type=float, default=1.5)
    parser.add_argument("--poisson-ratio", type=float, default=0.35)
    parser.add_argument("--fixed-area-scale", type=float, default=FIXED_AREA_SCALE)
    parser.add_argument("--force-radius-max", type=float, default=FORCE_RADIUS_MAX,
                        help="Maximum sampled force-patch radius in metres; minimum is 0.01")
    parser.add_argument("--fixed-margin", type=float, default=FIXED_MARGIN,
                        help="Fixed-to-force patch clearance in metres; 0 disables clearance")
    parser.add_argument("--relaxation-factor", type=float, default=.1,
                        help="Multiplier for each nonlinear displacement update")
    parser.add_argument("--nonlinear-max-iterations", type=int, default=100)
    parser.add_argument("--min-force-patches", type=int, choices=[1, 2, 3], default=2)
    parser.add_argument("--input-scale", type=float, default=1., help="Coordinate multiplier to metres; use 0.001 for mm")
    parser.add_argument("--workers", type=int, default=8)
    selected, _ = parser.parse_known_args(argv)
    if selected.preset:
        parser.set_defaults(**PRESETS[selected.preset], input_scale=.001, workers=8,
                            force_factor=3., poisson_ratio=.35, min_force_patches=2,
                            fixed_area_scale=2., fixed_margin=.02, force_radius_max=.07)
    args = parser.parse_args(argv)
    if selected.preset and args.samples is not None:
        parser.error("--samples cannot be used with --preset; use --deformations-per-liver")
    if args.samples is None and args.deformations_per_liver is None:
        args.samples = 20
    return args


def make_jobs(args):
    """Preserve sorted mesh order and the original per-sample seed sequence."""
    rng = random.Random(args.seed)
    meshes = sorted(args.input.glob(args.pattern))
    if args.deformations_per_liver is None:
        meshes = rng.sample(meshes, args.samples)
    else:
        args.samples = None
        meshes = [mesh for mesh in meshes for _ in range(args.deformations_per_liver)]
    return [{"sample": f"{index:03d}_{mesh.stem}", "input_mesh": str(mesh),
             "seed": rng.randrange(2 ** 31)} for index, mesh in enumerate(meshes, 1)]


def main():
    args = arguments()
    args.input, args.output = args.input.resolve(), args.output.resolve()
    args.output.mkdir(parents=True)
    jobs = make_jobs(args)
    manifest = {"argv": sys.argv[1:], "settings": {k: str(v) if isinstance(v, Path) else v
                for k, v in vars(args).items()}, "jobs": jobs}
    (args.output / "batch_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = "1"
    start, rows = time.monotonic(), []
    with ProcessPoolExecutor(max_workers=args.workers,
                             mp_context=multiprocessing.get_context("spawn")) as pool:
        futures = [pool.submit(run_sample, job, args) for job in jobs]
        with tqdm(total=len(jobs), desc="Deformations", unit="sample", smoothing=0) as progress:
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                summary = save_report(args.output, rows, time.monotonic() - start)
                progress.set_postfix(success=summary["successful_samples"],
                                     failed=summary["failed_samples"], refresh=False)
                progress.update(1)
                tqdm.write(f"[{len(rows)}/{len(jobs)}] {json.dumps(row)}")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
