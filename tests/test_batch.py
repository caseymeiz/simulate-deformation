import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from batch import arguments, make_jobs, rigid_metrics, save_report


class ReportTests(unittest.TestCase):
    def test_failure_counts_and_successful_tre(self):
        reasons = [None, "line search failed", "line search failed", "not_converged"]
        rows = [{"sample": str(index),
                 "status": "success" if reason is None else "solver_failed",
                 "failure_reason": reason, "solve_seconds": 10.,
                 "solver_attempts": 2,
                 "sample_seconds": 11., "mean_tre_mm": 2. if reason is None else None}
                for index, reason in enumerate(reasons)]
        with TemporaryDirectory() as directory:
            summary = save_report(Path(directory), rows, 30.)
            saved = json.loads((Path(directory) / "summary.json").read_text())
        self.assertEqual(saved, summary)
        self.assertEqual(summary["completed_attempts"], 4)
        self.assertEqual(summary["total_solver_attempts"], 8)
        self.assertEqual(summary["successful_samples"], 1)
        self.assertEqual(summary["failed_samples"], 3)
        self.assertEqual(summary["failure_reason_counts"],
                         {"line search failed": 2, "not_converged": 1})
        self.assertEqual(summary["mean_tre_mm"], 2.)


class RigidTests(unittest.TestCase):
    def test_recovers_known_rigid_motion(self):
        points = np.random.default_rng(7).normal(size=(100, 3)) * .1
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        translation = np.array([.2, -.4, .3])
        result = rigid_metrics(points, points @ rotation.T + translation)
        self.assertLess(result["mean_tre_mm"], 1e-10)
        np.testing.assert_allclose(result["rotation"], rotation, atol=1e-14)
        np.testing.assert_allclose(result["translation_mm"], translation * 1000)

    def test_does_not_fit_reflection_or_scaling(self):
        points = np.array([[1., 1., 1.], [1., -1., -1.],
                           [-1., 1., -1.], [-1., -1., 1.]]) * .001
        reflected = rigid_metrics(points, points * [-1, 1, 1])
        self.assertAlmostEqual(np.linalg.det(reflected["rotation"]), 1.)
        self.assertGreater(reflected["mean_tre_mm"], 0.)
        scaled = rigid_metrics(points, 2 * points)
        self.assertAlmostEqual(scaled["mean_tre_mm"], np.sqrt(3))


class PresetTests(unittest.TestCase):
    def test_original_job_counts_and_seed_sequences(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for index in range(131):
                (root / f"liver_{index}_hcp_tet4.vtk").touch()
            for preset, count, first, last in (
                ("liver18_1200", 1200, 749935133, 1161346008),
                ("all131_20", 2620, 907289777, 142363792),
            ):
                with self.subTest(preset=preset):
                    args = arguments([directory, str(root / "output"), "--preset", preset])
                    jobs = make_jobs(args)
                    self.assertEqual(len(jobs), count)
                    self.assertEqual(jobs[0]["seed"], first)
                    self.assertEqual(jobs[-1]["seed"], last)
                    self.assertEqual(len({job["sample"] for job in jobs}), count)
                    counts = {}
                    for job in jobs:
                        counts[job["input_mesh"]] = counts.get(job["input_mesh"], 0) + 1
                    self.assertEqual(set(counts.values()), {args.deformations_per_liver})
                    self.assertEqual(args.input_scale, .001)
                    self.assertEqual(args.force_factor, 3.)
                    self.assertEqual(args.fixed_area_scale, 2.)
                    self.assertFalse(args.line_search)
                    self.assertEqual(args.divergence_retries, 0)
                    self.assertEqual(args.template, "cg")

    def test_explicit_options_override_preset(self):
        args = arguments(["input", "output", "--seed", "7", "--preset", "all131_20",
                          "--deformations-per-liver", "2", "--workers", "1"])
        self.assertEqual(args.seed, 7)
        self.assertEqual(args.deformations_per_liver, 2)
        self.assertEqual(args.workers, 1)
        self.assertEqual(args.relaxation_factor, .05)
        self.assertEqual(args.nonlinear_max_iterations, 200)
