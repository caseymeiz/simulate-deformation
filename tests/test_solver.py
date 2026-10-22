import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from main import solve


DIVERGENCE_ERROR = "ERROR:: IterSolve: Numerical Error: System diverged over maximum tolerance."


class SolverTests(unittest.TestCase):
    def run_case(self, outcome):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        output = Path(temporary.name)
        (output / "result").mkdir()
        (output / "mesh.nodes").write_text("unchanged mesh")
        args = SimpleNamespace(output=output, relaxation_factor=.2,
                               solver="ElmerSolver")
        manifest = {"case_sif": "  Nonlinear System Relaxation Factor = 0.2\n"
                    "  Nonlinear System Convergence Tolerance = 1e-4\n"
                    "  Nonlinear System Max Iterations = 200\n"}
        (output / "case.sif").write_text(manifest["case_sif"])
        factors = []

        def run(command, cwd, env, stdout, stderr):
            self.assertEqual(command, ["ElmerSolver", "case.sif"])
            self.assertEqual(list((cwd / "result").iterdir()), [])
            self.assertEqual((cwd / "mesh.nodes").read_text(), "unchanged mesh")
            sif = (cwd / "case.sif").read_text()
            self.assertIn("Nonlinear System Max Iterations = 200", sif)
            factors.append(float(re.search(r"Relaxation Factor = (\S+)", sif)[1]))
            kind = outcome
            if kind == "divergence":
                stdout.write(DIVERGENCE_ERROR + "\n")
            elif kind == "iteration_limit":
                stdout.write("ERROR:: IterSolve: Numerical Error: Too many iterations were needed.\n")
            else:
                change = "0.1" if kind == "not_converged" else "1e-5"
                stdout.write(f"NS (ITER=200) (NRM,RELC): ( 1.0 {change} )\n")
            # Even a failed solver may leave an output file.
            (cwd / "result/case_t0001.vtu").write_text(kind)
            return SimpleNamespace(returncode=1 if kind in ("divergence", "iteration_limit") else 0)

        def save(*_):
            kind = outcome
            return {"status": "invalid_geometry" if kind == "invalid_geometry" else "success"}

        with patch("main.subprocess.run", side_effect=run), \
             patch("main.save_result", side_effect=save), \
             patch("main.time.monotonic", side_effect=range(20)):
            solve(args, None, None, manifest)
        self.assertEqual(json.loads((output / "manifest.json").read_text()), manifest)
        self.assertEqual(manifest["elapsed_seconds"], len(factors))
        return output, manifest, factors

    def test_one_solve_and_recorded_outcome(self):
        for outcome in ("success", "divergence", "iteration_limit", "not_converged", "invalid_geometry"):
            with self.subTest(outcome=outcome):
                output, manifest, factors = self.run_case(outcome)
                self.assertEqual(factors, [.2])
                self.assertEqual(manifest["final_relaxation_factor"], .2)
                self.assertEqual(len(manifest["attempts"]), 1)
                self.assertEqual(manifest["status"],
                                 "solver_failed" if outcome in ("divergence", "iteration_limit") else outcome)
                self.assertFalse((output / "attempts").exists())
