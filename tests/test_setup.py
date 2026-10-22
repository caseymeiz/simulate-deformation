"""Run with: pixi run test."""

from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch as mock_patch

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from boundaries import PatchSizeError, grow_patch, outside_fixed_margin, sample_boundaries, triangle_neighbors
from elmer_mesh import write_boundary, write_volume
from mesh_utils import boundary_faces


class SetupTests(unittest.TestCase):
    def test_original_force_radius_range(self):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(.2)
        sphere.SetThetaResolution(32)
        sphere.SetPhiResolution(32)
        sphere.Update()
        mesh = sphere.GetOutput()
        points = vtk_to_numpy(mesh.GetPoints().GetData()).astype(float)
        faces = vtk_to_numpy(mesh.GetPolys().GetConnectivityArray()).reshape(-1, 3)
        radii = []
        for seed in range(5):
            labels, fixed, loads = sample_boundaries(
                points, faces, random.Random(seed), 3., min_force_patches=2,
                fixed_area_scale=1., fixed_margin=0., force_radius_max=.15)
            self.assertEqual(fixed["area_scale"], 1.)
            self.assertEqual(fixed["margin_m"], 0.)
            triangles = points[faces]
            areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                            triangles[:, 2] - triangles[:, 0]), axis=1) / 2
            for load in loads:
                radius = load["sampled_radius_m"]
                radii.append(radius)
                self.assertGreaterEqual(radius, .01)
                self.assertLess(radius, .15)
                self.assertAlmostEqual(load["target_area_m2"], np.pi * radius ** 2)
                self.assertAlmostEqual(areas[labels == load["boundary_id"]].sum(), load["area_m2"])
                self.assertAlmostEqual(np.linalg.norm(load["traction_pa"]) * load["area_m2"],
                                       load["force_magnitude_n"], places=6)
        self.assertGreater(max(radii), .07)

    def test_surface_parents_and_native_files(self):
        points = np.array([[0., 0., 0.], [0.123456789, 0., 0.], [0., 1., 0.],
                           [0., 0., 1.], [0., 0., -1.]])
        tetra = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
        faces, parents = boundary_faces(tetra)
        self.assertEqual(len(faces), 6)
        for face, parent in zip(faces, parents):
            self.assertTrue(set(face) <= set(tetra[parent]))
            self.assertNotEqual(set(face), {0, 1, 2})
        labels = np.array([1, 2, 3, 3, 2, 1])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            write_volume(path, points, tetra, len(faces))
            write_boundary(path, faces, parents, labels)
            np.testing.assert_array_equal(np.loadtxt(path / "mesh.nodes")[:, 2:], points)
            np.testing.assert_array_equal(np.loadtxt(path / "mesh.elements", dtype=int)[:, 3:], tetra + 1)
            boundary = np.loadtxt(path / "mesh.boundary", dtype=int)
            np.testing.assert_array_equal(boundary[:, 1], labels)
            np.testing.assert_array_equal(boundary[:, 2], parents + 1)
            np.testing.assert_array_equal(boundary[:, 5:], faces + 1)
            self.assertEqual((path / "mesh.header").read_text(), "5 2 6\n2\n504 2\n303 6\n")

    def test_all_edges_and_unreachable_target(self):
        faces = np.array([[0, 1, 2], [1, 0, 3], [2, 1, 4], [0, 2, 5]])
        neighbors = triangle_neighbors(faces)
        self.assertEqual(set(neighbors[0]), {1, 2, 3})
        patch, info = grow_patch(faces, neighbors, np.ones(4), np.ones(4, dtype=bool), 3.5, random.Random(2))
        self.assertEqual(set(patch), {0, 1, 2, 3})
        self.assertEqual(info["area_m2"], 4)
        with self.assertRaises(PatchSizeError):
            grow_patch(faces, neighbors, np.ones(4), np.array([False, True, True, True]), 2, random.Random(2))
        with self.assertRaises(PatchSizeError):
            grow_patch(faces, neighbors, np.ones(4), np.zeros(4, dtype=bool), 1, random.Random(2))

    def test_margin_covers_whole_triangles(self):
        # The last triangle's centroid is far away, but one corner is within 2 cm.
        points = np.array([[-1., -1., 0.], [1., -1., 0.], [0., 1., 0.],
                           [0., 0., .019], [.001, 0., .019], [0., .001, .019],
                           [0., 0., .025], [.001, 0., .025], [0., .001, .025],
                           [0., 0., .01], [.1, 0., .1], [0., .1, .1]])
        faces = np.arange(12).reshape(-1, 3)
        allowed = outside_fixed_margin(points, faces, np.array([0]))
        np.testing.assert_array_equal(allowed, [False, False, True, False])

    def test_sampling_disjoint_connected_reproducible_and_force_scaling(self):
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(.2)
        sphere.SetThetaResolution(32)
        sphere.SetPhiResolution(32)
        sphere.Update()
        mesh = sphere.GetOutput()
        points = vtk_to_numpy(mesh.GetPoints().GetData()).astype(float)
        faces = vtk_to_numpy(mesh.GetPolys().GetConnectivityArray()).reshape(-1, 3)
        neighbors = triangle_neighbors(faces)
        tri = points[faces]
        areas = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
        for seed in range(5):
            labels, fixed, loads = sample_boundaries(points, faces, random.Random(seed), 1.5)
            doubled, same_fixed, bigger = sample_boundaries(points, faces, random.Random(seed), 3.)
            np.testing.assert_array_equal(labels, doubled)
            self.assertEqual(fixed, same_fixed)
            smaller, small_fixed, _ = sample_boundaries(
                points, faces, random.Random(seed), 1.5, fixed_area_scale=1.)
            self.assertAlmostEqual(small_fixed["target_area_m2"], fixed["target_area_m2"] / 3)
            self.assertEqual(small_fixed["seed_vertex_id"], fixed["seed_vertex_id"])
            self.assertEqual(small_fixed["area_scale"], 1.)
            self.assertTrue((labels[smaller == 2] == 2).all())
            small_eligible = outside_fixed_margin(points, faces, np.flatnonzero(smaller == 2))
            self.assertTrue(small_eligible[smaller >= 3].all())
            with mock_patch("boundaries.outside_fixed_margin") as clearance:
                touching, no_margin_fixed, touching_loads = sample_boundaries(
                    points, faces, random.Random(seed), 1.5,
                    fixed_area_scale=1., fixed_margin=0.)
            clearance.assert_not_called()
            self.assertEqual(no_margin_fixed["margin_m"], 0.)
            np.testing.assert_array_equal(touching == 2, smaller == 2)
            for load in touching_loads:
                self.assertAlmostEqual(areas[touching == load["boundary_id"]].sum(), load["area_m2"])
            eligible = outside_fixed_margin(points, faces, np.flatnonzero(labels == 2))
            self.assertTrue(eligible[labels >= 3].all())
            for info in [fixed, *loads]:
                patch = set(np.flatnonzero(labels == info["boundary_id"]))
                self.assertIn(info["seed_vertex_id"] - 1, faces[list(patch)])
                self.assertAlmostEqual(areas[list(patch)].sum(), info["area_m2"])
                self.assertGreaterEqual(info["area_m2"], info["target_area_m2"])
                self.assertLess(info["area_m2"] - info["target_area_m2"], areas.max())
                queue = [patch.pop()]
                while queue:
                    for other in neighbors[queue.pop()]:
                        if other in patch:
                            patch.remove(other)
                            queue.append(other)
                self.assertFalse(patch)
            for first, second in zip(loads, bigger):
                np.testing.assert_allclose(second["traction_pa"], np.array(first["traction_pa"]) * 2, atol=1.1e-6, rtol=0)
                self.assertAlmostEqual(np.linalg.norm(first["traction_pa"]) * first["area_m2"], first["force_magnitude_n"], places=6)


if __name__ == "__main__":
    unittest.main()
