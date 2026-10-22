"""Write Elmer's Tet4 mesh directly from NumPy coordinates and connectivity."""

import numpy as np


def write_volume(directory, points, tetra, surface_count):
    nodes = np.column_stack((np.arange(1, len(points) + 1), np.full(len(points), -1), points))
    np.savetxt(directory / "mesh.nodes", nodes, fmt=["%d", "%d", "%.17g", "%.17g", "%.17g"])
    elements = np.column_stack((np.arange(1, len(tetra) + 1), np.ones(len(tetra), dtype=int),
                                np.full(len(tetra), 504), tetra + 1))
    np.savetxt(directory / "mesh.elements", elements, fmt="%d")
    (directory / "mesh.header").write_text(
        f"{len(points)} {len(tetra)} {surface_count}\n2\n504 {len(tetra)}\n303 {surface_count}\n")


def write_boundary(directory, faces, parents, labels):
    count = len(faces)
    rows = np.column_stack((np.arange(1, count + 1), labels, parents + 1,
                            np.zeros(count, dtype=int), np.full(count, 303), faces + 1))
    np.savetxt(directory / "mesh.boundary", rows, fmt="%d")
