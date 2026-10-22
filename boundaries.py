"""Grow disjoint surface patches. Coordinates are metres, areas are square metres."""

# Force and target-area distributions adapted from Volume2SurfaceCNN.
# Portions copyright (c) 2020 National Center for Tumor Diseases,
# Division of Translational Surgical Oncology. See VOLUME2SURFACECNN_LICENSE.txt.

from collections import deque

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

FIXED_MARGIN = 0.02
FIXED_AREA_SCALE = 3.0
FORCE_RADIUS_MAX = 0.07


class PatchSizeError(ValueError):
    """The sampled connected patch cannot reach its target area."""


def triangle_neighbors(faces):
    neighbors = [[] for _ in faces]
    edges = {}
    for index, (a, b, c) in enumerate(faces):
        for edge in ((a, b), (b, c), (c, a)):
            key = tuple(sorted(edge))
            if key in edges:
                other = edges.pop(key)
                neighbors[index].append(other)
                neighbors[other].append(index)
            else:
                edges[key] = index
    return neighbors


def grow_patch(faces, neighbors, areas, allowed, target_area, rng):
    """Choose a surface vertex, then grow through eligible edge neighbors."""
    vertices = np.unique(faces[allowed])
    if not len(vertices):
        raise PatchSizeError("No eligible surface remains for a patch.")
    seed = int(rng.choice(vertices))
    incident = np.flatnonzero(allowed & np.any(faces == seed, axis=1))
    first = int(rng.choice(incident))
    queue, seen, patch = deque([first]), {first}, []
    area = 0.0
    while queue:
        face = queue.popleft()
        patch.append(face)
        area += float(areas[face])
        if area >= target_area:
            return np.array(patch), {"seed_vertex_id": seed + 1,
                                     "target_area_m2": target_area, "area_m2": area}
        for other in neighbors[face]:
            if allowed[other] and other not in seen:
                seen.add(other)
                queue.append(other)
    raise PatchSizeError(f"Patch at vertex {seed + 1} needs {target_area:.6g} m²; "
                         f"its eligible connected region has only {area:.6g} m².")


def outside_fixed_margin(points, faces, fixed, margin=FIXED_MARGIN):
    """Conservatively keep whole triangles at least margin from the fixed patch.

    Distance from a triangle centroid to the fixed surface, minus its bounding
    radius, is a lower bound on the distance of every point in that triangle.
    VTK is used only for this in-memory closest-surface query, not mesh export.
    """
    surface = vtk.vtkPolyData()
    vertices = vtk.vtkPoints()
    vertices.SetData(numpy_to_vtk(points, deep=True))
    surface.SetPoints(vertices)
    cells = vtk.vtkCellArray()
    cells.ImportLegacyFormat(numpy_to_vtkIdTypeArray(
        np.column_stack((np.full(len(fixed), 3), faces[fixed])).astype(np.int64).ravel(), deep=True))
    surface.SetPolys(cells)
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()
    triangles = points[faces]
    centers = triangles.mean(axis=1)
    radii = np.linalg.norm(triangles - centers[:, None], axis=2).max(axis=1)
    closest, cell, sub, squared = [0.0] * 3, vtk.reference(0), vtk.reference(0), vtk.reference(0.0)
    distances = np.empty(len(faces))
    for i, center in enumerate(centers):
        locator.FindClosestPoint(center, closest, cell, sub, squared)
        distances[i] = float(squared) ** 0.5
    return distances - radii >= margin


def sample_boundaries(points, faces, rng, force_factor, min_force_patches=1,
                      fixed_area_scale=FIXED_AREA_SCALE, fixed_margin=FIXED_MARGIN,
                      force_radius_max=FORCE_RADIUS_MAX):
    triangles = points[faces]
    areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                    triangles[:, 2] - triangles[:, 0]), axis=1) / 2
    neighbors = triangle_neighbors(faces)
    labels = np.ones(len(faces), dtype=int)
    radius = 0.025 + 0.05 * rng.random() * np.ptp(points, axis=0).max()
    fixed, fixed_info = grow_patch(faces, neighbors, areas, labels == 1,
                                   float(fixed_area_scale * np.pi * radius ** 2), rng)
    labels[fixed] = 2
    fixed_info.update(boundary_id=2, margin_m=fixed_margin, area_scale=fixed_area_scale)
    allowed = (outside_fixed_margin(points, faces, fixed, fixed_margin)
               if fixed_margin else np.ones(len(faces), dtype=bool))
    allowed[fixed] = False
    loads = []
    for boundary_id in range(3, 3 + rng.randint(min_force_patches, 3)):
        radius = 0.01 + (force_radius_max - 0.01) * rng.random()
        patch, info = grow_patch(faces, neighbors, areas, allowed, np.pi * radius ** 2, rng)
        labels[patch] = boundary_id
        allowed[patch] = False
        direction = np.array([rng.random() - 0.5 for _ in range(3)])
        direction /= np.linalg.norm(direction)
        magnitude = force_factor * (1 - rng.random() ** 2)
        traction = np.round(direction * magnitude / info["area_m2"], 6)
        info.update(boundary_id=boundary_id, sampled_radius_m=radius, force_magnitude_n=magnitude,
                    traction_pa=traction.tolist())
        loads.append(info)
    return labels, fixed_info, loads
