"""Export an Elmer mesh.boundary and its sibling mesh.nodes as a labeled VTK surface."""

import argparse
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray


def boundary_to_vtk(boundary_path, output_path):
    """Write triangle group IDs and vertex masks; coordinates remain in metres."""
    boundary_path = Path(boundary_path)
    nodes = np.loadtxt(boundary_path.with_suffix(".nodes"), ndmin=2)
    boundary = np.loadtxt(boundary_path, dtype=np.int64, ndmin=2)
    node_ids, connectivity = np.unique(boundary[:, 5:8], return_inverse=True)
    triangles = connectivity.reshape(-1, 3)
    nodes = nodes[np.argsort(nodes[:, 0])]
    coordinates = nodes[np.searchsorted(nodes[:, 0], node_ids), 2:5]
    surface = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(coordinates, deep=True))
    surface.SetPoints(points)
    cells = vtk.vtkCellArray()
    cells.ImportLegacyFormat(numpy_to_vtkIdTypeArray(
        np.column_stack((np.full(len(triangles), 3), triangles)).astype(np.int64).ravel(), deep=True))
    surface.SetPolys(cells)
    for name, values in (("boundary_id", boundary[:, 1]), ("triangle_id", boundary[:, 0]),
                         ("parent_tetra_id", boundary[:, 2])):
        array = numpy_to_vtk(values, deep=True)
        array.SetName(name)
        surface.GetCellData().AddArray(array)
    surface.GetCellData().SetActiveScalars("boundary_id")
    fixed, force = np.zeros(len(node_ids), dtype=np.uint8), np.zeros(len(node_ids), dtype=np.uint8)
    fixed[triangles[boundary[:, 1] == 2].ravel()] = 1
    force[triangles[boundary[:, 1] >= 3].ravel()] = 1
    for name, values in (("elmer_node_id", node_ids), ("is_fixed", fixed), ("is_force", force)):
        array = numpy_to_vtk(values, deep=True)
        array.SetName(name)
        surface.GetPointData().AddArray(array)
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(str(output_path))
    writer.SetInputData(surface)
    writer.SetFileTypeToBinary()
    writer.Write()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("boundary", type=Path, help="mesh.boundary; mesh.nodes must be alongside it")
    parser.add_argument("output", type=Path, help="Output .vtk surface")
    args = parser.parse_args()
    boundary_to_vtk(args.boundary, args.output)
