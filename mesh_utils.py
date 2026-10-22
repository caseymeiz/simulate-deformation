"""Mesh and result helpers copied from our MSD deformation scripts."""

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def boundary_faces(tetra):
    faces = tetra[:, ((0, 1, 2), (0, 3, 1), (0, 2, 3), (1, 3, 2))].reshape(-1, 3)
    _, inverse, counts = np.unique(np.sort(faces, axis=1), axis=0, return_inverse=True, return_counts=True)
    surface = np.flatnonzero(counts[inverse] == 1)
    return faces[surface], surface // 4


def write_legacy(path, mesh):
    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(mesh)
    writer.SetFileTypeToBinary()
    writer.Write()


def check_result_geometry(result_path, tetra):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(result_path))
    reader.Update()
    mesh = reader.GetOutput()
    points = vtk_to_numpy(mesh.GetPoints().GetData())
    displacement = vtk_to_numpy(mesh.GetPointData().GetArray("displacement"))
    finite = bool(np.isfinite(points).all() and np.isfinite(displacement).all())
    reference = (points - displacement)[tetra]
    deformed = points[tetra]
    original_volume = np.linalg.det(reference[:, 1:] - reference[:, :1])
    deformed_volume = np.linalg.det(deformed[:, 1:] - deformed[:, :1])
    inverted = int(np.count_nonzero(deformed_volume * original_volume <= 0))
    return {"finite_geometry": finite, "inverted_tetrahedra": inverted}
