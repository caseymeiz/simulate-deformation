"""Show liver boundary patches beside displacement arrows using ParaView."""

import argparse
import json
from math import ceil
from pathlib import Path

import numpy as np
from paraview import servermanager
from vtk.util.numpy_support import vtk_to_numpy

from paraview.simple import (
    AssignViewToLayout, Calculator, ColorBy, CreateLayout, CreateView, ExtractSurface, GetActiveCamera,
    GetColorTransferFunction, GetScalarBar, Glyph, OpenDataFile, PointSource, Render,
    RenameSource, ResetCamera, SaveScreenshot, SaveState, SetActiveView, Show, Text, Threshold, WarpByVector,
)

from boundary_to_vtk import boundary_to_vtk


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path, help="deformed.vtk from a successful simulation")
    parser.add_argument("--boundary", type=Path,
                        help="Elmer mesh.boundary; default: beside the deformed mesh")
    parser.add_argument("--output", type=Path, help="PNG path; default: displacement.png beside the mesh")
    parser.add_argument("--glyph-scale", type=float, default=1.,
                        help="Arrow length multiplier; 1 shows the actual displacement")
    parser.add_argument("--force-scale", type=float, default=.02,
                        help="Force arrow length in metres per newton (default: 0.02)")
    parser.add_argument("--max-glyphs", type=int, default=1000,
                        help="Approximate maximum number of surface arrows")
    args = parser.parse_args()
    output = (args.output or args.mesh.with_name("displacement.png")).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    deformed = OpenDataFile(str(args.mesh.resolve()))
    RenameSource("Deformed liver", deformed)
    # Elmer's exported coordinates already include the displacement.
    reference = WarpByVector(registrationName="Original liver", Input=deformed)
    reference.Vectors = ["POINTS", "displacement"]
    reference.ScaleFactor = -1.
    surface = ExtractSurface(registrationName="Original surface", Input=reference)
    magnitude = Calculator(registrationName="Displacement magnitude (mm)", Input=surface)
    magnitude.ResultArrayName = "displacement_mm"
    magnitude.Function = "1000 * mag(displacement)"
    magnitude.UpdatePipeline()

    layout = CreateLayout(name="Liver deformation")
    patch_view = CreateView("RenderView")
    view = CreateView("RenderView")
    AssignViewToLayout(patch_view, layout, hint=0)
    layout.SplitHorizontal(0, .5)
    AssignViewToLayout(view, layout, hint=2)
    layout.SetSize(2400, 1000)
    for panel in (patch_view, view):
        panel.UseColorPaletteForBackground = 0
        panel.BackgroundColorMode = "Single Color"
        panel.Background = [1., 1., 1.]
        panel.CameraParallelProjection = 1
        panel.DepthPeeling = 1

    boundary_path = args.boundary or args.mesh.with_name("mesh.boundary")
    loads = json.loads(boundary_path.with_name("manifest.json").read_text())["tractions"]
    boundary_surface = output.with_name(output.stem + "_boundaries.vtk")
    boundary_to_vtk(boundary_path, boundary_surface)
    boundaries = OpenDataFile(str(boundary_surface))
    RenameSource("Original boundary labels", boundaries)
    boundaries.UpdatePipeline()
    last_label = int(boundaries.CellData["boundary_id"].GetRange()[1])
    free = Threshold(registrationName="Free surface", Input=boundaries)
    free.Scalars = ["CELLS", "boundary_id"]
    free.LowerThreshold = 1
    free.UpperThreshold = 1
    free_display = Show(free, patch_view)
    free_display.ColorArrayName = ["POINTS", ""]
    free_display.Representation = "Surface With Edges"
    free_display.DiffuseColor = [.78, .8, .83]
    free_display.EdgeColor = [.35, .38, .42]
    free_display.Opacity = .2

    fixed = Threshold(registrationName="Fixed patch", Input=boundaries)
    fixed.Scalars = ["CELLS", "boundary_id"]
    fixed.LowerThreshold = 2
    fixed.UpperThreshold = 2
    fixed_display = Show(fixed, patch_view)
    fixed_display.ColorArrayName = ["POINTS", ""]
    fixed_display.DiffuseColor = [0., 0., 0.]

    patches = Threshold(registrationName="Force patches", Input=boundaries)
    patches.Scalars = ["CELLS", "boundary_id"]
    patches.LowerThreshold = 3
    patches.UpperThreshold = last_label
    patch_display = Show(patches, patch_view)
    ColorBy(patch_display, ("CELLS", "boundary_id"))
    patch_colors = GetColorTransferFunction("boundary_id")
    patch_colors.InterpretValuesAsCategories = 1
    patch_colors.Annotations = ["2", "Fixed"] + [
        value for load in loads for value in (
            str(load["boundary_id"]),
            f"Force {load['boundary_id'] - 2} ({load['force_magnitude_n']:.2f} N)")]
    patch_display.Opacity = .45
    palette = [(0., 0., 0.), (.9, .32, .12), (0., .45, .7), (.65, .25, .65)]
    patch_colors.IndexedColors = [value for color in palette[:last_label - 1] for value in color]
    patch_display.SetScalarBarVisibility(patch_view, True)
    patch_legend = GetScalarBar(patch_colors, patch_view)
    patch_legend.Title = "Boundary patches"
    patch_legend.ComponentTitle = ""
    patch_legend.TitleColor = [.15, .15, .15]
    patch_legend.LabelColor = [.15, .15, .15]
    patch_legend.DrawNanAnnotation = 0
    patch_caption = Text()
    patch_caption.Text = ("Original mesh + applied forces\n"
                          f"One resultant per patch; arrow scale: {args.force_scale * 1000:g} mm/N")
    patch_caption_display = Show(patch_caption, patch_view)
    patch_caption_display.WindowLocation = "Upper Left Corner"
    patch_caption_display.Color = [.15, .15, .15]
    patch_caption_display.FontSize = 18

    # Integrate the saved uniform traction over each patch to obtain its net force.
    boundary_data = servermanager.Fetch(boundaries)
    points = vtk_to_numpy(boundary_data.GetPoints().GetData())
    faces = vtk_to_numpy(boundary_data.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    labels = vtk_to_numpy(boundary_data.GetCellData().GetArray("boundary_id"))
    for load in loads:
        label = load["boundary_id"]
        triangles = points[faces[labels == label]]
        centers = triangles.mean(axis=1)
        areas = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                        triangles[:, 2] - triangles[:, 0]), axis=1) / 2
        center = np.average(centers, axis=0, weights=areas)
        anchor = centers[np.linalg.norm(centers - center, axis=1).argmin()]
        location = PointSource(registrationName=f"Force {label - 2} anchor")
        location.Center = anchor.tolist()
        location.NumberOfPoints = 1
        location.Radius = 0.
        force = np.array(load["traction_pa"]) * load["area_m2"]
        vectors = Calculator(registrationName=f"Force {label - 2} (N)", Input=location)
        vectors.ResultArrayName = "force_n"
        vectors.Function = " + ".join(f"({value:.17g})*{axis}" for value, axis in zip(force, ("iHat", "jHat", "kHat")))
        force_glyph = Glyph(registrationName=f"Force {label - 2} arrow", Input=vectors, GlyphType="Arrow")
        force_glyph.OrientationArray = ["POINTS", "force_n"]
        force_glyph.ScaleArray = ["POINTS", "force_n"]
        force_glyph.ScaleFactor = args.force_scale
        force_glyph.GlyphMode = "All Points"
        force_glyph.GlyphType.ShaftRadius = .04
        force_glyph.GlyphType.TipRadius = .12
        force_display = Show(force_glyph, patch_view)
        force_display.ColorArrayName = ["POINTS", ""]
        force_display.DiffuseColor = palette[label - 2]
        force_display.AmbientColor = palette[label - 2]
        force_display.Ambient = .5

    original_display = Show(surface, view)
    original_display.ColorArrayName = ["POINTS", ""]
    original_display.DiffuseColor = [0.65, 0.68, 0.72]
    original_display.Opacity = .3

    glyphs = Glyph(registrationName="Surface displacement arrows", Input=magnitude, GlyphType="Arrow")
    glyphs.OrientationArray = ["POINTS", "displacement"]
    glyphs.ScaleArray = ["POINTS", "displacement"]
    glyphs.ScaleFactor = args.glyph_scale
    glyphs.GlyphMode = "Every Nth Point"
    glyphs.Stride = max(1, ceil(magnitude.GetDataInformation().GetNumberOfPoints() / args.max_glyphs))
    arrows = Show(glyphs, view)
    ColorBy(arrows, ("POINTS", "displacement_mm"))
    colors = GetColorTransferFunction("displacement_mm")
    colors.ApplyPreset("Viridis (matplotlib)", True)
    colors.RescaleTransferFunction(*magnitude.PointData["displacement_mm"].GetRange())
    arrows.SetScalarBarVisibility(view, True)
    legend = GetScalarBar(colors, view)
    legend.Title = "Displacement (mm)"
    legend.ComponentTitle = ""
    legend.TitleColor = [.15, .15, .15]
    legend.LabelColor = [.15, .15, .15]

    caption = Text()
    caption.Text = f"Original liver + displacement\nArrow length multiplier: {args.glyph_scale:g}"
    caption_display = Show(caption, view)
    caption_display.WindowLocation = "Upper Left Corner"
    caption_display.Color = [.15, .15, .15]
    caption_display.FontSize = 18
    SetActiveView(view)
    ResetCamera(view)
    camera = GetActiveCamera()
    camera.Azimuth(30)
    camera.Elevation(20)
    ResetCamera(view)
    for name in ("CameraPosition", "CameraFocalPoint", "CameraViewUp", "CameraParallelScale"):
        setattr(patch_view, name, getattr(view, name))
    Render(patch_view)
    Render(view)
    SaveScreenshot(str(output), layout, ImageResolution=[2400, 1000])
    SaveState(str(output.with_suffix(".pvsm")))
    print(f"Saved {output} and {output.with_suffix('.pvsm')}")


if __name__ == "__main__":
    main()
