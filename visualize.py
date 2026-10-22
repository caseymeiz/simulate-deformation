"""Show liver boundary patches beside displacement arrows using ParaView."""

import argparse
from math import ceil
from pathlib import Path

from paraview.simple import (
    AssignViewToLayout, Calculator, ColorBy, CreateLayout, CreateView, ExtractSurface, GetActiveCamera,
    GetColorTransferFunction, GetScalarBar, Glyph, OpenDataFile, Render,
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

    patches = Threshold(registrationName="Fixed and force patches", Input=boundaries)
    patches.Scalars = ["CELLS", "boundary_id"]
    patches.LowerThreshold = 2
    patches.UpperThreshold = last_label
    patch_display = Show(patches, patch_view)
    ColorBy(patch_display, ("CELLS", "boundary_id"))
    patch_colors = GetColorTransferFunction("boundary_id")
    patch_colors.InterpretValuesAsCategories = 1
    patch_colors.Annotations = [value for label in range(2, last_label + 1)
                                for value in (str(label), "Fixed" if label == 2 else f"Force {label - 2}")]
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
    patch_caption.Text = "Original mesh + boundary patches"
    patch_caption_display = Show(patch_caption, patch_view)
    patch_caption_display.WindowLocation = "Upper Left Corner"
    patch_caption_display.Color = [.15, .15, .15]
    patch_caption_display.FontSize = 18

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
