"""Show the original liver surface with displacement arrows using ParaView."""

import argparse
from math import ceil
from pathlib import Path

from paraview.simple import (
    Calculator, ColorBy, ExtractSurface, GetActiveCamera, GetActiveViewOrCreate,
    GetColorTransferFunction, GetScalarBar, Glyph, OpenDataFile, Render,
    RenameSource, ResetCamera, SaveScreenshot, SaveState, Show, Text, WarpByVector,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path, help="deformed.vtk from a successful simulation")
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

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [1400, 1000]
    view.UseColorPaletteForBackground = 0
    view.BackgroundColorMode = "Single Color"
    view.Background = [1., 1., 1.]
    view.CameraParallelProjection = 1
    view.DepthPeeling = 1
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
    ResetCamera(view)
    camera = GetActiveCamera()
    camera.Azimuth(30)
    camera.Elevation(20)
    ResetCamera(view)
    Render(view)
    SaveScreenshot(str(output), view, ImageResolution=view.ViewSize)
    SaveState(str(output.with_suffix(".pvsm")))
    print(f"Saved {output} and {output.with_suffix('.pvsm')}")


if __name__ == "__main__":
    main()
