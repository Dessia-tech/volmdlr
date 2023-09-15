"""
Some conical surfaces models to perform unit test.
"""
import math

import volmdlr
from volmdlr import surfaces
conical_surface1 = surfaces.ConicalSurface3D(volmdlr.OXYZ, math.pi / 3)

frame_cone = volmdlr.Frame3D(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.X3D, -volmdlr.Y3D, -volmdlr.Z3D)
conical_surface2 = surfaces.ConicalSurface3D(frame_cone, 0.336674819387)

frame = volmdlr.Frame3D(volmdlr.Point3D(0.0, 0.0, -0.0022143719324716786), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
conical_surface3 = surfaces.ConicalSurface3D(frame, 1.0471975511966)
