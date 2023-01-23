"""
Some conical surfaces models to perfom unittest.
"""
import math

import volmdlr
import volmdlr.faces as vmf

conical_surface1 = vmf.ConicalSurface3D(volmdlr.OXYZ, math.pi / 3)

frame_cone = volmdlr.Frame3D(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.X3D, -volmdlr.Y3D, -volmdlr.Z3D)
conical_surface2 = vmf.ConicalSurface3D(frame_cone, 0.336674819387)

frame = volmdlr.Frame3D(volmdlr.Point3D(0.0, 0.0, -0.0022143719324716786), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
conical_surface3 = vmf.ConicalSurface3D(frame, 1.0471975511966)
