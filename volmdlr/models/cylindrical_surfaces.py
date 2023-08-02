"""
Some cylindrical surfaces models to perform unit test.
"""

import volmdlr
from volmdlr import surfaces

cylindrical_surface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 0.32)
cylindrical_surface2 = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 1.0)
frame = volmdlr.Frame3D(volmdlr.Point3D(-0.005829, 0.000765110438227, -0.0002349369830163),
                        volmdlr.Vector3D(-0.6607898454031987, 0.562158151695499, -0.4973278523210991),
                        volmdlr.Vector3D(-0.7505709694705869, -0.4949144228333324, 0.43783893597935386),
                        volmdlr.Vector3D(-0.0, 0.6625993710787045, 0.748974013865705))
cylindrical_surface3 = surfaces.CylindricalSurface3D(frame, 0.003)
cylindrical_surface4 = surfaces.CylindricalSurface3D(volmdlr.OXYZ, radius=0.03)
