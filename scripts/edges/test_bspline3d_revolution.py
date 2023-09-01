
import volmdlr
from volmdlr.models.edges import bspline_curve3d
bspline_curve_3d = bspline_curve3d()

face = bspline_curve_3d.revolution(volmdlr.O3D, volmdlr.X3D, volmdlr.TWO_PI)
face.babylonjs()
