# -*- coding: utf-8 -*-

import math
import pkg_resources
__version__ = pkg_resources.require("volmdlr")[0].version

from volmdlr.core_compiled import (Vector2D, Vector3D, Point2D, Point3D,
                            O2D, X2D, Y2D, OXY,
                            Basis2D, Basis3D, Frame2D, Frame3D,
                            O3D, X3D, Y3D, Z3D, Matrix22, Matrix33
                            )

TWO_PI = 2*math.pi

XYZ = Basis3D(X3D, Y3D, Z3D)
YZX = Basis3D(Y3D, Z3D, X3D)
ZXY = Basis3D(Z3D, X3D, Y3D)

OXYZ = Frame3D(O3D, X3D, Y3D, Z3D)
OYZX = Frame3D(O3D, Y3D, Z3D, X3D)
OZXY = Frame3D(O3D, Z3D, X3D, Y3D)