# -*- coding: utf-8 -*-

import math
import pkg_resources
__version__ = pkg_resources.require("volmdlr")[0].version

from volmdlr.core_compiled import (Vector2D, Vector3D, Point2D, Point3D,
                            O2D, X2D, Y2D, OXY,
                            Basis2D, Basis3D, Frame2D, Frame3D,
                            O3D, X3D, Y3D, Z3D,
                            )

two_pi = 2*math.pi