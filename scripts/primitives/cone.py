 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm

# from volmdlr.primitives2d import ClosedRoundedLineSegments2D, OpenedRoundedLineSegments2D
from volmdlr.primitives3d import Cone



cone = Cone(position=vm.O3D, axis=vm.Z3D, radius = 0.02, length=0.05)
cone.babylonjs()

cone2 = cone.translation(vm.X3D)
cone3 = cone.rotation(-vm.Y3D, vm.X3D, 0.3)

cone3.babylonjs()