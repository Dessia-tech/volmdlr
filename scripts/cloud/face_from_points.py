# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.cloud as vmc
import volmdlr.faces as faces

p1 = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)
p2 = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)
p3 = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)


tri = faces.Triangle3D(p1, p2, p3)
# tri.babylonjs()
# points_sub = tri.subdescription()
triangles3d = tri.subdescription_to_triangles(1e-2)

ax = tri.plot()
# # for pt in points_sub:
# #     pt.plot(ax=ax)
for tri in triangles3d:
    tri.color = (250, 0, 250)
    tri.alpha = 0.5
    tri.plot(ax=ax, color='r')

vol = vm.core.VolumeModel([tri] + triangles3d)
vol.babylonjs()
