"""
Showcase of BoundingBox / Triangle3D intersection
"""
from volmdlr.core import BoundingBox
from volmdlr.faces import Triangle3D
import volmdlr

bounding_box1 = BoundingBox(0, 1, 0, 1, 0, 1)
bounding_box2 = BoundingBox(1, 2, 1, 2, 1, 2)

# p0 = volmdlr.Point3D(5, -1, 0.5)
# p1 = volmdlr.Point3D(-1, 5, 0.5)
# p2 = volmdlr.Point3D(5, 5, 0.5)
# triangle = Triangle3D(p0, p1, p2)
#
# ax = bounding_box.plot()
# triangle.plot(ax=ax)
#
# print(bounding_box.triangle_intersection(triangle))
