"""
Showcase of BoundingBox / Triangle3D intersection
"""
from volmdlr.core import BoundingBox
from volmdlr.faces import Triangle3D
import volmdlr

bounding_box = BoundingBox(0.0, 2.0, 0.0, 2.0, 0.0, 2.0)

p0 = volmdlr.Point3D(5, -3, 0.5)
p1 = volmdlr.Point3D(-3, 5, 0.5)
p2 = volmdlr.Point3D(5, 5, 0.5)
triangle = Triangle3D(p0, p1, p2)

ax = bounding_box.plot()
triangle.plot(ax=ax)

print(bounding_box.is_intersecting_triangle(triangle))
Uni