"""Open rounded line segments model."""
import random

import volmdlr
from volmdlr import primitives3d

random.seed(2)

p1 = volmdlr.Point3D(0, 0, 0)
p2 = volmdlr.Point3D(-0.150, 0, 0)
p3 = volmdlr.Point3D(-0.150, 0.215, 0)
p4 = volmdlr.Point3D(-0.150, 0.215, -0.058)
p5 = volmdlr.Point3D(-0.220, 0.186, -0.042)

points = [p1, p2, p3, p4, p5]
radius = {1: 0.015, 2: 0.020, 3: 0.03}

current_point = p5

for i in range(6):
    current_point += volmdlr.Point3D.random(-0.1, 0.3, -0.1, 0.3, -0.1, 0.3)
    points.append(current_point)
    radius[4 + i] = 0.01 + 0.03 * random.random()


open_rounded_line_segements = primitives3d.OpenRoundedLineSegments3D(points, radius, adapt_radius=True, name='wire')
