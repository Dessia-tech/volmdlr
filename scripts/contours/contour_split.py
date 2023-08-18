

import matplotlib.pyplot as plt

import volmdlr as vm
from volmdlr.core import EdgeStyle
import volmdlr.primitives2d as p2d
from volmdlr import edges, curves
from volmdlr.utils.common_operations import random_color
plt, (ax1, ax2, ax3) = plt.subplots(1, 3)

u = vm.Vector2D.random(0, 1, 0, 1)
u = u.unit_vector()
v = u.normal_vector()

l = 0.05

p1 = v * l
p2 = p1 + l*u
p3 = p2 - 2*l*v
p4 = p3 + l*u
p5 = p4 + 2*l*u + 3*l*v
p6 = p5 + l*u
p7 = p6 - 4*l*v + l*v
p8 = p7 - l*u
p9 = 0.5*(p5 + p6) - l*v
p10 = p4 - l*v
p11 = p1 - 3*l*v

contour =p2d.ClosedRoundedLineSegments2D(
    [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11],
    {3:0.3*l})

ax = contour.plot(ax=ax1)
line = curves.Line2D(vm.O2D, u)
line2 = curves.Line2D(0.3*v, u+0.3*v)

line.plot(ax=ax1, edge_style=EdgeStyle(color='b'))
line2.plot(ax=ax1, edge_style=EdgeStyle(color='g'))

split_contours1 = contour.cut_by_line(line)
for c in split_contours1:
    c.plot(ax=ax2, edge_style=EdgeStyle(color=random_color()))
ax2.set_title('{} splitted contours'.format(len(split_contours1)))

split_contours2 = contour.cut_by_line(line2)
for c in split_contours2:
    c.plot(ax=ax3, edge_style=EdgeStyle(color='g'))
ax3.set_title('{} splitted contours'.format(len(split_contours2)))
