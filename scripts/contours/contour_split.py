

import volmdlr as vm
import volmdlr.edges as edges
import volmdlr.wires as wires
import volmdlr.primitives2d as p2d

u = vm.Vector2D.random(0, 1, 0, 1)
u.normalize()
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

ax = contour.plot()
line = edges.Line2D(vm.O2D, u)
# 
ax2 = line.plot(ax=ax, color='b')

split_contours = contour.cut_by_line(line)
for c in split_contours:
    c.plot(color='r')

