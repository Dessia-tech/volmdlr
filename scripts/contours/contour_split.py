

import volmdlr as vm
import volmdlr.edges as edges
import volmdlr.wires as wires
import volmdlr.primitives2d as p2d
import matplotlib.pyplot as plt

plt, (ax1, ax2, ax3) = plt.subplots(1, 3)

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

ax = contour.plot(ax=ax1)
line = edges.Line2D(vm.O2D, u)
line2 = edges.Line2D(0.3*v, u+0.3*v)

line.plot(ax=ax1, color='b')
line2.plot(ax=ax1, color='g')

split_contours1 = contour.cut_by_line(line)
for c in split_contours1:
    c.plot(ax=ax2, color='b')
ax2.set_title('{} splitted contours'.format(len(split_contours1)))

split_contours2 = contour.cut_by_line(line2)
for c in split_contours2:
    c.plot(ax=ax3, color='g')
ax3.set_title('{} splitted contours'.format(len(split_contours2)))