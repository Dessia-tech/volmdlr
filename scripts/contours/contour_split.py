

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

line_seg1 = edges.LineSegment2D(vm.Point2D(-0.5, -0.2), vm.O2D)
line_seg2 = edges.LineSegment2D(vm.O2D, vm.Point2D(0.3, 1))
line_seg3 = edges.LineSegment2D(vm.Point2D(0.3, 1), vm.Point2D(1, 1))
line_seg4 = edges.LineSegment2D(vm.Point2D(1, 1), vm.Point2D(1, -0.5))
line_seg5 = edges.LineSegment2D(vm.Point2D(1, -0.5), vm.Point2D(-0.5, -0.2))

contour = wires.Contour2D([line_seg1, line_seg2, line_seg3, line_seg4, line_seg5])
ax = contour.plot()

line1 = edges.Line2D(vm.Point2D(-0.5, 1), vm.O2D)
line1.plot(ax=ax, color='r')

cute_wire_line1 = contour.cut_by_line(line1)
for c in cute_wire_line1:
    c.plot()
assert len(cute_wire_line1) == 2

line2 = edges.Line2D(vm.Point2D(-0.5, -0.5), vm.O2D)
line2.plot(ax=ax, color='b')

cute_wire_line2 = contour.cut_by_line(line2)
for c in cute_wire_line2:
    c.plot()
assert len(cute_wire_line2) == 2
