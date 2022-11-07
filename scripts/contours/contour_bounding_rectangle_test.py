import matplotlib.pyplot as plt

import volmdlr as vm
from volmdlr.core import BoundingRectangle
import volmdlr.wires as vmw
import volmdlr.edges as vme

line_seg1 = vme.LineSegment2D(vm.Point2D(-0.5, -0.2), vm.O2D)
line_seg2 = vme.LineSegment2D(vm.O2D, vm.Point2D(0.3, 1))
line_seg3 = vme.LineSegment2D(vm.Point2D(0.3, 1), vm.Point2D(1, 1))
line_seg4 = vme.LineSegment2D(vm.Point2D(1, 1), vm.Point2D(1, -0.5))
line_seg5 = vme.LineSegment2D(vm.Point2D(1, -0.5), vm.Point2D(-0.5, -0.2))

contour1 = vmw.Contour2D([line_seg1, line_seg2, line_seg3, line_seg4, line_seg5])

bd_rectangle = contour1.bounding_rectangle()

xmin = bd_rectangle[0]
xmax = bd_rectangle[1]
ymin = bd_rectangle[2]
ymax = bd_rectangle[3]
b_rec = BoundingRectangle(xmin, xmax, ymin, ymax)

# Test plot sans contour
ax0 = b_rec.plot(color='b', linestyle='dotted')

#Test plot avec contour
ax = contour1.plot()
b_rec.plot(ax, linestyle='dotted')

# test area et centre
center = b_rec.center()
area = b_rec.area()
ax.scatter(center[0], center[1], color='r', marker='+')

line_fig2_seg1 = vme.LineSegment2D(vm.Point2D(0, 1), vm.Point2D(0.25, 0.5))
line_fig2_seg2 = vme.LineSegment2D(vm.Point2D(0.25, 0.5), vm.Point2D(1, 0.5))
line_fig2_seg3 = vme.LineSegment2D(vm.Point2D(1, 0.5), vm.Point2D(0.45, 0))
line_fig2_seg4 = vme.LineSegment2D(vm.Point2D(0.45, 0), vm.Point2D(1, -1))
line_fig2_seg5 = vme.LineSegment2D(vm.Point2D(1, -1), vm.Point2D(0, -0.5))
line_fig2_seg6 = vme.LineSegment2D(vm.Point2D(0, -0.5), vm.Point2D(-1, -1))
line_fig2_seg7 = vme.LineSegment2D(vm.Point2D(-1, -1), vm.Point2D(-0.45, 0))
line_fig2_seg8 = vme.LineSegment2D(vm.Point2D(-0.45, 0), vm.Point2D(-1, 0.5))
line_fig2_seg9 = vme.LineSegment2D(vm.Point2D(-1, 0.5), vm.Point2D(-0.25, 0.5))
line_fig2_seg10 = vme.LineSegment2D(vm.Point2D(-0.25, 0.5), vm.Point2D(0, 1))

contour2 = vmw.Contour2D([line_fig2_seg1, line_fig2_seg2, line_fig2_seg3, line_fig2_seg4,
                          line_fig2_seg5, line_fig2_seg6, line_fig2_seg7, line_fig2_seg8, line_fig2_seg9,
                          line_fig2_seg10])

bd_rectangle = contour2.bounding_rectangle()

xmin = bd_rectangle[0]
xmax = bd_rectangle[1]
ymin = bd_rectangle[2]
ymax = bd_rectangle[3]
b_rec2 = BoundingRectangle(xmin, xmax, ymin, ymax)

# Test plot sans contour
ax2 = b_rec2.plot()

## Test plot avec contour

ax3 = contour2.plot()
b_rec2.plot(ax3)
center2 = b_rec2.center()
ax3.scatter(center2[0], center2[1])

# test boundingbox intersection true
print(b_rec.b_rectangle_intersection(b_rec2))

# test boundingbox intersection flase
line_fig3_seg1 = vme.LineSegment2D(vm.Point2D(5, 1), vm.Point2D(5.25, 0.5))
line_fig3_seg2 = vme.LineSegment2D(vm.Point2D(5.25, 0.5), vm.Point2D(6, 0.5))
line_fig3_seg3 = vme.LineSegment2D(vm.Point2D(6, 0.5), vm.Point2D(5.45, 0))
line_fig3_seg4 = vme.LineSegment2D(vm.Point2D(5.45, 0), vm.Point2D(6, -1))
line_fig3_seg5 = vme.LineSegment2D(vm.Point2D(6, -1), vm.Point2D(5, -0.5))
line_fig3_seg6 = vme.LineSegment2D(vm.Point2D(5, -0.5), vm.Point2D(4, -1))
line_fig3_seg7 = vme.LineSegment2D(vm.Point2D(4, -1), vm.Point2D(4.55, 0))
line_fig3_seg8 = vme.LineSegment2D(vm.Point2D(4.55, 0), vm.Point2D(4, 0.5))
line_fig3_seg9 = vme.LineSegment2D(vm.Point2D(4, 0.5), vm.Point2D(4.75, 0.5))
line_fig3_seg10 = vme.LineSegment2D(vm.Point2D(4.75, 0.5), vm.Point2D(5, 1))

contour3 = vmw.Contour2D([line_fig3_seg1, line_fig3_seg2, line_fig3_seg3, line_fig3_seg4,
                          line_fig3_seg5, line_fig3_seg6, line_fig3_seg7, line_fig3_seg8, line_fig3_seg9,
                          line_fig3_seg10])
bd_rectangle = contour3.bounding_rectangle()

xmin = bd_rectangle[0]
xmax = bd_rectangle[1]
ymin = bd_rectangle[2]
ymax = bd_rectangle[3]
b_rec3 = BoundingRectangle(xmin, xmax, ymin, ymax)
ax4 = contour3.plot(ax3)
b_rec3.plot(ax4)
# Test intersection
print(b_rec.b_rectangle_intersection(b_rec2))
# test non intersection
print(b_rec.b_rectangle_intersection(b_rec3))
# test distance
print(b_rec.distance_to_b_rectangle(b_rec3))
# test distance avec intersection = 0
print(b_rec.distance_to_b_rectangle(b_rec2))  # d=0

# test distance to point
print(b_rec.distance_to_point([0.45, 0.45]))  #

print(b_rec.distance_to_point([2, 0.45]))  # d= 1

print(b_rec.distance_to_point([2, 2]))
