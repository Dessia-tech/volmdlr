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
bdr = BoundingRectangle(xmin, xmax, ymin, ymax)
# Test plot sans contour
ax0 = bdr.plot()

## Test plot avec contour
ax = contour1.plot()
bdr.plot(ax)

# test area et centre
center = bdr.center()
area = bdr.area()
ax.scatter(center[0], center[1])


line_fig2_seg1 = vme.LineSegment2D(vm.Point2D(0, 1), vm.Point2D(0.25, 0.5))
line_fig2_seg2 = vme.LineSegment2D(vm.Point2D(0.25, 0.5), vm.Point2D(1, 0.5))
line_fig2_seg3 = vme.LineSegment2D(vm.Point2D(1, 0.5), vm.Point2D(0.45, 0))
line_fig2_seg4 = vme.LineSegment2D(vm.Point2D(0.45, 0), vm.Point2D(1, -1))
line_fig2_seg5 = vme.LineSegment2D(vm.Point2D(1, -1), vm.Point2D(0, -0.5))
line_fig2_seg6 = vme.LineSegment2D(vm.Point2D(0, -0.5), vm.Point2D(-1, -1))
line_fig2_seg7 = vme.LineSegment2D(vm.Point2D(-1, -1),vm.Point2D(-0.45, 0))
line_fig2_seg8 = vme.LineSegment2D(vm.Point2D(-0.45, 0),vm.Point2D(-1, 0.5))
line_fig2_seg9 = vme.LineSegment2D(vm.Point2D(-1, 0.5),vm.Point2D(-0.25, 0.5))
line_fig2_seg10 = vme.LineSegment2D(vm.Point2D(-0.25, 0.5),vm.Point2D(0, 1))

contour2 = vmw.Contour2D([line_fig2_seg1, line_fig2_seg2, line_fig2_seg3, line_fig2_seg4,
                          line_fig2_seg5, line_fig2_seg6, line_fig2_seg7, line_fig2_seg8, line_fig2_seg9, line_fig2_seg10 ])

bd_rectangle = contour2.bounding_rectangle()

xmin = bd_rectangle[0]
xmax = bd_rectangle[1]
ymin = bd_rectangle[2]
ymax = bd_rectangle[3]
bdr2 = BoundingRectangle(xmin, xmax, ymin, ymax)
# Test plot sans contour
ax2 = bdr2.plot()

## Test plot avec contour

ax3 = contour2.plot()
ax3 = bdr2.plot(ax3)
center2 = bdr2.center()
ax3.scatter(center2[0], center2[1])
