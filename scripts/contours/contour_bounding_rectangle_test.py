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

#ax = contour1.plot()

#for point in [p1, p2, p3, p4]:
#    point.plot(ax, color = 'r')

bdr = BoundingRectangle(xmin, xmax, ymin, ymax)

bdr.plot()






