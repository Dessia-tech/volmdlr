import volmdlr
import volmdlr.wires as vmw
import volmdlr.edges as vme

line_seg1 = vme.LineSegment2D(volmdlr.Point2D(-0.5, -0.2), volmdlr.O2D)
line_seg2 = vme.LineSegment2D(volmdlr.O2D, volmdlr.Point2D(0.3, 1))
line_seg3 = vme.LineSegment2D(volmdlr.Point2D(0.3, 1), volmdlr.Point2D(1, 1))
line_seg4 = vme.LineSegment2D(volmdlr.Point2D(1, 1), volmdlr.Point2D(1, -0.5))
line_seg5 = vme.LineSegment2D(volmdlr.Point2D(1, -0.5), volmdlr.Point2D(-0.5, -0.2))

contour = vmw.Contour2D([line_seg1, line_seg2, line_seg3, line_seg4, line_seg5])
ax = contour.plot()

line1 = vme.Line2D(volmdlr.Point2D(-0.5, 1), volmdlr.O2D)
line1.plot(ax=ax, color='r')

cute_wire_line1 = contour.cut_by_line(line1)
for c in cute_wire_line1:
    c.plot()
assert len(cute_wire_line1) == 2

line2 = vme.Line2D(volmdlr.Point2D(-0.5, -0.5), volmdlr.O2D)
line2.plot(ax=ax, color='b')

cute_wire_line2 = contour.cut_by_line(line2)
for c in cute_wire_line2:
    c.plot()
assert len(cute_wire_line2) == 2

# for cntr in cute_wire_line2:
#     cut_by_line_results.extend(cntr.cut_by_line2(line1))
#     ax = line1.plot(color='r')
#     cntr.plot(ax=ax)
# for c in cut_by_line_results:
#     c.plot()
