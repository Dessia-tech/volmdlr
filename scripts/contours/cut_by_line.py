import volmdlr
import volmdlr.wires as vmw
import volmdlr.edges as vme
import matplotlib.pyplot as plt
import numpy as np

line_seg1 = vme.LineSegment2D(volmdlr.Point2D(-0.5, -0.2), volmdlr.O2D)
line_seg2 = vme.LineSegment2D(volmdlr.O2D, volmdlr.Point2D(0.3, 1))
line_seg3 = vme.LineSegment2D(volmdlr.Point2D(0.3, 1), volmdlr.Point2D(1, 1))
line_seg4 = vme.LineSegment2D(volmdlr.Point2D(1, 1), volmdlr.Point2D(1, -0.5))
line_seg5 = vme.LineSegment2D(volmdlr.Point2D(1, -0.5), volmdlr.Point2D(-0.5, -0.2))

contour1 = vmw.Contour2D([line_seg1, line_seg2, line_seg3, line_seg4, line_seg5])

line1 = vme.Line2D(volmdlr.Point2D(-0.5, 1), volmdlr.O2D)

cute_wire_line1 = contour1.cut_by_line(line1)
assert len(cute_wire_line1) == 2

line2 = vme.Line2D(volmdlr.Point2D(-0.5, -0.5), volmdlr.O2D)
cute_wire_line2 = contour1.cut_by_line(line2)
assert len(cute_wire_line2) == 2

line_segment1 = vme.LineSegment2D(volmdlr.Point2D(1, -1), volmdlr.Point2D(2, 2))
line_segment2 = vme.LineSegment2D(volmdlr.Point2D(2, 2), volmdlr.Point2D(-2, 1))
line_segment3 = vme.LineSegment2D(volmdlr.Point2D(-2, 1), volmdlr.Point2D(-2, 0.7))
lie_segment4 = vme.LineSegment2D(volmdlr.Point2D(-2, 0.7), volmdlr.Point2D(-1, 1))
points2d = [volmdlr.Point2D(-1, 1),
            volmdlr.Point2D(2, 2),
            volmdlr.Point2D(-2, -2),
            volmdlr.Point2D(1, -1)]
bspline = vme.BSplineCurve2D(3, points2d, knot_multiplicities=[4, 4],
                               knots=[0.0, 1.0])
bspline_middle_point = bspline.point_at_abscissa(bspline.length()*0.5)
bspline_tangent = bspline.tangent(0.5)
infinit_line1 = vme.Line2D(bspline_middle_point,
                           bspline_tangent)
infinit_line2 = vme.Line2D(bspline.point_at_abscissa(bspline.length()*0.73),
                           volmdlr.Point2D(-0.8, 1))
infinit_line3 = vme.Line2D(bspline_middle_point,
                           volmdlr.Point2D(2, 2))

contour2 = vmw.Contour2D([bspline, line_segment1, line_segment2, line_segment3, lie_segment4])

cut_contour_by_line1 = contour2.cut_by_line(infinit_line1)
assert len(cut_contour_by_line1) == 2


cut_contour_by_line2 = contour2.cut_by_line(infinit_line2)
assert len(cut_contour_by_line2) == 2

cut_contour_by_line3 = contour2.cut_by_line(infinit_line3)
assert len(cut_contour_by_line3) == 2

list_contours = [contour1, contour2, contour2, contour2]
lines = [line2, infinit_line1, infinit_line2, infinit_line3]
lists_cutted_contours = [cute_wire_line2, cut_contour_by_line1,
                         cut_contour_by_line2, cut_contour_by_line3]
fig, axs = plt.subplots(4, 3)
for i in range(0, 4):
    list_contours[i].plot(ax=axs[i][0])
    for j in range(0, 3):
        lines[i].plot(ax=axs[i][j], color='r')
        if j != 0:
            r, g, b = np.random.random(), np.random.random(), np.random.random()
            lists_cutted_contours[i][j-1].plot(ax=axs[i][j], color=(r, g, b))
