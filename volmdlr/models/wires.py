"""
Wires model.
"""

import volmdlr
from volmdlr import curves, edges, wires

point1 = volmdlr.Point3D(0.009988247648354, 0.021515325970206, 0.01)
point2 = volmdlr.Point3D(-0.005370502813076001, 0.030435502328536, 0.01)
point3 = volmdlr.Point3D(-0.01886730653873, 0.021842697827615002, 0.01)
point4 = volmdlr.Point3D(-0.007181504688291, 0.0036374839872410003, 0.01)

circle_frame = volmdlr.OXYZ.copy()
circle_frame.origin = volmdlr.Point3D(0.0, 0.022, 0.01)
circle = curves.Circle3D(circle_frame, 0.01)

arc = edges.Arc3D(circle, point1, point2)
linesegment = edges.LineSegment3D(point2, point3)

ellipse_frame = volmdlr.Frame3D(origin=volmdlr.Point3D(-0.01305821955846, 0.012718306973665, 0.01),
                                u=volmdlr.Vector3D(-0.5370502813077209, 0.8435502328535615, 0.0),
                                v=volmdlr.Vector3D(-0.8435502328535615, -0.5370502813077209, 0.0),
                                w=volmdlr.Vector3D(0.0, 0.0, 1.0))

ellipse = curves.Ellipse3D(frame=ellipse_frame, major_axis=0.010816653826392, minor_axis=0.009)

arc_ellipse = edges.ArcEllipse3D(ellipse, point3, point4)

control_points = [point4,
                  volmdlr.Point3D(-0.0008453667619050001, 0.006902088714286, 0.01),
                  volmdlr.Point3D(0.006404972476190001, 0.007439177428571, 0.01),
                  volmdlr.Point3D(0.010049324, 0.005650317, 0.01),
                  volmdlr.Point3D(0.010123425, 0.014348723, 0.01)]
knot_multiplicities = [4, 1, 4]

knots = [0., 0.5, 1.]

bspline = edges.BSplineCurve3D(3, control_points, knot_multiplicities, knots)

wire3d_all_edges = wires.Wire3D([arc, linesegment, arc_ellipse, bspline])