"""some volmdlr edges models."""
from geomdl import utilities
import volmdlr
from volmdlr import edges, curves

degree = 3
points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
knotvector = utilities.generate_knot_vector(degree, len(points))
knot_multiplicity = [1] * len(knotvector)
bspline1 = edges.BSplineCurve2D(degree, points, knot_multiplicity, knotvector, None, False)

lineseg = edges.LineSegment2D(volmdlr.Point2D(0, 0.2), volmdlr.Point2D(3, -0.2))

arc = edges.Arc2D(volmdlr.Point2D(0, 0.3), volmdlr.Point2D(1, -0.3), volmdlr.Point2D(2, 2))
u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
frame = volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector)
ellipse2d = curves.Ellipse2D(2, 1, frame)
arc_ellipse2d = edges.ArcEllipse2D(ellipse2d, start=volmdlr.Point2D(0.5, 1.5), end=volmdlr.Point2D(1.5, 0.5))

degree = 5
control_points = [volmdlr.Point3D(0, 3, 0),
                  volmdlr.Point3D(3, 2, 1),
                  volmdlr.Point3D(5, -1, 4),
                  volmdlr.Point3D(5, -4, 0),
                  volmdlr.Point3D(-1, -2, -3),
                  volmdlr.Point3D(-3, 4, 1)]
knots = [0.0, 1.0]
knot_multiplicities = [6, 6]
weights = None  # [1, 2, 1, 2, 1, 2]
bspline_curve3d = edges.BSplineCurve3D(degree=degree, control_points=control_points,
                                       knot_multiplicities=knot_multiplicities,
                                       knots=knots,
                                       weights=weights,
                                       periodic=False,
                                       name='B Spline Curve 3D 1')
