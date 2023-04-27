"""some volmdlr edges models."""
from geomdl import utilities
import volmdlr
from volmdlr import edges

degree = 3
points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
knotvector = utilities.generate_knot_vector(degree, len(points))
knot_multiplicity = [1] * len(knotvector)
bspline1 = edges.BSplineCurve2D(degree, points, knot_multiplicity, knotvector, None, False)

lineseg = edges.LineSegment2D(volmdlr.Point2D(0, 0.2), volmdlr.Point2D(3, -0.2))

arc = edges.Arc2D(volmdlr.Point2D(0, 0.3), volmdlr.Point2D(1, -0.3), volmdlr.Point2D(2, 2))

arc_ellipse2d = edges.ArcEllipse2D(start=10 * volmdlr.Point2D(-0.125, -0.08416500663326211),
                                   interior=10 * volmdlr.Point2D(-0.03543560762586048, -0.011930639375832372),
                                   end=10 * volmdlr.Point2D(0.0, 0.125), center=10 * volmdlr.Point2D(-0.15, 0.125),
                                   major_dir=volmdlr.Vector2D(0, 1))
