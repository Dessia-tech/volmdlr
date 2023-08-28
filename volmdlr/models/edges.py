"""some volmdlr edges models."""
from geomdl import utilities
import volmdlr
from volmdlr import edges, curves

DEGREE = 3
points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
knotvector = utilities.generate_knot_vector(DEGREE, len(points))
knot_multiplicity = [1] * len(knotvector)
bspline1 = edges.BSplineCurve2D(DEGREE, points, knot_multiplicity, knotvector, None)

lineseg = edges.LineSegment2D(volmdlr.Point2D(0, 0.2), volmdlr.Point2D(3, -0.2))

arc = edges.Arc2D.from_3_points(volmdlr.Point2D(0, 0.3), volmdlr.Point2D(1, -0.3), volmdlr.Point2D(2, 2))
u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
frame = volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector)
ellipse2d = curves.Ellipse2D(2, 1, frame)
arc_ellipse2d = edges.ArcEllipse2D(ellipse2d, start=volmdlr.Point2D(0.5, 1.5), end=volmdlr.Point2D(1.5, 0.5))


def bspline_curve3d():
    """Bspline curve3d model."""
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
    return edges.BSplineCurve3D(degree=degree, control_points=control_points,
                                knot_multiplicities=knot_multiplicities,
                                knots=knots,
                                weights=weights,
                                name='B Spline Curve 3D 1')


vector1 = volmdlr.Vector3D(1, 1, 1)
vector1 = vector1.unit_vector()
vector2 = vector1.deterministic_unit_normal_vector()
vector3 = vector1.cross(vector2)


def arc3d():
    """Arc 3d model."""
    circle3d = curves.Circle3D(volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3), 1)
    return edges.Arc3D(circle3d, start=volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                       end=volmdlr.Point3D(-0.9855985596534886, -0.11957315586905026, -0.11957315586905026))


def arc_ellipse3d():
    """Arc ellipse 3d model."""
    ellipse3d = curves.Ellipse3D(2, 1, volmdlr.Frame3D(volmdlr.Point3D(1, 2, 1), vector3, vector1, vector2))
    return edges.ArcEllipse3D(ellipse3d,
                              start=volmdlr.Point3D(0.42264973081037405, 1.4226497308103743, 0.42264973081037427),
                              end=volmdlr.Point3D(1.577350269189626, 2.5773502691896253, 1.5773502691896257))


def linesegment3d():
    """Linesegment 3d model."""
    return edges.LineSegment3D(volmdlr.Point3D(1, 2, 4), volmdlr.Point3D(-1, 5, -3))
