"""
Module to translate objects in OCP to Volmdlr.
"""
from OCP.Geom import Geom_Circle, Geom_Line, Geom_BSplineCurve, Geom_Ellipse
from OCP.TColStd import TColStd_Array2OfReal, TColStd_Array1OfReal

import volmdlr
from volmdlr import curves, edges, surfaces


def point2d_from_occt(occt_point):
    """
    Instanciates a volmdlr Point2D, from occt object.

    :param occt_point: OCCT point.
    :return: volmdlr Point2D.
    """
    return volmdlr.Point2D(occt_point.X(), occt_point.Y())


def vector2d_from_occt(occt_vector):
    """
    Instanciates a volmdlr Vector2D, from occt object.

    :param occt_vector: OCCT Vector.
    :return: volmdlr Vector2D.
    """
    return volmdlr.Vector2D(occt_vector.X(), occt_vector.Y())


def point3d_from_occt(occt_point):
    """
    Instanciates a volmdlr Point3D, from occt object.

    :param occt_point: OCCT point.
    :return: volmdlr Point3D.
    """
    return volmdlr.Point3D(occt_point.X(), occt_point.Y(), occt_point.Z())


def vector3d_from_occt(occt_vector):
    """
    Instanciates a volmdlr Vector3D, from occt object.

    :param occt_vector: OCCT Vector.
    :return: volmdlr Vector3D.
    """
    return volmdlr.Vector3D(occt_vector.X(), occt_vector.Y(), occt_vector.Z())


def frame3d_from_occt_ax3(ax2):
    """
    Instanciates a volmdlr Frame3D, from occt object.

    :param ax2: OCCT ax2.
    :return: volmdlr Frame3D.
    """
    origin = point3d_from_occt(ax2.Location())
    u = vector3d_from_occt(ax2.XDirection())
    v = vector3d_from_occt(ax2.YDirection())
    return volmdlr.Frame3D(origin, u, v, u.cross(v))


def frame2d_from_occt_ax22d(ax22d):
    """
    Instanciates a volmdlr Frame2D, from occt object.

    :param ax22d: OCCT ax2 2d.
    :return: volmdlr Frame2D.
    """
    origin = point2d_from_occt(ax22d.Location())
    u = vector2d_from_occt(ax22d.XDirection())
    v = vector2d_from_occt(ax22d.YDirection())
    return volmdlr.Frame2D(origin, u, v)


# Curves


def line2d_from_occt(occt_line):
    """
    Instanciates a volmdlr Line2D, from occt object.

    :param occt_line: OCCT Line.
    :return: volmdlr Line2D.
    """
    vector = vector2d_from_occt(occt_line.Direction())
    point = point2d_from_occt(occt_line.Location())
    return curves.Line2D.from_point_and_vector(point, vector)


def line3d_from_occt(occt_line):
    """
    Instanciates a volmdlr Line3D, from occt object.

    :param occt_line: OCCT Line 3D.
    :return: volmdlr Line3D.
    """
    position = occt_line.Position()
    return curves.Line3D.from_point_and_vector(point3d_from_occt(position.Location()),
                                               vector3d_from_occt(position.Direction()))


def circle2d_from_occt(curve):
    """
    Instanciates a volmdlr Circle2D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr Circle2D.
    """
    frame = frame2d_from_occt_ax22d(curve.Circ2d().Axis())
    return curves.Circle2D(frame, curve.Radius())


def circle3d_from_occt(curve):
    """
    Instanciates a volmdlr Circle3D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr Circle3D.
    """
    frame = frame3d_from_occt_ax3(curve.Position())
    return curves.Circle3D(frame, curve.Radius())


def ellipse2d_from_occt(curve):
    """
    Instanciates a volmdlr Ellipse2D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr Ellipse2D.
    """
    frame = frame2d_from_occt_ax22d(curve.Position())
    return curves.Ellipse2D(frame, curve.MajorRadius(), curve.MinorRadius())


def ellipse3d_from_occt(curve):
    """
    Instanciates a volmdlr Ellipse3D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr Ellipse3D.
    """
    frame = frame3d_from_occt_ax3(curve.Position())
    return curves.Ellipse3D(curve.MajorRadius(), curve.MinorRadius(), frame)


def hyperbola3d_from_occt(curve):
    """
    Instanciates a volmdlr Hyperbola3D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr Hyperbola3D.
    """
    frame = frame3d_from_occt_ax3(curve.Position())
    return curves.Hyperbola3D(frame, curve.MajorRadius(), curve.MinorRadius())


def parabola3d_from_occt(curve):
    """
    Instanciates a volmdlr Parabola3D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr Parabola3D.
    """
    frame = frame3d_from_occt_ax3(curve.Position())
    frame = volmdlr.Frame3D(frame.origin, frame.v, frame.u, frame.w)
    return curves.Parabola3D(frame, curve.Focal())


# Edges


def bsplinecurve2d_from_occt(curve):
    """
    Instanciates a volmdlr BSplineCurve2D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr BSplineCurve2D.
    """
    control_points = [volmdlr.Point2D(point.X(), point.Y()) for point in curve.Poles()]
    knots = list(curve.Knots())
    multiplicities = list(curve.Multiplicities())
    weigths = None
    if curve.IsRational():
        curve.Weights(weights_array := TColStd_Array1OfReal(1, len(control_points)))
        weigths = list(weights_array)
    return edges.BSplineCurve2D(curve.Degree(), control_points, multiplicities, knots, weigths)


def bsplinecurve3d_from_occt(curve):
    """
    Instanciates a volmdlr BSplineCurve3D, from occt object.

    :param curve: OCCT curve.
    :return: volmdlr BSplineCurve3D.
    """
    control_points = [volmdlr.Point3D(point.X(), point.Y(), point.Z()) for point in curve.Poles()]
    knots = list(curve.Knots())
    multiplicities = list(curve.Multiplicities())
    weigths = None
    if curve.IsRational():
        curve.Weights(weights_array := TColStd_Array1OfReal(1, len(control_points)))
        weigths = list(weights_array)
    return edges.BSplineCurve3D(curve.Degree(), control_points, multiplicities, knots, weigths)


OCCT_TO_VOLMDLR = {Geom_Line: line3d_from_occt,
                   Geom_Circle: circle3d_from_occt,
                   Geom_Ellipse: ellipse3d_from_occt,
                   Geom_BSplineCurve: bsplinecurve3d_from_occt}


def volmdlr_edge_from_occt_curve(occt_curve, first, last, orientation):
    """
    Instanciate a volmdlr edge form an occt curve.

    :param occt_curve: occt curve.
    :param first: first point.
    :param last: last point.
    :param orientation: orientation of the curve to be considered.
    :return: Volmdlr trimmed edge.
    """
    start = point3d_from_occt(occt_curve.Value(first))
    end = point3d_from_occt(occt_curve.Value(last))
    function = OCCT_TO_VOLMDLR[occt_curve.__class__]
    curve = function(occt_curve)
    return curve.trim(start, end, same_sense=orientation == 0)


def trimmedcurve3d_from_occt(occt_curve):
    """
    Intanciate an edge from a trimmed curve in OCCT.

    :param occt_curve: occt trimmed curve.
    :return: Volmdlr trimmed edge.
    """
    start_point = point3d_from_occt(occt_curve.StartPoint())
    end_point = point3d_from_occt(occt_curve.EndPoint())
    occt_basis_curve = occt_curve.BasisCurve()
    volmdlr_curve = OCCT_TO_VOLMDLR[occt_basis_curve.__class__](occt_basis_curve)
    return volmdlr_curve.trim(start_point, end_point)


# Surfaces


def sphericalsurface_from_occt(occt_surface):
    """
    Instanciates a volmdlr SphericalSurface3D, from occt object.

    :param occt_surface: OCCT surface.
    :return: volmdlr SphericalSurface3D.
    """
    frame = frame3d_from_occt_ax3(occt_surface.Position())
    radius = occt_surface.Radius()
    return surfaces.SphericalSurface3D(frame, radius)


def cylindricalsurface_from_occt(occt_surface):
    """
    Instanciates a volmdlr CylindricalSurface3D, from occt object.

    :param occt_surface: OCCT surface.
    :return: volmdlr CylindricalSurface3D.
    """
    frame = frame3d_from_occt_ax3(occt_surface.Position())
    radius = occt_surface.Radius()
    return surfaces.CylindricalSurface3D(frame, radius)


def plane_from_occt(occt_surface):
    """
    Instanciates a volmdlr Plane3D, from occt object.

    :param occt_surface: OCCT surface.
    :return: volmdlr Plane3D.
    """
    frame = frame3d_from_occt_ax3(occt_surface.Position())
    return surfaces.Plane3D(frame)


def toroidalsurface_from_occt(occt_surface):
    """
    Instanciates a volmdlr ToroidalSurface3D, from occt object.

    :param occt_surface: OCCT surface.
    :return: volmdlr ToroidalSurface3D.
    """
    frame = frame3d_from_occt_ax3(occt_surface.Position())
    return surfaces.ToroidalSurface3D(frame, occt_surface.MajorRadius(), occt_surface.MinorRadius())


def conicalsurface_from_occt(occt_surface):
    """
    Instanciates a volmdlr ConicalSurface3D, from occt object.

    :param occt_surface: OCCT surface.
    :return: volmdlr ConicalSurface3D.
    """
    frame = frame3d_from_occt_ax3(occt_surface.Position())
    radius = occt_surface.RefRadius()
    semi_angle = occt_surface.SemiAngle()
    return surfaces.ConicalSurface3D(frame, semi_angle, radius)


def bsplinesurface_from_occt(occt_surface):
    """
    Instanciates a volmdlr BSplineSurface3D, from occt object.

    :param occt_surface: OCCT surface.
    :return: volmdlr BSplineSurface3D.
    """
    array = occt_surface.Poles()
    nb_v = occt_surface.NbVPoles()
    nb_u = occt_surface.NbUPoles()

    control_points = [point3d_from_occt(array.Value(i + 1, j + 1)) for i in range(nb_u) for j in range(nb_v)]
    u_knots = list(occt_surface.UKnots())
    u_multiplicities = list(occt_surface.UMultiplicities())
    v_knots = list(occt_surface.VKnots())
    v_multiplicities = list(occt_surface.VMultiplicities())
    weights = None
    if occt_surface.IsURational() or occt_surface.IsVRational():
        occt_surface.Weights(weights_array := TColStd_Array2OfReal(1, nb_u, 1, nb_v))

        weights = [weights_array.Value(i + 1, j + 1) for i in range(nb_u) for j in range(nb_v)]
    return surfaces.BSplineSurface3D(occt_surface.UDegree(), occt_surface.VDegree(), control_points, nb_u, nb_v,
                                     u_multiplicities, v_multiplicities, u_knots, v_knots, weights)
