"""
Module to translate objects in Volmdlr to OCP.
"""
from OCP.Geom import (Geom_BSplineSurface, Geom_CylindricalSurface, Geom_ConicalSurface, Geom_ToroidalSurface,
                      Geom_SphericalSurface, Geom_Plane, Geom_BSplineCurve, Geom_Line)
from OCP.Geom2d import Geom2d_BSplineCurve
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
# from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt, TColgp_Array1OfPnt2d
from OCP.gp import gp_Pnt, gp_Ax3, gp_Dir, gp_Pnt2d, gp_Dir2d


def list_to_tcolstd_array1ofinteger(list_of_int):
    """
    Helper function.

    """
    array = TColStd_Array1OfInteger(1, len(list_of_int))
    for i, value in enumerate(list_of_int):
        array.SetValue(i + 1, value)
    return array


def list_to_tcolstd_array1ofreal(list_of_real):
    """
    Helper function.

    """
    array = TColStd_Array1OfReal(1, len(list_of_real))
    for i, value in enumerate(list_of_real):
        array.SetValue(i + 1, value)
    return array


def list_to_tcolgp_array10fpnt(list_of_points):
    """
    Helper function.

    """
    array = TColgp_Array1OfPnt(1, len(list_of_points))
    for i, point in enumerate(list_of_points):
        array.SetValue(i + 1, point3d_to_ocp(point))
    return array


def list_to_tcolgp_array10fpnt2d(list_of_points):
    """
    Helper function.

    """
    array = TColgp_Array1OfPnt2d(1, len(list_of_points))
    for i, point in enumerate(list_of_points):
        array.SetValue(i + 1, point2d_to_ocp(point))
    return array


def point2d_to_ocp(point2d):
    """
    Create an OCP Point2D from a voldmlr Point2D.

    :param point2d: volmdlr Point2D.
    :return: OCP Point2D
    """
    return gp_Pnt2d(*point2d)


def point3d_to_ocp(point):
    """
    Create an OCP Point3D from a voldmlr Point3D.

    :param point: volmdlr Point3D.
    :return: OCP Point3D
    """
    return gp_Pnt(*point)


def vector3d_to_ocp(vector):
    """
    Create an OCP Vector3D from a voldmlr Vector3D.

    :param vector: volmdlr Vector3D.
    :return: OCP Vector3D
    """
    return gp_Dir(*vector)


def vector2d_to_ocp(vector):
    """
    Create an OCP Vector2D from a voldmlr Vector2D.

    :param vector: volmdlr Vector2D.
    :return: OCP Vector2D
    """
    return gp_Dir2d(*vector)


def frame3d_to_ocp(frame):
    """
    Create an OCP Frame3D from a voldmlr Frame3D.

    :param frame: volmdlr Frame3D.
    :return: OCP Frame3D.
    """
    point = point3d_to_ocp(frame.origin)
    z_vector = vector3d_to_ocp(frame.w)
    x_vector = vector3d_to_ocp(frame.u)
    return gp_Ax3(point, z_vector, x_vector)


def line3d_to_ocp(line):
    """
    Create an OCP Line from a voldmlr Line3D.

    :param line: volmdlr Line3D.
    :return: OCP Frame3D.
    """
    point = point3d_to_ocp(line.point1)
    direction = vector3d_to_ocp(line.unit_direction_vector())
    return Geom_Line(point, direction)


def bsplinecurve3d_to_ocp(bsplinecurve):
    """
    Creates a Bspline Curve 3D from a volmdlr object.

    :param bsplinecurve: volmdlr BSpline Curve 3D.
    :return:
    """
    if bsplinecurve.weights is None:
        return Geom_BSplineCurve(list_to_tcolgp_array10fpnt(bsplinecurve.control_points),
                                 list_to_tcolstd_array1ofreal(bsplinecurve.knots),
                                 list_to_tcolstd_array1ofinteger(bsplinecurve.knot_multiplicities),
                                 bsplinecurve.degree)
    return Geom_BSplineCurve(list_to_tcolgp_array10fpnt(bsplinecurve.control_points),
                             list_to_tcolstd_array1ofreal(bsplinecurve.weights),
                             list_to_tcolstd_array1ofreal(bsplinecurve.knots),
                             list_to_tcolstd_array1ofinteger(bsplinecurve.knot_multiplicities),
                             bsplinecurve.degree)


def bsplinecurve2d_to_ocp(bsplinecurve):
    """
    Creates a Bspline Curve 3D from a volmdlr object.

    :param bsplinecurve: volmdlr BSpline Curve 2D.
    :return:
    """
    if not bsplinecurve.weights:
        return Geom2d_BSplineCurve(list_to_tcolgp_array10fpnt2d(bsplinecurve.control_points),
                                   list_to_tcolstd_array1ofreal(bsplinecurve.knots),
                                   list_to_tcolstd_array1ofinteger(bsplinecurve.knot_multiplicities),
                                   bsplinecurve.degree)
    return Geom2d_BSplineCurve(list_to_tcolgp_array10fpnt2d(bsplinecurve.control_points),
                               list_to_tcolstd_array1ofreal(bsplinecurve.weights),
                               list_to_tcolstd_array1ofreal(bsplinecurve.knots),
                               list_to_tcolstd_array1ofinteger(bsplinecurve.knot_multiplicities),
                               bsplinecurve.degree)


def plane_to_ocp(surface):
    """
    Create an OCP Plane from a voldmlr Plane3D.

    :param surface: volmdlr Plane3D.
    :return: OCP Plane.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_Plane(gp_ax3)


def cylindricalsurface_to_ocp(surface):
    """
    Create an OCP CylindricalSurface from a voldmlr CylindricalSurface3D.

    :param surface: volmdlr CylindricalSurface3D.
    :return: OCP CylindricalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_CylindricalSurface(gp_ax3, surface.radius)


def conicalsurface_to_ocp(surface):
    """
    Create an OCP ConicalSurface from a voldmlr ConicalSurface3D.

    :param surface: volmdlr ConicalSurface3D.
    :return: OCP ConicalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_ConicalSurface(gp_ax3, surface.semi_angle, surface.ref_radius)


def sphericalsurface_to_ocp(surface):
    """
    Create an OCP SphericalSurface from a voldmlr SphericalSurface3D.

    :param surface: volmdlr SphericalSurface3D.
    :return: OCP SphericalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_SphericalSurface(gp_ax3, surface.radius)


def toroidalsurface_to_ocp(surface):
    """
    Create an OCP ToroidalSurface from a voldmlr ToroidalSurface3D.

    :param surface: volmdlr ToroidalSurface3D.
    :return: OCP ToroidalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_ToroidalSurface(gp_ax3, surface.major_radius, surface.minor_radius)


def bsplinesurface_to_ocp(surface):
    """
    Create an OCP BSplineSurface3D from a voldmlr FraBSplineSurface3Dme3D.

    :param surface: volmdlr BSplineSurface3D.
    :return: OCP BSplineSurface.
    """
    points = surface.control_points_table
    u_deg, v_deg = surface.degree_u, surface.degree_v

    poles = TColgp_Array2OfPnt(1, len(points), 1, len(points[0]))
    for i, row in enumerate(points):
        for j, point in enumerate(row):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*point))

    uknots = TColStd_Array1OfReal(1, len(surface.u_knots))
    for i, value in enumerate(surface.u_knots):
        uknots.SetValue(i + 1, value)

    vknots = TColStd_Array1OfReal(1, len(surface.v_knots))
    for i, value in enumerate(surface.v_knots):
        vknots.SetValue(i + 1, value)

    umult = TColStd_Array1OfInteger(1, len(surface.u_multiplicities))
    for i, value in enumerate(surface.u_multiplicities):
        umult.SetValue(i + 1, value)

    vmult = TColStd_Array1OfInteger(1, len(surface.v_multiplicities))
    for i, value in enumerate(surface.v_multiplicities):
        vmult.SetValue(i + 1, value)

    return Geom_BSplineSurface(poles, uknots, vknots, umult, vmult, u_deg, v_deg, False, False)
