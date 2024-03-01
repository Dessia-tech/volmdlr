from OCP.Geom import Geom_Plane, Geom_Line
from OCP.gp import gp_Pnt, gp_Ax3, gp_Dir, gp_Pnt2d, gp_Dir2d
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge
from OCP.Precision import Precision
from OCP.GeomAPI import GeomAPI_IntSS
# from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt, TColgp_Array1OfPnt2d
from OCP.Geom import (Geom_BSplineSurface, Geom_CylindricalSurface, Geom_ConicalSurface, Geom_ToroidalSurface,
                      Geom_SphericalSurface, Geom_Plane, Geom_BSplineCurve, Geom_Line)
from OCP.Geom2d import Geom2d_BSplineCurve


def list_to_tcolstd_array1ofinteger(list_of_int):
    array = TColStd_Array1OfInteger(1, len(list_of_int))
    for i, value in enumerate(list_of_int):
        array.SetValue(i + 1, value)
    return array


def list_to_tcolstd_array1ofreal(list_of_real):
    array = TColStd_Array1OfReal(1, len(list_of_real))
    for i, value in enumerate(list_of_real):
        array.SetValue(i + 1, value)
    return array


def list_to_tcolgp_array10fpnt(list_of_points):
    array = TColgp_Array1OfPnt(1, len(list_of_points))
    for i, point in enumerate(list_of_points):
        array.SetValue(i + 1, volmdlr_point3d_to_occt(point))
    return array


def list_to_tcolgp_array10fpnt2d(list_of_points):
    array = TColgp_Array1OfPnt2d(1, len(list_of_points))
    for i, point in enumerate(list_of_points):
        array.SetValue(i + 1, volmdlr_point2d_to_occt(point))
    return array


def volmdlr_point2d_to_occt(volmdlr_point2d):
    """
    Create an OCCT Point2D from a voldmlr Point2D.

    :param volmdlr_point2d: volmdlr Point2D.
    :return: OCCT Point2D
    """
    return gp_Pnt2d(*volmdlr_point2d)


def volmdlr_point3d_to_occt(volmdlr_point):
    """
    Create an OCCT Point3D from a voldmlr Point3D.

    :param volmdlr_point: volmdlr Point3D.
    :return: OCCT Point3D
    """
    return gp_Pnt(*volmdlr_point)


def volmdlr_vector3d_to_occt(volmdlr_vector):
    """
    Create an OCCT Vector3D from a voldmlr Vector3D.

    :param volmdlr_vector: volmdlr Vector3D.
    :return: OCCT Vector3D
    """
    return gp_Dir(*volmdlr_vector)


def volmdlr_vector2d_to_occt(volmdlr_vector):
    """
    Create an OCCT Vector2D from a voldmlr Vector2D.

    :param volmdlr_vector: volmdlr Vector2D.
    :return: OCCT Vector2D
    """
    return gp_Dir2d(*volmdlr_vector)


def volmdlr_frame3d_to_occt(volmdlr_frame):
    """
    Create an OCCT Frame3D from a voldmlr Frame3D.

    :param volmdlr_frame: volmdlr Frame3D.
    :return: OCCT Frame3D.
    """
    point = volmdlr_point3d_to_occt(volmdlr_frame.origin)
    z_vector = volmdlr_vector3d_to_occt(volmdlr_frame.w)
    x_vector = volmdlr_vector3d_to_occt(volmdlr_frame.u)
    return gp_Ax3(point, z_vector, x_vector)


def volmdlr_line3d_to_occt(volmdlr_line):
    """
    Create an OCCT Line from a voldmlr Line3D.

    :param volmdlr_line: volmdlr Line3D.
    :return: OCCT Frame3D.
    """
    point = volmdlr_point3d_to_occt(volmdlr_line.point1)
    direction = volmdlr_vector3d_to_occt(volmdlr_line.unit_direction_vector())
    return Geom_Line(point, direction)


def volmdlr_bsplinecurve3d_to_occt(volmdlr_bsplinecurve):
    if volmdlr_bsplinecurve.weights is None:
        return Geom_BSplineCurve(list_to_tcolgp_array10fpnt(volmdlr_bsplinecurve.control_points),
                                 list_to_tcolstd_array1ofreal(volmdlr_bsplinecurve.knots),
                                 list_to_tcolstd_array1ofinteger(volmdlr_bsplinecurve.knot_multiplicities),
                                 volmdlr_bsplinecurve.degree)
    return Geom_BSplineCurve(list_to_tcolgp_array10fpnt(volmdlr_bsplinecurve.control_points),
                             list_to_tcolstd_array1ofreal(volmdlr_bsplinecurve.weights),
                             list_to_tcolstd_array1ofreal(volmdlr_bsplinecurve.knots),
                             list_to_tcolstd_array1ofinteger(volmdlr_bsplinecurve.knot_multiplicities),
                             volmdlr_bsplinecurve.degree)


def volmdlr_bsplinecurve2d_to_occt(volmdlr_bsplinecurve):
    if not volmdlr_bsplinecurve.weights:
        return Geom2d_BSplineCurve(list_to_tcolgp_array10fpnt2d(volmdlr_bsplinecurve.control_points),
                                   list_to_tcolstd_array1ofreal(volmdlr_bsplinecurve.knots),
                                   list_to_tcolstd_array1ofinteger(volmdlr_bsplinecurve.knot_multiplicities),
                                   volmdlr_bsplinecurve.degree)
    return Geom2d_BSplineCurve(list_to_tcolgp_array10fpnt2d(volmdlr_bsplinecurve.control_points),
                               list_to_tcolstd_array1ofreal(volmdlr_bsplinecurve.weights),
                               list_to_tcolstd_array1ofreal(volmdlr_bsplinecurve.knots),
                               list_to_tcolstd_array1ofinteger(volmdlr_bsplinecurve.knot_multiplicities),
                               volmdlr_bsplinecurve.degree)


def volmdlr_plane_to_occt(volmdlr_surface):
    """
    Create an OCCT Plane from a voldmlr Plane3D.

    :param volmdlr_surface: volmdlr Plane3D.
    :return: OCCT Plane.
    """
    gp_ax3 = volmdlr_frame3d_to_occt(volmdlr_surface.frame)
    return Geom_Plane(gp_ax3)


def volmdlr_cylindricalsurface_to_occt(volmdlr_surface):
    """
    Create an OCCT CylindricalSurface from a voldmlr CylindricalSurface3D.

    :param volmdlr_surface: volmdlr CylindricalSurface3D.
    :return: OCCT CylindricalSurface.
    """
    gp_ax3 = volmdlr_frame3d_to_occt(volmdlr_surface.frame)
    return Geom_CylindricalSurface(gp_ax3, volmdlr_surface.radius)


def volmdlr_conicalsurface_to_occt(volmdlr_surface):
    """
    Create an OCCT ConicalSurface from a voldmlr ConicalSurface3D.

    :param volmdlr_surface: volmdlr ConicalSurface3D.
    :return: OCCT ConicalSurface.
    """
    gp_ax3 = volmdlr_frame3d_to_occt(volmdlr_surface.frame)
    return Geom_ConicalSurface(gp_ax3, volmdlr_surface.semi_angle, 0)


def volmdlr_sphericalsurface_to_occt(volmdlr_surface):
    """
    Create an OCCT SphericalSurface from a voldmlr SphericalSurface3D.

    :param volmdlr_surface: volmdlr SphericalSurface3D.
    :return: OCCT SphericalSurface.
    """
    gp_ax3 = volmdlr_frame3d_to_occt(volmdlr_surface.frame)
    return Geom_SphericalSurface(gp_ax3, volmdlr_surface.radius)


def volmdlr_toroidalsurface_to_occt(volmdlr_surface):
    """
    Create an OCCT ToroidalSurface from a voldmlr ToroidalSurface3D.

    :param volmdlr_surface: volmdlr ToroidalSurface3D.
    :return: OCCT ToroidalSurface.
    """
    gp_ax3 = volmdlr_frame3d_to_occt(volmdlr_surface.frame)
    return Geom_ToroidalSurface(gp_ax3, volmdlr_surface.major_radius, volmdlr_surface.minor_radius)


def volmdlr_bsplinesurface_to_occt(volmdlr_surface):
    """
    Create an OCCT BSplineSurface3D from a voldmlr FraBSplineSurface3Dme3D.

    :param volmdlr_surface: volmdlr BSplineSurface3D.
    :return: OCCT BSplineSurface.
    """
    points = volmdlr_surface.control_points_table
    u_deg, v_deg = volmdlr_surface.degree_u, volmdlr_surface.degree_v

    poles = TColgp_Array2OfPnt(1, len(points), 1, len(points[0]))
    for i, row in enumerate(points):
        for j, point in enumerate(row):
            poles.SetValue(i + 1, j + 1, gp_Pnt(*point))

    uknots = TColStd_Array1OfReal(1, len(volmdlr_surface.u_knots))
    for i, value in enumerate(volmdlr_surface.u_knots):
        uknots.SetValue(i + 1, value)

    vknots = TColStd_Array1OfReal(1, len(volmdlr_surface.v_knots))
    for i, value in enumerate(volmdlr_surface.v_knots):
        vknots.SetValue(i + 1, value)

    umult = TColStd_Array1OfInteger(1, len(volmdlr_surface.u_multiplicities))
    for i, value in enumerate(volmdlr_surface.u_multiplicities):
        umult.SetValue(i + 1, value)

    vmult = TColStd_Array1OfInteger(1, len(volmdlr_surface.v_multiplicities))
    for i, value in enumerate(volmdlr_surface.v_multiplicities):
        vmult.SetValue(i + 1, value)

    return Geom_BSplineSurface(poles, uknots, vknots, umult, vmult, u_deg, v_deg, False, False)
