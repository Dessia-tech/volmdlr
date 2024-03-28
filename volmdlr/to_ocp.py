"""
Module to translate objects in Volmdlr to OCP.
"""
# pylint: disable=no-name-in-module
from typing import List, Union
from OCP.Geom import (Geom_BSplineSurface, Geom_CylindricalSurface, Geom_ConicalSurface, Geom_ToroidalSurface,
                      Geom_SphericalSurface, Geom_SurfaceOfLinearExtrusion, Geom_Plane, Geom_BSplineCurve, Geom_Line,
                      Geom_Circle, Geom_Ellipse)
from OCP.Geom2d import Geom2d_BSplineCurve, Geom2d_Circle, Geom2d_Line, Geom2d_Ellipse
from OCP.TopoDS import TopoDS_Edge, TopoDS_Wire
from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
# from OCP.TColStd import TColStd_Array1OfReal, TColStd_Array1OfInteger
from OCP.TColgp import TColgp_Array2OfPnt, TColgp_Array1OfPnt, TColgp_Array1OfPnt2d
from OCP.gp import gp_Pnt, gp_Ax2, gp_Ax3, gp_Vec, gp_Vec2d, gp_Dir, gp_Pnt2d, gp_Dir2d, gp_Ax2d
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge2d

import volmdlr


def list_to_tcolstd_array1ofinteger(list_of_int: List[int]) -> TColStd_Array1OfInteger:
    """
    Helper function.

    """
    array = TColStd_Array1OfInteger(1, len(list_of_int))
    for i, value in enumerate(list_of_int):
        array.SetValue(i + 1, value)
    return array


def list_to_tcolstd_array1ofreal(list_of_real: List[Union[int, float]]) -> TColStd_Array1OfReal:
    """
    Helper function.

    """
    array = TColStd_Array1OfReal(1, len(list_of_real))
    for i, value in enumerate(list_of_real):
        array.SetValue(i + 1, value)
    return array


def list_to_tcolgp_array10fpnt(list_of_points: List[volmdlr.Point3D]) -> TColgp_Array1OfPnt:
    """
    Helper function.

    """
    array = TColgp_Array1OfPnt(1, len(list_of_points))
    for i, point in enumerate(list_of_points):
        array.SetValue(i + 1, point3d_to_ocp(point))
    return array


def list_to_tcolgp_array10fpnt2d(list_of_points: List[volmdlr.Point3D]) -> TColgp_Array1OfPnt2d:
    """
    Helper function.

    """
    array = TColgp_Array1OfPnt2d(1, len(list_of_points))
    for i, point in enumerate(list_of_points):
        array.SetValue(i + 1, point2d_to_ocp(point))
    return array


def point2d_to_ocp(point2d: volmdlr.Point2D) -> gp_Pnt2d:
    """
    Create an OCP Point2D from a voldmlr Point2D.

    :param point2d: volmdlr Point2D.
    :return: OCP Point2D
    """
    return gp_Pnt2d(*point2d)


def point3d_to_ocp(point: volmdlr.Point3D) -> gp_Pnt:
    """
    Create an OCP Point3D from a voldmlr Point3D.

    :param point: volmdlr Point3D.
    :return: OCP Point3D
    """
    return gp_Pnt(*point)


def vector3d_to_ocp(vector: volmdlr.Vector3D, unit_vector: bool = False) -> gp_Vec:
    """
    Create an OCP Vector3D from a voldmlr Vector3D.

    :param vector: volmdlr Vector3D.
    :param unit_vector: If set to True, returns an OCP.gp.gp_Dir object.
    :return: OCP Vector3D
    """
    if unit_vector:
        return gp_Dir(*vector)
    return gp_Vec(*vector)


def vector2d_to_ocp(vector: volmdlr.Vector2D, unit_vector: bool = False) -> gp_Vec2d:
    """
    Create an OCP Vector2D from a voldmlr Vector2D.

    :param vector: volmdlr Vector2D.
    :param unit_vector: If set to True, returns an OCP.gp.gp_Dir2d object.
    :return: OCP Vector2D
    """
    if unit_vector:
        return gp_Dir2d(*vector)
    return gp_Vec2d(*vector)


def frame3d_to_ocp(frame: volmdlr.Frame3D, right_handed: bool = False) -> Union[gp_Ax2, gp_Ax3]:
    """
    Create an OCP Frame3D from a voldmlr Frame3D.

    :param frame: volmdlr Frame3D.
    :param right_handed: (Optional) If set to True returns a gp_Ax2, that is right-handed coordinate system. If
    set to False (default) returns a gp_Ax3 can be right-handed ("direct sense") or left-handed ("indirect sense").
    :return: OCP Frame3D.
    """
    point = point3d_to_ocp(frame.origin)
    z_vector = vector3d_to_ocp(frame.w, unit_vector=True)
    x_vector = vector3d_to_ocp(frame.u, unit_vector=True)
    if right_handed:
        return gp_Ax2(point, z_vector, x_vector)
    return gp_Ax3(point, z_vector, x_vector)


def frame2d_to_ocp(frame: volmdlr.Frame2D) -> gp_Ax2d:
    """
    Create an OCP Frame2D from a voldmlr Frame2D.

    :param frame: volmdlr Frame2D.
    :return: OCP Frame2D.
    """
    point = point2d_to_ocp(frame.origin)
    x_vector = vector2d_to_ocp(frame.u, unit_vector=True)
    return gp_Ax2d(point, x_vector)


def line3d_to_ocp(line) -> Geom_Line:
    """
    Create an OCP Line from a voldmlr Line3D.

    :param line: volmdlr Line3D.
    :return: OCP Geom_Line.
    """
    point = point3d_to_ocp(line.point1)
    direction = vector3d_to_ocp(line.unit_direction_vector(), unit_vector=True)
    return Geom_Line(point, direction)


def circle3d_to_ocp(circle) -> Geom_Circle:
    """
    Create an OCP Circle from a voldmlr Circle3D.

    :param circle: volmdlr Circle3D.
    :return: OCP Geom_Circle.
    """
    frame_ax2 = frame3d_to_ocp(circle.frame, right_handed=True)
    return Geom_Circle(frame_ax2, circle.radius)


def ellipse3d_to_ocp(ellipse) -> Geom_Ellipse:
    """
    Create an OCP ellipse from a voldmlr Ellipse3D.

    :param ellipse: volmdlr Ellipse3D.
    :return: OCP Geom_Ellipse.
    """
    frame_ax2 = frame3d_to_ocp(frame=ellipse.frame, right_handed=True)
    return Geom_Ellipse(frame_ax2, ellipse.major_axis, ellipse.minor_axis)


def bsplinecurve3d_to_ocp(bsplinecurve) -> Geom_BSplineCurve:
    """
    Creates a Bspline Curve 3D from a volmdlr object.

    :param bsplinecurve: volmdlr BSpline Curve 3D.
    :return: Geom_BSplineCurve
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


def line2d_to_ocp(line) -> Geom2d_Line:
    """
    Create an OCP Line from a voldmlr Line3D.

    :param line: volmdlr Line3D.
    :return: OCP Geom_Line.
    """
    point = point2d_to_ocp(line.point1)
    direction = vector2d_to_ocp(line.unit_direction_vector(), unit_vector=True)
    return Geom2d_Line(point, direction)


def circle2d_to_ocp(circle) -> Geom2d_Circle:
    """
    Create an OCP Circle from a voldmlr Circle3D.

    :param circle: volmdlr Circle3D.
    :return: OCP Geom_Circle.
    """
    frame_ax2d = frame2d_to_ocp(circle.frame)
    return Geom2d_Circle(frame_ax2d, circle.radius, circle.is_trigo)


def ellipse2d_to_ocp(ellipse) -> Geom2d_Ellipse:
    """
    Create an OCP ellipse from a voldmlr Ellipse3D.

    :param ellipse: volmdlr Ellipse3D.
    :return: OCP Geom_Ellipse.
    """
    frame_ax2d = frame2d_to_ocp(ellipse.frame)
    return Geom2d_Ellipse(frame_ax2d, ellipse.major_axis, ellipse.minor_axis)


def bsplinecurve2d_to_ocp(bsplinecurve) -> Geom2d_BSplineCurve:
    """
    Creates a Bspline Curve 3D from a volmdlr object.

    :param bsplinecurve: volmdlr BSpline Curve 2D.
    :return: Geom2d_BSplineCurve
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


def edge2d_to_ocp(edge2d, ocp_surface=None) -> TopoDS_Edge:
    """
    Creates a OCCT edge from a volmdlr object.

    :param edge2d: volmdlr Contour3D.
    :param ocp_surface: (Optional) OCCT surface if the edge is a parametric representaion of the surface.
    :return:
    """
    volmdlr_curve = edge2d.curve()
    curve = globals()[volmdlr_curve.__class__.__name__.lower() + '_to_ocp'](volmdlr_curve)
    start = point2d_to_ocp(edge2d.start)
    end = point2d_to_ocp(edge2d.end)
    if ocp_surface:
        return BRepBuilderAPI_MakeEdge(curve, ocp_surface, volmdlr_curve.abscissa(edge2d.start),
                                       volmdlr_curve.abscissa(edge2d.end)).Edge()
    return BRepBuilderAPI_MakeEdge2d(curve, start, end).Edge()


def curve_ocp_from_edge_volmdlr(edge3d):
    """
    Gets the OCP curve from a volmdlr edge.
    """
    volmdlr_curve = edge3d.curve()
    return globals()[volmdlr_curve.__class__.__name__.lower() + '_to_ocp'](volmdlr_curve)


def edge3d_to_ocp(edge3d) -> TopoDS_Edge:
    """
    Creates a OCCT edge from a volmdlr object.

    :param edge3d: volmdlr edge 3D.
    :return:
    """
    curve = curve_ocp_from_edge_volmdlr(edge3d)
    start = point3d_to_ocp(edge3d.start)
    end = point3d_to_ocp(edge3d.end)
    return BRepBuilderAPI_MakeEdge(curve, start, end).Edge()


def contour3d_to_ocp(contour3d) -> TopoDS_Wire:
    """
    Creates a OCCT wire from a volmdlr object.

    :param contour3d: volmdlr Contour3D.
    :return:
    """
    builder = BRepBuilderAPI_MakeWire()
    for primitive in contour3d.primitives:
        builder.Add(edge3d_to_ocp(primitive))
    return builder.Wire()


def contour2d_to_ocp(contour2d, ocp_surface=None) -> TopoDS_Wire:
    """
    Creates a OCCT wire from a volmdlr object.

    :param contour2d: volmdlr Contour2D.
    :param ocp_surface: (Optional) OCCT surface if the edge is a parametric representaion of the surface.
    :return:
    """
    builder = BRepBuilderAPI_MakeWire()
    for primitive in contour2d.primitives:
        builder.Add(edge2d_to_ocp(primitive, ocp_surface=ocp_surface))

    return builder.Wire()


def plane_to_ocp(surface) -> Geom_Plane:
    """
    Create an OCP Plane from a voldmlr Plane3D.

    :param surface: volmdlr Plane3D.
    :return: OCP Plane.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_Plane(gp_ax3)


def cylindricalsurface_to_ocp(surface) -> Geom_CylindricalSurface:
    """
    Create an OCP CylindricalSurface from a voldmlr CylindricalSurface3D.

    :param surface: volmdlr CylindricalSurface3D.
    :return: OCP CylindricalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_CylindricalSurface(gp_ax3, surface.radius)


def conicalsurface_to_ocp(surface) -> Geom_ConicalSurface:
    """
    Create an OCP ConicalSurface from a voldmlr ConicalSurface3D.

    :param surface: volmdlr ConicalSurface3D.
    :return: OCP ConicalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_ConicalSurface(gp_ax3, surface.semi_angle, surface.ref_radius)


def sphericalsurface_to_ocp(surface) -> Geom_SphericalSurface:
    """
    Create an OCP SphericalSurface from a voldmlr SphericalSurface3D.

    :param surface: volmdlr SphericalSurface3D.
    :return: OCP SphericalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_SphericalSurface(gp_ax3, surface.radius)


def toroidalsurface_to_ocp(surface) -> Geom_ToroidalSurface:
    """
    Create an OCP ToroidalSurface from a voldmlr ToroidalSurface3D.

    :param surface: volmdlr ToroidalSurface3D.
    :return: OCP ToroidalSurface.
    """
    gp_ax3 = frame3d_to_ocp(surface.frame)
    return Geom_ToroidalSurface(gp_ax3, surface.major_radius, surface.minor_radius)


def bsplinesurface_to_ocp(surface) -> Geom_BSplineSurface:
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


def extrusionsurface_to_ocp(surface) -> Geom_SurfaceOfLinearExtrusion:
    """
    Create an OCP ExtrusionSurface from a voldmlr ExtrsuionSurface3D.

    :param surface: volmdlr ExtrsuionSurface3D.
    :return: OCP ExtrusionSurface.
    """
    curve = curve_ocp_from_edge_volmdlr(surface.edge)
    direction = vector3d_to_ocp(surface.direction, unit_vector=True)
    return Geom_SurfaceOfLinearExtrusion(curve, direction)


VOLMDLR_TO_OCP = {"Line3D": line3d_to_ocp,
                  "BSplineCurve3D": bsplinecurve3d_to_ocp,
                  "Circle3D": circle3d_to_ocp,
                  "Ellipse3D": ellipse3d_to_ocp,
                  "Line2D": line2d_to_ocp,
                  "BSplineCurve2D": bsplinecurve2d_to_ocp,
                  "Circle2D": circle2d_to_ocp,
                  "Ellipse2D": ellipse2d_to_ocp
                  }
