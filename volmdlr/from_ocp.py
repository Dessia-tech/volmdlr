"""
Module to translate objects in OCP to Volmdlr.
"""
# pylint: disable=no-name-in-module
from typing import Any
from OCP.TColStd import TColStd_Array2OfReal, TColStd_Array1OfReal
from OCP.BRep import BRep_Tool
from OCP.BRepTools import BRepTools, BRepTools_WireExplorer
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX, TopAbs_WIRE, TopAbs_SHELL, TopAbs_ShapeEnum,
                        TopAbs_SOLID, TopAbs_COMPSOLID, TopAbs_COMPOUND)
from OCP.ShapeFix import ShapeFix_Shape
from OCP.TopoDS import (TopoDS_Face, TopoDS_Shell, TopoDS_Wire, TopoDS, TopoDS_Shape)

from OCP.TopExp import TopExp_Explorer

import volmdlr


def point2d_from_ocp(occt_point):
    """
    Instanciates a volmdlr Point2D, from occt object.

    :param occt_point: OCCT point.
    :return: volmdlr Point2D.
    """
    return volmdlr.Point2D(occt_point.X(), occt_point.Y())


def vector2d_from_ocp(occt_vector):
    """
    Instanciates a volmdlr Vector2D, from occt object.

    :param occt_vector: OCCT Vector.
    :return: volmdlr Vector2D.
    """
    return volmdlr.Vector2D(occt_vector.X(), occt_vector.Y())


def point3d_from_ocp(occt_point):
    """
    Instanciates a volmdlr Point3D, from occt object.

    :param occt_point: OCCT point.
    :return: volmdlr Point3D.
    """
    return volmdlr.Point3D(occt_point.X(), occt_point.Y(), occt_point.Z())


def vector3d_from_ocp(occt_vector):
    """
    Instanciates a volmdlr Vector3D, from occt object.

    :param occt_vector: OCCT Vector.
    :return: volmdlr Vector3D.
    """
    return volmdlr.Vector3D(occt_vector.X(), occt_vector.Y(), occt_vector.Z())


def frame3d_from_ocp_ax3(ax2):
    """
    Instanciates a volmdlr Frame3D, from occt object.

    :param ax2: OCCT ax2.
    :return: volmdlr Frame3D.
    """
    origin = point3d_from_ocp(ax2.Location())
    u = vector3d_from_ocp(ax2.XDirection())
    v = vector3d_from_ocp(ax2.YDirection())
    return volmdlr.Frame3D(origin, u, v, u.cross(v))


def frame2d_from_ocp_ax22d(ax22d):
    """
    Instanciates a volmdlr Frame2D, from occt object.

    :param ax22d: OCCT ax2 2d.
    :return: volmdlr Frame2D.
    """
    origin = point2d_from_ocp(ax22d.Location())
    u = vector2d_from_ocp(ax22d.XDirection())
    v = vector2d_from_ocp(ax22d.YDirection())
    return volmdlr.Frame2D(origin, u, v)


# Curves


def line2d_from_ocp(cls, occt_line):
    """
    Instanciates a volmdlr Line2D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_line: OCCT Line.
    :return: volmdlr Line2D.
    """
    vector = vector2d_from_ocp(occt_line.Direction())
    point = point2d_from_ocp(occt_line.Location())
    return cls.from_point_and_vector(point, vector)


def line3d_from_ocp(cls, occt_line):
    """
    Instanciates a volmdlr Line3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_line: OCCT Line 3D.
    :return: volmdlr Line3D.
    """
    position = occt_line.Position()
    return cls.from_point_and_vector(point3d_from_ocp(position.Location()),
                                     vector3d_from_ocp(position.Direction()))


def circle2d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr Circle2D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param curve: OCCT curve.
    :return: volmdlr Circle2D.
    """
    frame = frame2d_from_ocp_ax22d(curve.Circ2d().Axis())
    return cls(frame, curve.Radius())


def circle3d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr Circle3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param curve: OCCT curve.
    :return: volmdlr Circle3D.
    """
    frame = frame3d_from_ocp_ax3(curve.Position())
    return cls(frame, curve.Radius())


def ellipse2d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr Ellipse2D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param curve: OCCT curve.
    :return: volmdlr Ellipse2D.
    """
    frame = frame2d_from_ocp_ax22d(curve.Position())
    return cls(curve.MajorRadius(), curve.MinorRadius(), frame)


def ellipse3d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr Ellipse3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param curve: OCCT curve.
    :return: volmdlr Ellipse3D.
    """
    frame = frame3d_from_ocp_ax3(curve.Position())
    return cls(curve.MajorRadius(), curve.MinorRadius(), frame)


def hyperbola3d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr Hyperbola3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param curve: OCCT curve.
    :return: volmdlr Hyperbola3D.
    """
    frame = frame3d_from_ocp_ax3(curve.Position())
    return cls(frame, curve.MajorRadius(), curve.MinorRadius())


def parabola3d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr Parabola3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param curve: OCCT curve.
    :return: volmdlr Parabola3D.
    """
    frame = frame3d_from_ocp_ax3(curve.Position())
    frame = volmdlr.Frame3D(frame.origin, frame.v, frame.u, frame.w)
    return cls(frame, curve.Focal())


# Edges


def bsplinecurve2d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr BSplineCurve2D, from occt object.

    :param cls: volmdlr class to be instanciated.
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
    return cls(curve.Degree(), control_points, multiplicities, knots, weigths)


def bsplinecurve3d_from_ocp(cls, curve):
    """
    Instanciates a volmdlr BSplineCurve3D, from occt object.

    :param cls: volmdlr class to be instanciated.
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
    return cls(curve.Degree(), control_points, multiplicities, knots, weigths)


OCCT_TO_VOLMDLR = {"Geom_Line": line3d_from_ocp,
                   "Geom_Circle": circle3d_from_ocp,
                   "Geom_Ellipse": ellipse3d_from_ocp,
                   "Geom_BSplineCurve": bsplinecurve3d_from_ocp,
                   "Geom_Parabola": parabola3d_from_ocp,
                   "Geom_Hyperbola": hyperbola3d_from_ocp,
                   "Geom2d_Line": line2d_from_ocp,
                   "Geom2d_Circle": circle2d_from_ocp,
                   "Geom2d_Ellipse": ellipse2d_from_ocp,
                   "Geom2d_BSplineCurve": bsplinecurve2d_from_ocp
                   }


def volmdlr_edge_from_ocp_curve(cls, occt_curve, first, last, orientation):
    """
    Instanciate a volmdlr edge form an occt curve.

    :param occt_curve: occt curve.
    :param first: first point.
    :param last: last point.
    :param orientation: orientation of the curve to be considered.
    :return: Volmdlr trimmed edge.
    """
    start = point3d_from_ocp(occt_curve.Value(first))
    end = point3d_from_ocp(occt_curve.Value(last))
    function = OCCT_TO_VOLMDLR[occt_curve.get_type_name_s()]
    curve = function(cls, occt_curve)
    return curve.trim(start, end, same_sense=orientation == 0)


def volmdlr_edge2d_from_ocp_curve(cls, occt_curve, first, last, orientation):
    """
    Instanciate a volmdlr edge form an occt curve.

    :param occt_curve: occt curve.
    :param first: first point.
    :param last: last point.
    :param orientation: orientation of the curve to be considered.
    :return: Volmdlr trimmed edge.
    """
    start = point2d_from_ocp(occt_curve.Value(first))
    end = point2d_from_ocp(occt_curve.Value(last))
    function = OCCT_TO_VOLMDLR[occt_curve.get_type_name_s()]
    curve = function(cls, occt_curve)
    return curve.trim(start, end, same_sense=orientation == 0)


def trimmedcurve3d_from_ocp(cls, occt_curve):
    """
    Intanciate an edge from a trimmed curve in OCCT.

    :param occt_curve: occt trimmed curve.
    :return: Volmdlr trimmed edge.
    """
    start_point = point3d_from_ocp(occt_curve.StartPoint())
    end_point = point3d_from_ocp(occt_curve.EndPoint())
    occt_basis_curve = occt_curve.BasisCurve()
    volmdlr_curve = OCCT_TO_VOLMDLR[occt_basis_curve.get_type_name_s()](cls, occt_basis_curve)
    return volmdlr_curve.trim(start_point, end_point)


# Surfaces


def sphericalsurface_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr SphericalSurface3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr SphericalSurface3D.
    """
    frame = frame3d_from_ocp_ax3(occt_surface.Position())
    radius = occt_surface.Radius()
    return cls(frame, radius)


def cylindricalsurface_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr CylindricalSurface3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr CylindricalSurface3D.
    """
    frame = frame3d_from_ocp_ax3(occt_surface.Position())
    radius = occt_surface.Radius()
    return cls(frame, radius)


def plane_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr Plane3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr Plane3D.
    """
    frame = frame3d_from_ocp_ax3(occt_surface.Position())
    return cls(frame)


def toroidalsurface_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr ToroidalSurface3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr ToroidalSurface3D.
    """
    frame = frame3d_from_ocp_ax3(occt_surface.Position())
    return cls(frame, occt_surface.MajorRadius(), occt_surface.MinorRadius())


def conicalsurface_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr ConicalSurface3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr ConicalSurface3D.
    """
    frame = frame3d_from_ocp_ax3(occt_surface.Position())
    radius = occt_surface.RefRadius()
    semi_angle = occt_surface.SemiAngle()
    return cls(frame, semi_angle, radius)


def bsplinesurface_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr BSplineSurface3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr BSplineSurface3D.
    """
    array = occt_surface.Poles()
    nb_v = occt_surface.NbVPoles()
    nb_u = occt_surface.NbUPoles()

    control_points = [point3d_from_ocp(array.Value(i + 1, j + 1)) for i in range(nb_u) for j in range(nb_v)]
    u_knots = list(occt_surface.UKnots())
    u_multiplicities = list(occt_surface.UMultiplicities())
    v_knots = list(occt_surface.VKnots())
    v_multiplicities = list(occt_surface.VMultiplicities())
    weights = None
    if occt_surface.IsURational() or occt_surface.IsVRational():
        occt_surface.Weights(weights_array := TColStd_Array2OfReal(1, nb_u, 1, nb_v))

        weights = [weights_array.Value(i + 1, j + 1) for i in range(nb_u) for j in range(nb_v)]
    return cls(occt_surface.UDegree(), occt_surface.VDegree(), control_points, nb_u, nb_v,
               u_multiplicities, v_multiplicities, u_knots, v_knots, weights)


def surfaceoflinearextrusion_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr ExtrusionSurface3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr ExtrusionSurface3D.
    """
    occt_to_volmdlr_classes = kwargs.get("occt_to_volmdlr")

    occt_curve = occt_surface.BasisCurve()
    curve = OCCT_TO_VOLMDLR[occt_curve.get_type_name_s()](occt_to_volmdlr_classes[occt_curve.get_type_name_s()],
                                                          occt_curve)
    direction = vector3d_from_ocp(occt_surface.Direction())
    return cls(curve, direction)


def surfaceofrevolution_from_ocp(cls, occt_surface, **kwargs):
    """
    Instanciates a volmdlr RevolutionSurface3D, from occt object.

    :param cls: volmdlr class to be instanciated.
    :param occt_surface: OCCT surface.
    :return: volmdlr RevolutionSurface3D.
    """
    occt_to_volmdlr_classes = kwargs.get("occt_to_volmdlr")
    occt_curve = occt_surface.BasisCurve()
    curve = OCCT_TO_VOLMDLR[occt_curve.get_type_name_s()](occt_to_volmdlr_classes[occt_curve.get_type_name_s()],
                                                          occt_curve)
    axis_point = point3d_from_ocp(occt_surface.Axis().Location())
    axis_direction = vector3d_from_ocp(occt_surface.Axis().Direction())
    if curve.__class__.__name__ == "Circle3D":
        start_end = curve.center + curve.frame.u * curve.radius
        curve = curve.trim(point1=start_end, point2=start_end)
    elif curve.__class__.__name__ == "Ellipse3D":
        start_end = curve.center + curve.frame.u * curve.major_axis
        curve = curve.trim(point1=start_end, point2=start_end)
    return cls(curve, axis_point, axis_direction)


def face_from_ocp(cls, occt_face: TopoDS_Shape, occt_to_volmdlr_lut: dict, surface2d_class):
    """
    Translate an OCCT face retrieving its 2D wires into a volmdlr face.
    """
    occt_surface = BRep_Tool().Surface_s(occt_face)
    if occt_surface.get_type_name_s() == 'Geom_RectangularTrimmedSurface':
        occt_surface = occt_surface.BasisSurface()
    surface_function = globals()[occt_surface.get_type_name_s().lower()[5:] + '_from_ocp']
    surface_class = occt_to_volmdlr_lut[occt_surface.get_type_name_s()]
    surface = surface_function(surface_class, occt_surface)
    surface2d = surface2d_class.from_ocp_face(occt_face)
    return cls(surface, surface2d)


def get_wires(shape):
    """
    Get wires from a shape.

    :param OCP.TopoDS.TopoDS_Shape shape: The shape.

    :return: Wires of shape.
    :rtype: list[OCP.TopoDS.TopoDS_Wire]
    """
    if isinstance(shape, TopoDS_Wire):
        return [shape]

    exp = TopExp_Explorer(shape, TopAbs_WIRE)
    wires = []
    while exp.More():
        wire_shape = exp.Current()
        wire = TopoDS.Wire_s(wire_shape)
        wires.append(wire)
        exp.Next()
    return wires


def get_faces(shape):
    """
    Get faces from a shape.

    :param OCP.TopoDS.TopoDS_Shape shape: The shape.

    :return: Faces of shape.
    :rtype: list[OCP.TopoDS.TopoDS_Face]
    """
    if isinstance(shape, TopoDS_Face):
        return [shape]

    exp = TopExp_Explorer(shape, TopAbs_FACE)
    faces = []
    while exp.More():
        face_shape = exp.Current()
        face = TopoDS.Face_s(face_shape)
        faces.append(face)
        exp.Next()
    return faces


def get_shells(shape):
    """
    Get shells from a shape.

    :param OCP.TopoDS.TopoDS_Shape shape: The shape.

    :return: Faces of shape.
    :rtype: list[OCP.TopoDS.TopoDS_Face]
    """
    if isinstance(shape, TopoDS_Shell):
        return [shape]

    exp = TopExp_Explorer(shape, TopAbs_SHELL)
    shells = []
    while exp.More():
        shell_shape = exp.Current()
        shell = TopoDS.Shell_s(shell_shape)
        shells.append(shell)
        exp.Next()
    return shells


def get_wires_from_face(face):
    """
    Returns faces outer wire and a list of inner wires.
    """
    face_wires = get_wires(face)
    if len(face_wires) > 1:
        outer_wire = BRepTools.OuterWire_s(face)
        inner_wires = [wire for wire in face_wires if not outer_wire.IsPartner(wire)]
    else:
        outer_wire = face_wires[0]
        inner_wires = []
    return outer_wire, inner_wires


def get_contour2d_from_face_wire(contour2d_class, wire, face, occt_to_volmdlr):
    """
    Get parametric representation of the face's wires.
    """
    exp = BRepTools_WireExplorer(wire, face)
    list_edges = []
    while exp.More():
        u_start, u_end = BRep_Tool().Range_s(exp.Current(), face)
        crv = BRep_Tool().CurveOnSurface_s(exp.Current(), face, u_start, u_end, False)
        orientation = exp.Current().Orientation()
        if orientation == 1:
            u_start, u_end = u_end, u_start
        if crv.IsInstance("Geom2d_TrimmedCurve"):
            crv = crv.BasisCurve()
        list_edges.append(volmdlr_edge2d_from_ocp_curve(occt_to_volmdlr[crv.get_type_name_s()], crv,
                                                         u_start, u_end, orientation))
        exp.Next()
    if not contour2d_class(list_edges).is_ordered(1e-2):
        print("Contour not ordered")
    return contour2d_class(list_edges)


def surface2d_from_ocp_face(cls, contour2d_class, face, occt_to_volmdlr):
    """
    Builds a surface 2D (Boundary representation of a face in volmdlr) from an OCP face.
    """
    outer_wire, inner_wires = get_wires_from_face(face)
    outer_contour2d = get_contour2d_from_face_wire(contour2d_class, outer_wire, face, occt_to_volmdlr)
    inner_contours2d = [get_contour2d_from_face_wire(contour2d_class, inner_wire, face, occt_to_volmdlr)
                        for inner_wire in inner_wires]
    return cls(outer_contour2d, inner_contours2d)


downcast_LUT = {
    TopAbs_VERTEX: TopoDS.Vertex_s,
    TopAbs_EDGE: TopoDS.Edge_s,
    TopAbs_WIRE: TopoDS.Wire_s,
    TopAbs_FACE: TopoDS.Face_s,
    TopAbs_SHELL: TopoDS.Shell_s,
    TopAbs_SOLID: TopoDS.Solid_s,
    TopAbs_COMPSOLID: TopoDS.CompSolid_s,
    TopAbs_COMPOUND: TopoDS.Compound_s,
}


def shapetype(obj: TopoDS_Shape) -> TopAbs_ShapeEnum:
    """
    Returns a number from 0 to 7, representing the type of the shape.

    COMPOUND = 0
    COMPSOLID = 1
    SHELL = 2
    FACE = 3
    WIRE = 4
    EDGE = 5
    VERTEX = 6
    SHAPE = 7
    """
    if obj.IsNull():
        raise ValueError("Null TopoDS_Shape object")

    return obj.ShapeType()


def downcast(obj: TopoDS_Shape) -> TopoDS_Shape:
    """
    Downcasts a TopoDS object to suitable specialized type.
    """

    f_downcast: Any = downcast_LUT[shapetype(obj)]

    return f_downcast(obj)


def fix(obj: TopoDS_Shape) -> TopoDS_Shape:
    """
    Fix a TopoDS object to suitable specialized type.
    """

    shape_fixer = ShapeFix_Shape(obj)
    shape_fixer.Perform()

    return downcast(shape_fixer.Shape())
