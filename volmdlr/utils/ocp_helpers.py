"""
Utils to help reading and writing OCP objects.
"""
import matplotlib.pyplot as plt
from typing import Any
import volmdlr
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Curve2d
from OCP.BRepLProp import BRepLProp_CLProps
from OCP.GCPnts import GCPnts_AbscissaPoint, GCPnts_QuasiUniformAbscissa
from OCP.GeomAbs import GeomAbs_CurveType
from OCP.BRepTools import BRepTools_WireExplorer
from OCP import TopAbs
from OCP.TopoDS import TopoDS, TopoDS_Shape

from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import (TopAbs_EDGE, TopAbs_FACE, TopAbs_VERTEX, TopAbs_WIRE, TopAbs_SHELL, TopAbs_ShapeEnum,
                        TopAbs_SOLID, TopAbs_COMPSOLID, TopAbs_COMPOUND)
from OCP.ShapeFix import ShapeFix_Shape, ShapeFix_Solid, ShapeFix_Face

from volmdlr import curves, edges, surfaces


OCCT_TO_VOLMDLR = {"Geom_SphericalSurface": surfaces.SphericalSurface3D,
                   "Geom_CylindricalSurface": surfaces.CylindricalSurface3D,
                   "Geom_Plane": surfaces.Plane3D,
                   "Geom_ToroidalSurface": surfaces.ToroidalSurface3D,
                   "Geom_ConicalSurface": surfaces.ConicalSurface3D,
                   "Geom_BSplineSurface": surfaces.BSplineSurface3D,
                   'Geom_SurfaceOfLinearExtrusion': surfaces.ExtrusionSurface3D,
                   "Geom_SurfaceOfRevolution": surfaces.RevolutionSurface3D,
                   "Geom_Line": curves.Line3D,
                   "Geom_Circle": curves.Circle3D,
                   "Geom_Ellipse": curves.Ellipse3D,
                   "Geom_BSplineCurve": edges.BSplineCurve3D,
                   "Geom_Parabola": curves.Parabola3D,
                   "Geom_Hyperbola": curves.Hyperbola3D,
                   "Geom2d_Line": curves.Line2D,
                   "Geom2d_Circle": curves.Circle2D,
                   "Geom2d_Ellipse": curves.Ellipse2D,
                   "Geom2d_BSplineCurve": edges.BSplineCurve2D
                   }


def discretize_edge(edge, number_points: int = 16, max_deflection: float = 0.01):
    """Uniformly samples an edge with specified number_points."""
    curve = BRepAdaptor_Curve(edge)
    try:
        gcpnts = GCPnts_QuasiUniformAbscissa(curve, number_points)
    except:
        return []
    curve_props = BRepLProp_CLProps(curve, 1, 1e-6)
    pts = []
    for i in range(number_points):
        point = gcpnts.Parameter(i + 1)
        curve_props.SetParameter(point)
        vpt = curve_props.Value()
        pts.append(volmdlr.Point3D(vpt.X(), vpt.Y(), vpt.Z()))
    return pts


def plot_edge(edge, ax=None, edge_style=volmdlr.core.EdgeStyle()):
    """
    Matplotlib plot method for a BSpline Curve 3D.

    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

    curve = BRepAdaptor_Curve(edge)
    if curve.GetType() == GeomAbs_CurveType.GeomAbs_Line:
        number_points = 2
    elif curve.GetType() in (GeomAbs_CurveType.GeomAbs_BezierCurve, GeomAbs_CurveType.GeomAbs_BSplineCurve):
        number_points = max(100, curve.NbPoles())
    else:
        number_points = int(72 * abs(curve.LastParameter() - curve.FirstParameter()) / volmdlr.TWO_PI)
    points = discretize_edge(edge, number_points=number_points)
    x = [point.x for point in points]
    y = [point.y for point in points]
    z = [point.z for point in points]
    ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha)
    if edge_style.edge_ends:
        ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha)
    return ax


def plot_wire(wire, ax=None, edge_style=volmdlr.core.EdgeStyle()):
    """OCP wire 3D plot using Matplotlib."""
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    exp = BRepTools_WireExplorer(wire)
    while exp.More():
        edge = exp.Current()
        plot_edge(edge, ax=ax, edge_style=edge_style)
        exp.Next()
    return ax


def plot_face(face, ax=None, color="k", alpha=1, edge_details=False):
    """Plots the face."""
    if not ax:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})

    exp = TopExp_Explorer(face, TopAbs.TopAbs_WIRE)
    while exp.More():
        wire_shape = exp.Current()
        wire = TopoDS.Wire_s(wire_shape)
        plot_wire(wire=wire,
        ax = ax, edge_style = volmdlr.core.EdgeStyle(color=color, alpha=alpha, edge_ends=edge_details,
                                                     edge_direction=edge_details)
        )
        exp.Next()
    return ax


def discretize_edge2d(edge, number_points: int = 16):
    """Discretizes a 2D edge with specified deflection."""
    curve = BRepAdaptor_Curve2d(edge)
    try:
        gcpnts = GCPnts_QuasiUniformAbscissa(curve, number_points)
    except:
        return []
    curve_props = BRepLProp_CLProps(curve, 1, 1e-6)
    pts = []
    for i in range(number_points):
        point = gcpnts.Parameter(i + 1)
        curve_props.SetParameter(point)
        vpt = curve_props.Value()
        pts.append(volmdlr.Point2D(vpt.X(), vpt.Y(), vpt.Z()))
    return pts


def plot_edge2d(edge, ax=None, edge_style=volmdlr.core.EdgeStyle()):
    """
    Matplotlib plot method for a 2D edge.

    """
    if ax is None:
        fig, ax = plt.subplots()

    # curve = BRepAdaptor_Curve2d(edge)
    # if curve.GetType() == GeomAbs_CurveType.GeomAbs_Line:
    #     number_points = 2
    # elif curve.GetType() in (GeomAbs_CurveType.GeomAbs_BezierCurve, GeomAbs_CurveType.GeomAbs_BSplineCurve):
    #     number_points = max(100, curve.NbPoles())
    # else:
    #     number_points = int(72 * abs(curve.LastParameter() - curve.FirstParameter()) / volmdlr.TWO_PI)
    points = discretize_edge2d(edge)
    x = [point.x for point in points]
    y = [point.y for point in points]
    ax.plot(x, y, color=edge_style.color, alpha=edge_style.alpha)
    if edge_style.edge_ends:
        ax.plot(x, y, color=edge_style.color, alpha=edge_style.alpha)
    return ax


def plot_wire2d(wire, ax=None, edge_style=volmdlr.core.EdgeStyle()):
    """OCP wire 2D plot using Matplotlib."""
    if ax is None:
        fig, ax = plt.subplots()
    exp = BRepTools_WireExplorer(wire)
    while exp.More():
        edge = exp.Current()
        plot_edge2d(edge, ax=ax, edge_style=edge_style)
        exp.Next()
    return ax


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
