"""
Utils to help reading and writing OCP objects.
"""
# pylint: disable=no-name-in-module
import matplotlib.pyplot as plt

from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.BRepLProp import BRepLProp_CLProps
from OCP.GCPnts import GCPnts_QuasiUniformAbscissa
from OCP.GeomAbs import GeomAbs_CurveType
from OCP.BRepTools import BRepTools_WireExplorer
from OCP import TopAbs
from OCP.TopoDS import TopoDS
from OCP.TopExp import TopExp_Explorer

import volmdlr
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


def discretize_edge(edge, number_points: int = 16):
    """Uniformly samples an edge with specified number_points."""
    curve = BRepAdaptor_Curve(edge)

    gcpnts = GCPnts_QuasiUniformAbscissa(curve, number_points)

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


def plot_shell(shell, ax=None, color="k", alpha=1, edge_details=False):
    """Plots the shell."""
    if not ax:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})

    exp = TopExp_Explorer(shell, TopAbs.TopAbs_FACE)
    while exp.More():
        face_shape = exp.Current()
        face = TopoDS.Face_s(face_shape)
        plot_face(face=face,
        ax = ax, color=color, alpha=alpha, edge_details=edge_details)
        exp.Next()
    return ax


def plot_solid(solid, ax=None, color="k", alpha=1, edge_details=False):
    """Plots the solid."""
    if not ax:
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})

    exp = TopExp_Explorer(solid, TopAbs.TopAbs_SHELL)
    while exp.More():
        shell_shape = exp.Current()
        shell = TopoDS.Shell_s(shell_shape)
        plot_shell(shell=shell,
        ax = ax, color=color, alpha=alpha, edge_details=edge_details)
        exp.Next()
    return ax
