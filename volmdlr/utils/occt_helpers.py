"""
Utils to help reading and writing OCP objects.
"""
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
