
# import matplotlib.pyplot as plt
# import mplcyberpunk
# plt.style.use("cyberpunk")

import volmdlr
from volmdlr import edges, curves, surfaces, faces, wires, shells
from volmdlr.core import EdgeStyle


ellipse = curves.Ellipse3D(4, 3, frame=volmdlr.Frame3D(volmdlr.O3D, volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D))

ellipse_points = ellipse.discretization_points(number_points=36)
bspline_ellipse = edges.BSplineCurve3D.from_points_interpolation(ellipse_points + [ellipse_points[0]], 5)

#half ellipsoides
revolution_surface_left = surfaces.RevolutionSurface3D(bspline_ellipse, volmdlr.O3D, volmdlr.X3D)
revolution_face_left = faces.RevolutionFace3D.from_surface_rectangular_cut(revolution_surface_left, 0, volmdlr.TWO_PI, 0, 0.23*bspline_ellipse.length())
revolution_face_right = faces.RevolutionFace3D.from_surface_rectangular_cut(revolution_surface_left, 0, volmdlr.TWO_PI, 0.77*bspline_ellipse.length(), bspline_ellipse.length())
revolution_face_right = revolution_face_right.translation(volmdlr.X3D*20)

#connectors
connector1_point1 = bspline_ellipse.point_at_abscissa(0.23*bspline_ellipse.length())
connector2_point2 = bspline_ellipse.point_at_abscissa(0.77*bspline_ellipse.length())
connector2_point2 = connector2_point2.translation(volmdlr.X3D*20)
conector_radius = connector1_point1.z
conector1_center = volmdlr.Point3D(connector1_point1.x, connector1_point1.y, 0)
conector2_center = volmdlr.Point3D(connector2_point2.x, connector2_point2.y, 0)

cyl_surface = surfaces.CylindricalSurface3D(volmdlr.Frame3D(volmdlr.O3D, volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D),
                                             conector_radius)

connector1 = faces.CylindricalFace3D.from_surface_rectangular_cut(cyl_surface, 0, volmdlr.TWO_PI,
                                                                  connector1_point1.x, connector1_point1.x-1)
connector2 = faces.CylindricalFace3D.from_surface_rectangular_cut(cyl_surface, 0, volmdlr.TWO_PI,
                                                                  connector2_point2.x, connector2_point2.x+1)
closing_face_left = faces.PlaneFace3D(
    surfaces.Plane3D(volmdlr.Frame3D(conector1_center.translation(volmdlr.X3D*(-1)), volmdlr.Y3D, volmdlr.Z3D,
                                     volmdlr.X3D)),
    surfaces.Surface2D(wires.Contour2D.from_circle(curves.Circle2D(volmdlr.O2D, conector_radius)), []))
closing_face_right = faces.PlaneFace3D(
    surfaces.Plane3D(volmdlr.Frame3D(conector2_center.translation(volmdlr.X3D*1),
                                     volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D)),
    surfaces.Surface2D(wires.Contour2D.from_circle(curves.Circle2D(volmdlr.O2D, conector_radius)), []))

#cylinder body
cylindrial_surface = surfaces.CylindricalSurface3D(volmdlr.Frame3D(volmdlr.O3D, volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D), 4)
cylindrial_face = faces.CylindricalFace3D.from_surface_rectangular_cut(cylindrial_surface, 0, volmdlr.TWO_PI, 0, 20)


closedshell3d = shells.ClosedShell3D([revolution_face_left, cylindrial_face, revolution_face_right,
                                      connector1, connector2, closing_face_left, closing_face_right])
closedshell3d.babylonjs()
