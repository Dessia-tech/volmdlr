# -*- coding: utf-8 -*-
"""
"""

import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.primitives3d as p3d

points = [vm.Point3D(-0.51883, -0.39865, 0.5391699999999999),
          vm.Point3D(-0.5502100000000001, -0.25545000000000007, 0.5352),
          vm.Point3D(-0.42862035022394307, -0.23097502716802853, 0.45694632230872445),
           vm.Point3D(-0.3779416142154201, -0.16043402099609375, 0.5017910926712679),
           vm.Point3D(-0.38599218750000003, -0.16043402099609375, 0.5256149291992187),
           vm.Point3D(-0.3110091552734375, -0.10281935119628904, 0.5621776123046875),
           vm.Point3D(-0.2779786987304688, -0.04520467758178713, 0.6191567993164062),
           vm.Point3D(-0.27749807039051816, -0.020624128269254388, 0.5906719135439563),
           vm.Point3D(-0.18039562917956825, 0.11453519976962154, 0.6765962482987753)]#,
           # vm.Point3D(-0.1670280102373628, 0.10003019655834136, 0.6798995263156462),
           # vm.Point3D(-0.1598302, 0.10245, 0.6613972)]

bezier = vme.BezierCurve3D(degree = 3, control_points = points)

ax = bezier.plot()
for pt in points :
    pt.plot(ax=ax)
    
tangents = []
for k, pt in enumerate(bezier.points) :
    position = k/(len(bezier.points)-1)
    tangents.append(bezier.tangent(position))
    
circles = []
for pt, tan in zip(bezier.points, tangents):
    circles.append(vmw.Circle3D.from_center_normal(center = pt,
                                                   normal = tan,
                                                   radius = 10e-3))

polys = [vmw.ClosedPolygon3D(c.tessellation_points()) for c in circles]

for p in polys :
    p.plot(ax=ax, color='r')
    
points_3d = []
for poly in polys :
    points_3d.extend(poly.points)
    points_3d.append(poly.points[0])
    
size_v = len(polys[0].points)+1
size_u = len(polys)
degree_u = 3
degree_v = 3
    
bezier_surface3d = vmf.BezierSurface3D(degree_u, 
                                       degree_v,
                                       points_3d,
                                       size_u,
                                       size_v)

outer_contour = vmw.Contour2D([vme.LineSegment2D(vm.O2D, vm.X2D),
                               vme.LineSegment2D(vm.X2D, vm.X2D + vm.Y2D),
                               vme.LineSegment2D(vm.X2D + vm.Y2D, vm.Y2D),
                               vme.LineSegment2D(vm.Y2D, vm.O2D)])
surf2d = vmf.Surface2D(outer_contour, [])

bsface3d = vmf.BSplineFace3D(bezier_surface3d, surf2d)


sphere1 = p3d.Sphere(points[0], 5e-3)
sphere2 = p3d.Sphere(points[-1], 5e-3)
sphere1.color = (0,250,0)
sphere2.color = (0,250,0)

# bsface3d.babylonjs()

vol = vm.core.VolumeModel([sphere1, sphere2, bsface3d])
vol.babylonjs()

ax = bsface3d.plot()
bezier.plot(ax=ax, color='r')