# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:19:52 2020

@author: Mack Pro
"""


import volmdlr as volmdlr
import volmdlr.primitives3d as p3d
# import volmdlr.primitives2d as primitives2D
import volmdlr.faces as vmf
import volmdlr.edges as vme
import volmdlr.wires as vmw
import random
import math

#### Cyl Cyl
rmin, rmax = 10, 100
mini, maxi = -0.5, 0.5

r1, r2 = random.randrange(rmin, rmax, 1)/1000, random.randrange(rmin, rmax, 1)/1000 #Choose the radius
c1 = volmdlr.Point3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center
c2 = volmdlr.Point3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center

n1 = volmdlr.Vector3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center
n2 = volmdlr.Vector3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center

n1.normalize() #Normalize the normal if it is not the case
n2.normalize()
plane1, plane2 = vmf.Plane3D.from_normal(c1, n1), vmf.Plane3D.from_normal(c2, n2) #Create a plane to give us two others vector

cylsurface1 = vmf.CylindricalSurface3D(plane1.frame, r1) 
cylsurface2 = vmf.CylindricalSurface3D(plane2.frame, r2)

hmin, hmax = 1, 100

h1, h2 = random.randrange(hmin, hmax, 5)/100, random.randrange(hmin, hmax, 5)/100 #Height of cylinder

# center2d = c1.To2D(c1, plane1.vectors[0], plane1.vectors[1])
center2d = volmdlr.Point2D(0,0)
segbh = vme.LineSegment2D(center2d, center2d + volmdlr.Point2D(0,h1)) #### Minus Pt2D because of Step adaptation
circlestart = vme.LineSegment2D(segbh.end, segbh.end+volmdlr.Point2D(2*math.pi*3/4,0)) #You can change 2*pi by an other angle
seghb = vme.LineSegment2D(circlestart.end,circlestart.end-segbh.end)
circlend = vme.LineSegment2D(seghb.end,segbh.start)
edges = [segbh, circlestart, seghb, circlend]
# points = [edges[0].start, edges[0].end]
# contours =  [vmw.Contour2D(edges)]
surface2d_1 = vmf.Surface2D(outer_contour = vmw.Contour2D(edges),
                            inner_contours = [])

# center2d2 = c2.To2D(c2, plane2.vectors[0], plane2.vectors[1])
center2d2 = volmdlr.Point2D(0,0)
segbh2 = vme.LineSegment2D(center2d2, center2d2 + volmdlr.Point2D(0,h2)) #### Minus Pt2D because of Step adaptation
circlestart2 = vme.LineSegment2D(segbh2.end, segbh2.end+volmdlr.Point2D(2*math.pi,0)) #You can change 2*pi by an other angle
seghb2 = vme.LineSegment2D(circlestart2.end,circlestart2.end-segbh2.end)
circlend2 = vme.LineSegment2D(seghb2.end,segbh2.start)
edges2 = [segbh2, circlestart2, seghb2, circlend2]
# points2 = [edges2[0].start, edges2[0].end] 
# contours2 =  [vmw.Contour2D(edges2)]
surface2d_2 = vmf.Surface2D(outer_contour = vmw.Contour2D(edges2),
                            inner_contours = [])


cyl1 = vmf.CylindricalFace3D(cylsurface1, surface2d_1)
cyl2 = vmf.CylindricalFace3D(cylsurface2, surface2d_2)

p1, p2 = cyl1.minimum_distance_points_cyl(cyl2)
print(p1.point_distance(p2))

ax = p1.plot(color='r')
cyl1.plot(ax=ax, color='r')
p2.plot(ax=ax, color='b')
cyl2.plot(ax=ax, color='b')


spheres = [p3d.Sphere(p1, 5e-3, color = (250,0,0)), 
           p3d.Sphere(p2, 5e-3, color = (0,0,250))]
vol = volmdlr.core.VolumeModel([cyl1,cyl2] + spheres)
vol.babylonjs()
