# -*- coding: utf-8 -*-
"""

"""

import volmdlr
import volmdlr.stl as vms
import volmdlr.faces as vmf
import volmdlr.cloud
import math
import volmdlr.primitives3d as p3d

from scipy.optimize import bisect, minimize

stl = vms.Stl.from_file('C:\\Users\\Mack_Pro\\Documents\\git\\Renault\\DynamicLoop\Mise a jour GMP S30\\jeu2\\forbi_disc\\new_HR18 FDU piece chaude.stl')
shell3d = stl.to_closed_shell()

points1 = stl.extract_points()
mid_pt1 = volmdlr.O3D
for pt in points1 :
    mid_pt1 += pt
center1 = mid_pt1/(len(points1))


frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D*1.5, volmdlr.Y3D*1.5, volmdlr.Z3D*1.5)
new_shell = shell3d.frame_mapping(frame, 'old')
new_shell.alpha = 0.4
new_shell.color = (250,0,0)

center2 = center1.frame_mapping(frame, 'old')
new_shell_displaced = new_shell.translation(center1-center2)

vol = volmdlr.core.VolumeModel([shell3d, new_shell_displaced])
vol.babylonjs()

dmin = shell3d.faces[0].point1.point_distance(new_shell_displaced.faces[0].point1)
for face1, face2 in zip(shell3d.faces, new_shell_displaced.faces):
    
    p1, p2, p3 = face1.point1, face1.point2, face1.point3
    p12, p22, p32 = face2.point1, face2.point2, face2.point3
    
    if p1.point_distance(p12) < dmin :
        dmin = p1.point_distance(p12)
        
    if p2.point_distance(p22) < dmin :
        dmin = p2.point_distance(p22)
        
    if p3.point_distance(p32) < dmin :
        dmin = p3.point_distance(p32)
        
        
print(dmin)

objectif = 100e-3

# def homothetie(x) :
#     frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D*x[0], volmdlr.Y3D*x[1], volmdlr.Z3D*x[2])
#     new_shell = shell3d.frame_mapping(frame, 'old')
    
#     center2 = center1.frame_mapping(frame, 'old')
#     new_shell_displaced = new_shell.translation(center1-center2) 
    
#     dmin = shell3d.faces[0].point1.point_distance(new_shell_displaced.faces[0].point1)
#     dmin_x = abs(shell3d.faces[0].point1.dot(volmdlr.X3D)-new_shell_displaced.faces[0].point1.dot(volmdlr.X3D))
#     dmin_y = abs(shell3d.faces[0].point1.dot(volmdlr.Y3D)-new_shell_displaced.faces[0].point1.dot(volmdlr.Y3D))
#     dmin_z = abs(shell3d.faces[0].point1.dot(volmdlr.Z3D)-new_shell_displaced.faces[0].point1.dot(volmdlr.Z3D))
#     for face1, face2 in zip(shell3d.faces, new_shell_displaced.faces):
#         p1, p2, p3 = face1.point1, face1.point2, face1.point3
#         p12, p22, p32 = face2.point1, face2.point2, face2.point3
        
#         if p1.point_distance(p12) < dmin :
#             dmin = p1.point_distance(p12)
#             if abs(p1.dot(volmdlr.X3D)-p12.dot(volmdlr.X3D)) < dmin_x:
#                 dmin_x = abs(p1.dot(volmdlr.X3D)-p12.dot(volmdlr.X3D))
#             if abs(p1.dot(volmdlr.Y3D)-p12.dot(volmdlr.Y3D)) < dmin_y:
#                 dmin_y = abs(p1.dot(volmdlr.Y3D)-p12.dot(volmdlr.Y3D))
#             if abs(p1.dot(volmdlr.Z3D)-p12.dot(volmdlr.Z3D)) < dmin_z:
#                 dmin_z = abs(p1.dot(volmdlr.Z3D)-p12.dot(volmdlr.Z3D))
            
            
#         if p2.point_distance(p22) < dmin :
#             dmin = p2.point_distance(p22)
#             if abs(p2.dot(volmdlr.X3D)-p22.dot(volmdlr.X3D)) < dmin_x:
#                 dmin_x = abs(p2.dot(volmdlr.X3D)-p22.dot(volmdlr.X3D))
#             if abs(p2.dot(volmdlr.Y3D)-p22.dot(volmdlr.Y3D)) < dmin_y:
#                 dmin_y = abs(p2.dot(volmdlr.Y3D)-p22.dot(volmdlr.Y3D))
#             if abs(p2.dot(volmdlr.Z3D)-p22.dot(volmdlr.Z3D)) < dmin_z:
#                 dmin_z = abs(p2.dot(volmdlr.Z3D)-p22.dot(volmdlr.Z3D))
            
#         if p3.point_distance(p32) < dmin :
#             dmin = p3.point_distance(p32)
#             if abs(p3.dot(volmdlr.X3D)-p32.dot(volmdlr.X3D)) < dmin_x:
#                 dmin_x = abs(p3.dot(volmdlr.X3D)-p32.dot(volmdlr.X3D))
#             if abs(p3.dot(volmdlr.Y3D)-p32.dot(volmdlr.Y3D)) < dmin_y:
#                 dmin_y = abs(p3.dot(volmdlr.Y3D)-p32.dot(volmdlr.Y3D))
#             if abs(p3.dot(volmdlr.Z3D)-p32.dot(volmdlr.Z3D)) < dmin_z:
#                 dmin_z = abs(p3.dot(volmdlr.Z3D)-p32.dot(volmdlr.Z3D))
    
#         if ((p1+p2+p3)/3).point_distance((p12+p22+p32/3))< dmin :
#             dmin = ((p1+p2+p3)/3).point_distance((p12+p22+p32/3))
#             if abs(((p1+p2+p3)/3).dot(volmdlr.X3D)-((p12+p22+p32/3)).dot(volmdlr.X3D)) < dmin_x:
#                 dmin_x = abs(((p1+p2+p3)/3).dot(volmdlr.X3D)-((p12+p22+p32/3)).dot(volmdlr.X3D))
#             if abs(((p1+p2+p3)/3).dot(volmdlr.Y3D)-((p12+p22+p32/3)).dot(volmdlr.Y3D)) < dmin_y:
#                 dmin_y = abs(((p1+p2+p3)/3).dot(volmdlr.Y3D)-((p12+p22+p32/3)).dot(volmdlr.Y3D))
#             if abs(((p1+p2+p3)/3).dot(volmdlr.Z3D)-((p12+p22+p32/3)).dot(volmdlr.Z3D)) < dmin_z:
#                 dmin_z = abs(((p1+p2+p3)/3).dot(volmdlr.Z3D)-((p12+p22+p32/3)).dot(volmdlr.Z3D))
    
#     return abs(3*objectif-dmin_x-dmin_y-dmin_z)

def homothetie_x(x) :
    frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D*x[0], volmdlr.Y3D*x[0], volmdlr.Z3D*x[0])
    new_shell = shell3d.frame_mapping(frame, 'old')
    
    center2 = center1.frame_mapping(frame, 'old')
    new_shell_displaced = new_shell.translation(center1-center2) 
    
    dmin_x = abs(shell3d.faces[0].point1.dot(volmdlr.X3D)-new_shell_displaced.faces[0].point1.dot(volmdlr.X3D))
    
    for face1, face2 in zip(shell3d.faces, new_shell_displaced.faces):
        p1, p2, p3 = face1.point1, face1.point2, face1.point3
        p12, p22, p32 = face2.point1, face2.point2, face2.point3
        mid1, mid12 = (p1+p2+p3)/3, (p12+p22+p32)/3
        
        dist = []
        for pt1, pt2 in zip([p1,p2,p3,mid1], [p12,p22,p32,mid12]):
            dist.append(pt1.point_distance(pt2))
            
        if min(dist) < dmin_x :
            dmin_x = min(dist)
    
    return abs(objectif-dmin_x)

def homothetie_y(x) :
    frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D*x[0], volmdlr.Y3D*x[0], volmdlr.Z3D*x[0])
    new_shell = shell3d.frame_mapping(frame, 'old')
    
    center2 = center1.frame_mapping(frame, 'old')
    new_shell_displaced = new_shell.translation(center1-center2) 
    
    dmin_y = abs(shell3d.faces[0].point1.dot(volmdlr.Y3D)-new_shell_displaced.faces[0].point1.dot(volmdlr.Y3D))
    for face1, face2 in zip(shell3d.faces, new_shell_displaced.faces):
        p1, p2, p3 = face1.point1, face1.point2, face1.point3
        p12, p22, p32 = face2.point1, face2.point2, face2.point3
        mid1, mid12 = (p1+p2+p3)/3, (p12+p22+p32)/3
        
        dist = []
        for pt1, pt2 in zip([p1,p2,p3,mid1], [p12,p22,p32,mid12]):
            dist.append(pt1.point_distance(pt2))
            
        if min(dist) < dmin_y :
            dmin_y = min(dist)
    
    return abs(objectif-dmin_y)

def homothetie_z(x) :
    frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D*x[0], volmdlr.Y3D*x[0], volmdlr.Z3D*x[0])
    new_shell = shell3d.frame_mapping(frame, 'old')
    
    center2 = center1.frame_mapping(frame, 'old')
    new_shell_displaced = new_shell.translation(center1-center2) 
    
    dmin_z = abs(shell3d.faces[0].point1.dot(volmdlr.Z3D)-new_shell_displaced.faces[0].point1.dot(volmdlr.Z3D))
    for face1, face2 in zip(shell3d.faces, new_shell_displaced.faces):
        p1, p2, p3 = face1.point1, face1.point2, face1.point3
        p12, p22, p32 = face2.point1, face2.point2, face2.point3
        mid1, mid12 = (p1+p2+p3)/3, (p12+p22+p32)/3
        
        dist = []
        for pt1, pt2 in zip([p1,p2,p3,mid1], [p12,p22,p32,mid12]):
            dist.append(pt1.point_distance(pt2))
            
        if min(dist) < dmin_z :
            dmin_z = min(dist)
            
    return abs(objectif-dmin_z)



res_x = minimize(homothetie_x, (2),
               options={'eps': 1e-6})

res_y = minimize(homothetie_y, (2),
                options={'eps': 1e-6})

res_z = minimize(homothetie_z, (2),
                options={'eps': 1e-6})

print(res_x)
print(res_y)
print(res_z)

spheres = []
for pt in points1:
    spheres.append(p3d.Sphere(pt, objectif))
    

frameres = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D*res_x.x[0], volmdlr.Y3D*res_y.x[0], volmdlr.Z3D*res_z.x[0])
new_shell = shell3d.frame_mapping(frameres, 'old')
new_shell.alpha = 0.4
new_shell.color = (0,0,0)

center2 = center1.frame_mapping(frameres, 'old')
new_shell_displaced = new_shell.translation(center1-center2)




frameresx = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D*res_x.x[0], volmdlr.Y3D, volmdlr.Z3D)
new_shell = shell3d.frame_mapping(frameresx, 'old')
new_shell.alpha = 0.4
new_shell.color = (250,0,0)

center2 = center1.frame_mapping(frameresx, 'old')
new_shell_displaced_x = new_shell.translation(center1-center2)


frameresy = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D*res_y.x[0], volmdlr.Z3D)
new_shell = shell3d.frame_mapping(frameresy, 'old')
new_shell.alpha = 0.4
new_shell.color = (0,250,0)

center2 = center1.frame_mapping(frameresy, 'old')
new_shell_displaced_y = new_shell.translation(center1-center2)


frameresz = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D*res_z.x[0])
new_shell = shell3d.frame_mapping(frameresz, 'old')
new_shell.alpha = 0.4
new_shell.color = (0,0,250)

center2 = center1.frame_mapping(frameresz, 'old')
new_shell_displaced_z = new_shell.translation(center1-center2)

vol = volmdlr.core.VolumeModel([shell3d, new_shell_displaced, new_shell_displaced_x, new_shell_displaced_y,new_shell_displaced_z]+spheres)
vol.babylonjs()