# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:16:18 2020

@author: Mack Pro
"""

import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random

######## SWEEP TEST RANDOM

nb_point1, nb_point2 = 6, 5

radius_circle1, radius_circle2 = 0.008, 0.01
c1 = volmdlr.Circle2D(volmdlr.Point2D((0,0)), radius_circle1)
contour1 = volmdlr.Contour2D([c1])

c2 = volmdlr.Circle2D(volmdlr.Point2D((0,0)), radius_circle2)
contour2 = volmdlr.Contour2D([c2])

mini, maxi = -1, 1
pts1 = []
for k in range (nb_point1):
    a1, a2, a3 = random.randint(mini, maxi), random.randint(mini, maxi), random.randint(mini, maxi)
    c1, c2, c3 = random.randrange(0,100,1), random.randrange(0,100,1), random.randrange(0,100,1)
    pts1.append(volmdlr.Point3D((a1*c1/100, a2*c2/100, a3*c3/100)))

radius1 = {1: 0.03, 2: 0.01, 3: 0.07, 4: 0.01}#, 5: 0.07, 6: 0.02, 7: 0.03, 8: 0.04}
rl1 = primitives3D.OpenedRoundedLineSegments3D(pts1, radius1, adapt_radius=True, name='wire1')
sweep1 = primitives3D.Sweep(contour1, rl1, name = 'pipe1')


pts2 = []
for k in range (nb_point2):
    a1, a2, a3 = random.randint(mini, maxi), random.randint(mini, maxi), random.randint(mini, maxi)
    c1, c2, c3 = random.randrange(0,100,1), random.randrange(0,100,1), random.randrange(0,100,1)
    pts2.append(volmdlr.Point3D((a1*c1/100, a2*c2/100, a3*c3/100)))

radius2 = {1: 0.01, 2: 0.05, 3: 0.06}#, 4: 0.02, 5: 0.01, 6: 0.03}
rl2 = primitives3D.OpenedRoundedLineSegments3D(pts2, radius2, adapt_radius=True, name='wire2')
sweep2 = primitives3D.Sweep(contour2, rl2, name = 'pipe2')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for prim1 in rl1.primitives :
    prim1.MPLPlot(ax=ax)
for prim2 in rl2.primitives :
    prim2.MPLPlot(ax=ax)
    
minimum_distance = rl1.minimum_distance(rl2)
print('minimum_distance', minimum_distance)

distance_sweep = sweep1.distance_to_shell(sweep2)
print('distance_sweep', distance_sweep.distance)

pt1, pt2 = distance_sweep.points

dmin, p1 ,p2 = sweep1.faces[0].minimum_distance(sweep2.faces[0], return_points=True)
face_min1, face_min2 = sweep1.faces[0], sweep2.faces[0]
# print('dmin', dmin)
for face1 in sweep1.faces :
    for face2 in sweep2.faces :
        dtest, ptest1, ptest2 = face1.minimum_distance(face2, return_points=True)
        # print('dtest', dtest)
        if dtest < dmin :
            p1, p2 = ptest1, ptest2
            dmin = dtest
            face_min1, face_min2 = face1, face2

print('distancep1p2', (p1-p2).Norm())        
 

pt1.MPLPlot(ax=ax, color='g')
pt2.MPLPlot(ax=ax, color='g')
p1.MPLPlot(ax=ax, color='r')
p2.MPLPlot(ax=ax, color='r')

# model=volmdlr.VolumeModel([sweep1, sweep2])

# model.babylonjs()

vol = volmdlr.VolumeModel([sweep1, sweep2, p1, p2])
vol.babylonjs_from_script()

vol = volmdlr.VolumeModel([sweep1, sweep2, pt1, pt2])
vol.babylonjs_from_script()