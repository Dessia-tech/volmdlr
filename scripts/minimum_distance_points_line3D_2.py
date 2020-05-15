#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:14:06 2020

@author: masfaraud
"""


import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random



radius_circle = 0.008
c = volmdlr.Circle2D(volmdlr.Point2D((0,0)), radius_circle)
contour = volmdlr.Contour2D([c])
pt0 = volmdlr.Point3D((0.01, 0.04, 0.16))
pt1 = volmdlr.Point3D((0.03, 0, 0.2))
pt2 = volmdlr.Point3D((0.45, 0.01, 0.1))
pt3 = volmdlr.Point3D((0.45, 0, -0.1))
pt4 = volmdlr.Point3D((0.3, 0.04, -0.02))
pts = [pt0, pt1, pt2, pt3, pt4]
radius = {1: 0.03, 2: 0.01, 3: 0.07}
rl = primitives3D.OpenedRoundedLineSegments3D(pts, radius, adapt_radius=True, name='wire')
sweep = primitives3D.Sweep(contour, rl, name = 'pipe')


pt10 = volmdlr.Point3D((0.02, 0.22, 0.25))
pt11 = volmdlr.Point3D((0.02, 0.24, 0.25))
pt12 = volmdlr.Point3D((0.6, 0.24, 0.20))
pt13 = volmdlr.Point3D((0.40, 0.17, 0.13))
pts1 = [pt10, pt11, pt12, pt13]
radius1 = {1: 0.01, 2: 0.05}

rl1 = primitives3D.OpenedRoundedLineSegments3D(pts1, radius1, adapt_radius=True, name='wire1')
sweep1 = primitives3D.Sweep(contour, rl1, name = 'pipe1')
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# for prim in rl.primitives :
#     prim.MPLPlot(ax=ax)
# for prim1 in rl1.primitives :
#     prim1.MPLPlot(ax=ax) 
    
l1 = rl.primitives[2]
l2 = rl1.primitives[2]

p1, p2 = l1.Matrix_distance(l2)



mes = volmdlr.Measure3D(p1, p2)
ll = primitives3D.OpenedRoundedLineSegments3D([p1, p2], {}, name='mesure')


# mes.MPLPlot(ax=ax)

model = volmdlr.VolumeModel([rl1, rl, ll])
# model.FreeCADExport('lines')



#Cas 1 

pt1 = volmdlr.Point3D((0,0,10))
pt2 = volmdlr.Point3D((0,0,4))
ptmid = ( pt1 + pt2 )/2
pt3 = volmdlr.Point3D((4,0,5))
pt4 = volmdlr.Point3D((-2,0,2))
ptmid2 = (pt3 + pt4)/2

LS1 = volmdlr.LineSegment3D(pt1, pt2)
LS2 = volmdlr.LineSegment3D(pt3, pt4)

p1, p2 = LS1.Matrix_distance(LS2)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
pt1.MPLPlot(ax=ax)
pt2.MPLPlot(ax=ax, color='r')
LS1.MPLPlot(ax=ax)

pt3.MPLPlot(ax=ax, color='g')
pt4.MPLPlot(ax=ax, color='b')
LS2.MPLPlot(ax=ax)
ptmid.MPLPlot(ax=ax)
ptmid2.MPLPlot(ax=ax)

p1.MPLPlot(ax=ax, color='m')
p2.MPLPlot(ax=ax, color='m')

d_min = LS1.minimum_distance(LS2)
# p1, p2 = LS1.Matrix_distance(LS2)
# d_min = (p1-p2).Norm()
print(d_min)
ll2 = primitives3D.OpenedRoundedLineSegments3D([p1, p2], {}, name='mesure')

model2 = volmdlr.VolumeModel([LS1, LS2, ll2])
#model2.MPLPlot()
#model2.FreeCADExport('lines2')

### Cas random
# mini, maxi = -5, 5

# pt1 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
# pt2 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
# ptmid = ( pt1 + pt2 )/2
# pt3 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
# pt4 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
# ptmid2 = (pt3 + pt4)/2

# LS1 = volmdlr.LineSegment3D(pt1, pt2)
# LS2 = volmdlr.LineSegment3D(pt3, pt4)


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# pt1.MPLPlot(ax=ax)
# pt2.MPLPlot(ax=ax, color='r')
# LS1.MPLPlot(ax=ax)

# pt3.MPLPlot(ax=ax, color='g')
# pt4.MPLPlot(ax=ax, color='b')
# LS2.MPLPlot(ax=ax)
# ptmid.MPLPlot(ax=ax)
# ptmid2.MPLPlot(ax=ax)

# d_min = LS1.minimum_distance(LS2)
# print(d_min)

### Cas LS Orthogonaux

# pt1 = volmdlr.Point3D((0,0,-3))
# pt2 = volmdlr.Point3D((4,0,1))
# ptmid = ( pt1 + pt2 )/2
# pt3 = volmdlr.Point3D((3,-1,0))
# pt4 = volmdlr.Point3D((3,2,0))
# ptmid2 = (pt3 + pt4)/2

# LS1 = volmdlr.LineSegment3D(pt1, pt2)
# LS2 = volmdlr.LineSegment3D(pt3, pt4)


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# pt1.MPLPlot(ax=ax)
# pt2.MPLPlot(ax=ax, color='r')
# LS1.MPLPlot(ax=ax)

# pt3.MPLPlot(ax=ax, color='g')
# pt4.MPLPlot(ax=ax, color='b')
# LS2.MPLPlot(ax=ax)
# ptmid.MPLPlot(ax=ax)
# ptmid2.MPLPlot(ax=ax)

# d_min = LS1.minimum_distance(LS2)
# print(d_min)

### Cas LS parallele

# pt1 = volmdlr.Point3D((2,0,5))
# pt2 = volmdlr.Point3D((2,0,0))
# ptmid = ( pt1 + pt2 )/2
# pt3 = volmdlr.Point3D((6,4,0))
# pt4 = volmdlr.Point3D((6,4,-5))
# ptmid2 = (pt3 + pt4)/2

# LS1 = volmdlr.LineSegment3D(pt1, pt2)
# LS2 = volmdlr.LineSegment3D(pt3, pt4)


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# pt1.MPLPlot(ax=ax)
# pt2.MPLPlot(ax=ax, color='r')
# LS1.MPLPlot(ax=ax)

# pt3.MPLPlot(ax=ax, color='g')
# pt4.MPLPlot(ax=ax, color='b')
# LS2.MPLPlot(ax=ax)
# ptmid.MPLPlot(ax=ax)
# ptmid2.MPLPlot(ax=ax)

# d_min = LS1.minimum_distance(LS2)
# print(d_min)

######## SWEEP TEST RANDOM

# nb_point1, nb_point2 = 6, 5

# radius_circle1, radius_circle2 = 0.008, 0.01
# c1 = volmdlr.Circle2D(volmdlr.Point2D((0,0)), radius_circle1)
# contour1 = volmdlr.Contour2D([c1])

# c2 = volmdlr.Circle2D(volmdlr.Point2D((0,0)), radius_circle2)
# contour2 = volmdlr.Contour2D([c2])

# mini, maxi = -1, 1
# pts1 = []
# for k in range (nb_point1):
#     a1, a2, a3 = random.randint(mini, maxi), random.randint(mini, maxi), random.randint(mini, maxi)
#     c1, c2, c3 = random.randrange(0,100,1), random.randrange(0,100,1), random.randrange(0,100,1)
#     pts1.append(volmdlr.Point3D((a1*c1/100, a2*c2/100, a3*c3/100)))

# radius1 = {1: 0.03, 2: 0.01, 3: 0.07, 4: 0.01}#, 5: 0.07, 6: 0.02, 7: 0.03, 8: 0.04}
# rl1 = primitives3D.OpenedRoundedLineSegments3D(pts1, radius1, adapt_radius=True, name='wire1')
# sweep1 = primitives3D.Sweep(contour1, rl1, name = 'pipe1')


# pts2 = []
# for k in range (nb_point2):
#     a1, a2, a3 = random.randint(mini, maxi), random.randint(mini, maxi), random.randint(mini, maxi)
#     c1, c2, c3 = random.randrange(0,100,1), random.randrange(0,100,1), random.randrange(0,100,1)
#     pts2.append(volmdlr.Point3D((a1*c1/100, a2*c2/100, a3*c3/100)))

# radius2 = {1: 0.01, 2: 0.05, 3: 0.06}#, 4: 0.02, 5: 0.01, 6: 0.03}
# rl2 = primitives3D.OpenedRoundedLineSegments3D(pts2, radius2, adapt_radius=True, name='wire2')
# sweep2 = primitives3D.Sweep(contour2, rl2, name = 'pipe2')

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# for prim1 in rl1.primitives :
#     prim1.MPLPlot(ax=ax)
# for prim2 in rl2.primitives :
#     prim2.MPLPlot(ax=ax)
    
# minimum_distance = rl1.minimum_distance(rl2)
# print('minimum_distance', minimum_distance)

# model=volmdlr.VolumeModel([sweep1, sweep2])

# model.babylonjs()

### Cas arc/LS

# mini, maxi = -5, 5

# pt1 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
# pt2 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
# ptmid = ( pt1 + pt2 )/2
# pt_midmid = pt1 + (pt2-pt1)/4
# pt_midmid2 = pt2 + (pt1-pt2)/4
# LS1 = volmdlr.LineSegment3D(pt1, pt2)

# pt = volmdlr.Point3D((random.randint(2*mini, 2*maxi),random.randint(2*mini, 2*maxi),random.randint(2*mini, 2*maxi)))
# radius = 2
# start, interior, end = pt, pt + volmdlr.Point3D((0,-radius,radius)),pt + volmdlr.Point3D((0,-radius,-radius))
# arc = volmdlr.Arc3D(start, interior, end)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# pt1.MPLPlot(ax=ax)
# pt2.MPLPlot(ax=ax, color='r')
# LS1.MPLPlot(ax=ax)
# start.MPLPlot(ax=ax,color='g')
# interior.MPLPlot(ax=ax,color='b')
# end.MPLPlot(ax=ax,color='y')
# arc.MPLPlot(ax=ax)
# ptmid.MPLPlot(ax=ax)

# pta1, pta2 = arc.Matrix_distance(LS1)
# pta1.MPLPlot(ax=ax, color='m')
# pta2.MPLPlot(ax=ax, color='m')

# print('int',(interior-pt2).Norm(), (interior-pt1).Norm())
# print('start',(start-pt2).Norm(), (start-pt1).Norm())
# print('end',(end-pt2).Norm(), (end-pt1).Norm())

# d_min = LS1.minimum_distance(arc)
# print('d_min',d_min)