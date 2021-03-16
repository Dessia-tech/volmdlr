#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:05:03 2019

@author: ringhausen
"""


import volmdlr as vm
import math
import matplotlib.pyplot as plt


def circle_3_segments(segment1, segment2, segment3):
    
    
    return vm.Circle2D

def circle_1_point_2_segments(point, line1, line2):
    
    # point will be called I(x_I, y_I)
    # semgent1 will be [AB]
    # segment2 will be [CD]

    I = vm.Vector2D((point[0], point[1]))
    A = vm.Vector2D((line1.points[0][0], line1.points[0][1]))
    B = vm.Vector2D((line1.points[1][0], line1.points[1][1]))
    C = vm.Vector2D((line2.points[0][0], line2.points[0][1]))
    D = vm.Vector2D((line2.points[1][0], line2.points[1][1]))
    
    # CHANGEMENT DE REPAIRE
    new_u = vm.Vector2D((B-A))
    new_u.Normalize()
    new_v = new_u.NormalVector(unit=True)
    new_basis = vm.Frame2D(I, new_u, new_v)
    
    new_A = new_basis.NewCoordinates(A)
    new_B = new_basis.NewCoordinates(B)
    new_C = new_basis.NewCoordinates(C)
    new_D = new_basis.NewCoordinates(D)
    
# =============================================================================
# LES SEGMENTS DECRIVENT UNE SEULE ET MEME DROITE
#   => AUCUNE SOLUTION   
# =============================================================================
    if new_C[1] == 0 and new_D[1] == 0:

        return None, None
    
# =============================================================================
# LES SEGMENTS SONT PARALLELES
#   => 1 SOLUTION
# =============================================================================
    elif math.isclose(line1.DirectionVector(unit=True).Dot(line2.NormalVector(unit=True)), 0, abs_tol=1e-06):

        segments_distance = abs(new_C[1] - new_A[1])
        r = segments_distance / 2
        new_circle_center = vm.Point2D((0, r))
        circle_center = new_basis.OldCoordinates(new_circle_center)
        circle = vm.Circle2D(circle_center, r)
        
        return circle, None
# =============================================================================
# LES SEGMENTS SONT PERPENDICULAIRES
#   => 2 SOLUTIONS
# =============================================================================
    elif math.isclose(line1.DirectionVector(unit=True).Dot(line2.DirectionVector(unit=True)), 0, abs_tol=1e-06):

        line_AB = vm.Line2D(vm.Point2D(new_A), vm.Point2D(new_B))
        line_CD = vm.Line2D(vm.Point2D(new_C), vm.Point2D(new_D))
        new_pt_K = vm.Point2D.LinesIntersection(line_AB ,line_CD)
        
        r = abs(new_pt_K[0])
        new_circle_center1 = vm.Point2D((0, r))
        new_circle_center2 = vm.Point2D((0, -r))
        circle_center1 = new_basis.OldCoordinates(new_circle_center1)
        circle_center2 = new_basis.OldCoordinates(new_circle_center2)
        circle1 = vm.Circle2D(circle_center1, r)
        circle2 = vm.Circle2D(circle_center2, r)
        
        return circle1, circle2
    
# =============================================================================
# LES SEGMENTS SONT QUELCONQUES
#   => 2 SOLUTIONS
# =============================================================================
    else:

        line_AB = vm.Line2D(vm.Point2D(new_A), vm.Point2D(new_B))
        line_CD = vm.Line2D(vm.Point2D(new_C), vm.Point2D(new_D))
        new_pt_K = vm.Point2D.LinesIntersection(line_AB ,line_CD)
        pt_K = vm.Point2D(new_basis.OldCoordinates(new_pt_K))

        # CHANGEMENT DE REPERE:
        new_u2 = vm.Vector2D(pt_K-I)
        new_u2.Normalize()
        new_v2 = new_u2.NormalVector(unit=True)
        new_basis2 = vm.Frame2D(I, new_u2, new_v2)
        
        new_A = new_basis2.NewCoordinates(A)
        new_B = new_basis2.NewCoordinates(B)
        new_C = new_basis2.NewCoordinates(C)
        new_D = new_basis2.NewCoordinates(D)
        new_pt_K = new_basis2.NewCoordinates(pt_K)

        teta1 = math.atan2(new_C[1], new_C[0] - new_pt_K[0])
        teta2 = math.atan2(new_D[1], new_D[0] - new_pt_K[0])
        
        if teta1 < 0:
            teta1 += math.pi
        if teta2 < 0:
            teta2 += math.pi
            
        if not math.isclose(teta1, teta2, abs_tol=1e-08):
            if math.isclose(teta1, math.pi, abs_tol=1e-08) or math.isclose(teta1, 0., abs_tol=1e-08):
                teta = teta2 
            elif math.isclose(teta2, math.pi, abs_tol=1e-08) or math.isclose(teta2, 0., abs_tol=1e-08):
                teta = teta1
        else:
            teta = teta1
            
        r1 = new_pt_K[0] * math.sin(teta) / (1 + math.cos(teta))
        r2 = new_pt_K[0] * math.sin(teta) / (1 - math.cos(teta))
        
        new_circle_center1 = vm.Point2D((0, -r1))
        new_circle_center2 = vm.Point2D((0, r2))
        
        circle_center1 = new_basis2.OldCoordinates(new_circle_center1)
        circle_center2 = new_basis2.OldCoordinates(new_circle_center2)
        
        if new_basis.NewCoordinates(circle_center1)[1] > 0:
            circle1 = vm.Circle2D(circle_center1, r1)
            circle2 = vm.Circle2D(circle_center2, r2)
        else:
            circle1 = vm.Circle2D(circle_center2, r2)
            circle2 = vm.Circle2D(circle_center1, r1)
        
        return circle1, circle2
    


def points_inside_circle(points, circle):
    point_inside_the_circle = False
    for point in points:
        if point.point_distance(circle.center) < circle.radius:
            point_inside_the_circle = True
            break
    return point_inside_the_circle


def rolling_circle_in_polygon(polygon, interpoints_distance=0.001):
    # discrétisation du polygon
    polygon_mesh = []
    # on partcourt les arrêtes
    for (vertice1, vertice2) in zip(polygon.points, polygon.points[1:]+[polygon.points[0]]):
        side_direction = vm.Vector2D((vertice2[0] - vertice1[0], vertice2[1] - vertice1[1]))
        normalized_side_direction = vm.Vector2D((vertice2[0] - vertice1[0], vertice2[1] - vertice1[1]))
        normalized_side_direction.Normalize()
        pt_number = 0
        # on ajoute les points un par un sans dépasser la longueur du côté
        segment_mesh = []
        while interpoints_distance * pt_number < side_direction.Norm():
            side_point = vertice1 + interpoints_distance * pt_number * normalized_side_direction
            segment_mesh.append(side_point)
            pt_number += 1
        polygon_mesh.append(segment_mesh)

    
    # prendre un point quelconque 
    # construire le plus grand cercle possible
    min_radius = 1e+10
    for index1, segment1 in enumerate(polygon_mesh):
        for index2, segment2 in enumerate(polygon_mesh):
            if index2-index1 >= 2 and index2-index1 < len(polygon_mesh)-1:
                
                for point1 in segment1[1:]:
                    seg1 = vm.LineSegment2D(segment1[0], segment1[-1])
                    seg2 = vm.LineSegment2D(segment2[0], segment2[-1])
                    
                    circle1, circle2  = seg2.CreateTangentCircle(point1, seg1)
                    
                    polygon_mesh_modified = polygon_mesh[:]
                    polygon_mesh_modified.pop(index1)
                    polygon_mesh_modified = [p for seg in polygon_mesh_modified for p in seg]
                    if circle1 is not None and not points_inside_circle(polygon_mesh_modified, circle1) and circle1.radius < min_radius:# and circle1.radius != 0:
                        min_radius = circle1.radius
                        min_circle = circle1
                        other_circle = circle2
                        min_point = point1
                        min_seg1 = index1
                        min_seg2 = index2
                
                for point2 in segment2[1:]:
                    seg1 = vm.LineSegment2D(segment1[0], segment1[-1])
                    seg2 = vm.LineSegment2D(segment2[0], segment2[-1])
                    
                    circle1, circle2  = seg1.CreateTangentCircle(point2, seg2)
                    
                    polygon_mesh_modified = polygon_mesh[:]
                    polygon_mesh_modified.pop(index2)
                    polygon_mesh_modified = [p for seg in polygon_mesh_modified for p in seg]
                    if circle1 is not None and not points_inside_circle(polygon_mesh_modified, circle1) and circle1.radius < min_radius:# and circle1.radius != 0:
                        min_radius = circle1.radius
                        min_circle = circle1
                        other_circle = circle2
                        min_point = point2
                        min_seg1 = index2
                        min_seg2 = index1
    
    
    return min_radius, min_circle, other_circle, min_point, min_seg1, min_seg2



# =============================================================================
# TEST
# =============================================================================

#point = vm.Point2D((0,2))
##point = vm.Point2D((1,0.5))
#
##segment1 = vm.LineSegment2D(vm.Point2D((0,0.5)), vm.Point2D((0,2)))
##segment1 = vm.LineSegment2D(vm.Point2D((0,2)), vm.Point2D((0,0.5)))
#segment1 = vm.LineSegment2D(vm.Point2D((0,2)), vm.Point2D((1,0.5)))
#
##segment2 = vm.LineSegment2D(vm.Point2D((-3,2)), vm.Point2D((-1,0)))
##segment2 = vm.LineSegment2D(vm.Point2D((-1,0)), vm.Point2D((-3,2)))
#segment2 = vm.LineSegment2D(vm.Point2D((-2,2)), vm.Point2D((-4,-1)))
#
#f = plt.figure()
##axe = plt.Axes(f, [-5, -5, 10, 10])
#axe = f.add_subplot(111)
#axe.set_aspect('equal')
#
#cercle1, cercle2 = circle_1_point_2_segments(point, segment1, segment2)
#
#cercle1.MPLPlot(axe, color='r')
#cercle2.MPLPlot(axe, color='b')
#segment1.MPLPlot(axe)
#segment2.MPLPlot(axe)
#point.MPLPlot(axe)


pts=[[-0.01,        0.23222834],
[0.01 ,      0.23222834],
[0.012  ,    0.23222834],
[0.024   ,   0.23222834],
[0.026  ,   0.2469398],
[0.07067182 ,0.2469398 ],
[0.19392024, 0.23222834],
[0.20092024, 0.23222834],
[0.20092024, 0.22722834],
[0.19392024 ,0.22722834],
[0.07067182, 0.2419398 ],
[0.026 ,    0.2419398],
[0.024 ,     0.22722834],
[0.012  ,    0.22722834],
[0.01   ,    0.22722834],
[-0.01   ,     0.22722834]]

points=[]
for pt in pts:
    points.append(vm.Point2D(tuple(pt)))
polygon = vm.Polygon2D(points)
min_radius, min_circle, other_circle, min_point, min_seg1, min_seg2 = rolling_circle_in_polygon(polygon)
print(min_radius)

fig, axe = plt.subplots()
axe.set_aspect('equal')
polygon.MPLPlot(axe)
min_seg1.MPLPlot(axe, style='r')
min_seg2.MPLPlot(axe, style='r')
min_point.MPLPlot(axe)
min_circle.MPLPlot(axe)

















