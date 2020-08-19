#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 16:14:08 2020

@author: joly
"""

import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random
import math

with_borders = False

xmax, xmin, ymax, ymin = 10, -10, 10, -10
c_min, c_max = 0.5, 6
h_min, h_max = 5, 20
extrusion_vector = vm.Z3D
x, y = vm.X3D, vm.Y3D
origin = vm.Point3D((0,0,0))
basis_plane = vm.Plane3D(origin, x, y)
alpha = 0.8

nb_step = 1
nb_components = 20
thickness = 0.5
screw_holes_diameter = 0.3
screw_holes_clearance = 0.4
n_screws = 25
thickness_min = 0.3
height_belt = thickness
minimum_stepfloor, delta_surface = 0.8, 0.8

li_m = []

class Component :
    def __init__(self, center, compo_side, vector1, vector2, height, plane) :
        self.center = center
        self.compo_side = compo_side
        self.vectors = [vector1, vector2]
        self.height = height
        self.plane = plane
        
        self.points = self.compo_points()
        self.primitives = [vm.LineSegment2D(self.points[0], self.points[1]), vm.LineSegment2D(self.points[1], self.points[2]),
                           vm.LineSegment2D(self.points[2], self.points[3]), vm.LineSegment2D(self.points[3], self.points[0])]
        self.contour = vm.Contour2D(self.primitives)
        self.solid = self.compo_solid()
        
    def compo_points(self) :
        pt1 = self.center + self.vectors[0]*self.compo_side*x_size + self.vectors[1]*self.compo_side*y_size
        pt2 = self.center + self.vectors[0]*self.compo_side*x_size - self.vectors[1]*self.compo_side*y_size
        pt3 = self.center - self.vectors[0]*self.compo_side*x_size - self.vectors[1]*self.compo_side*y_size
        pt4 = self.center - self.vectors[0]*self.compo_side*x_size + self.vectors[1]*self.compo_side*y_size
        return [pt1, pt2, pt3, pt4]
    
    def compo_solid(self) :
        extrusion_vector = self.plane.vectors[0].Cross(self.plane.vectors[1])
        return primitives3D.ExtrudedProfile(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1], 
                                            self.contour, [], extrusion_vector*self.height, color=(0.44313725, 0.27058824, 0.12156863))
        
    def update(self, new_height, new_side) :
        return Component(self.center, new_side, self.vectors[0], self.vectors[1], new_height, self.plane)
        
        self.points = self.compo_points()
        self.primitives = [vm.LineSegment2D(self.points[0], self.points[1]), vm.LineSegment2D(self.points[1], self.points[2]),
                           vm.LineSegment2D(self.points[2], self.points[3]), vm.LineSegment2D(self.points[3], self.points[0])]
        self.contour = vm.Contour2D(self.primitives)
        self.solid = self.compo_solid()
        return 
    
    def MPLPlot(self, color_center='k', color_points='k') :
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        self.center.MPLPlot(ax=ax, color=color_center)
        [pt.MPLPlot(ax=ax, color=color_points) for pt in self.points]
        [prim.MPLPlot(ax=ax) for prim in self.primitives]
        
def generate_param_component(xmin, xmax, ymin, ymax, c_min, c_max, h_min, h_max) :
    x, y = random.randrange(xmin*100, xmax*100, 1)/100, random.randrange(ymin*100, ymax*100, 1)/100
    c = random.randrange(c_min*100, c_max*100, 1)/100 
    x_vec, y_vec = random.randrange(xmin*100, xmax*100, 1)/100, random.randrange(ymin*100, ymax*100, 1)/100
    vec1 = vm.Vector2D((x_vec, y_vec))
    vec1.Normalize()
    vec2 = vec1.deterministic_unit_normal_vector()
    
    center = vm.Point2D((x, y))
    height = random.randrange(h_min*100, h_max*100, 1)/100
    return center, c, vec1, vec2, height





random_change, list_component_init = [], []
all_solid_init = []
for k in range(0, nb_components) :
    percent_side, percent_height = random.randrange(-30, 30, 1)/100, random.randrange(-30, 30, 1)/100
    random_change.append([percent_side, percent_height])
    
    x_size, y_size = random.randrange(10, 100, 1)/100, random.randrange(10, 100, 1)/100
    center, c, vec1, vec2, height = generate_param_component(xmin, xmax, ymin, ymax, c_min, c_max, h_min, h_max)
    if k==0 or nb_components%k == 0 :
        vec1, vec2 = vm.Vector2D((1,0)), vm.Vector2D((0,1))
    component = Component(center, c, vec1*x_size, vec2*y_size, height, basis_plane)
    list_component_init.append(component)
    all_solid_init.append(component.solid)

# for k in range(0, 3) :
#     list_component_init[k].MPLPlot(color_center='r', color_points='b')


for step in range(0, nb_step) :
    list_component, all_solid = [], []
    all_points, all_height = [], []
    for k in range(0, nb_components) :
        if step == 0 :
            new_component = list_component_init[k]
        else :
            component = list_component_init[k]
            new_side = component.compo_side*(1-random_change[k][0]*step/(nb_step-1))
            new_height = component.height*(1-random_change[k][1]*step/(nb_step-1))
            new_component = component.update(new_height, new_side)
            
        list_component.append(new_component)
        all_solid.append(new_component.solid)
        all_points.extend(new_component.points)
        all_height.append(new_component.height)
        
    height_sorted = sorted(all_height)
    height, max_height = height_sorted[0], height_sorted[-1]
      
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    for component in list_component :
        [pt.MPLPlot(ax=ax) for pt in component.points]
        component.center.MPLPlot(ax=ax, color='r')
        [prim.MPLPlot(ax=ax) for prim in component.primitives]
    
    poly = vm.Polygon2D.points_convex_hull(all_points)
    poly.MPLPlot(ax=ax)
    [pt.MPLPlot(ax=ax, color='m') for pt in poly.points]
    
    height_bot_belt = height
    initial_height = height_bot_belt + height_belt
    initial_area = poly.Area()
    
    list_polyfloor, list_height_polyfloor = [poly], [initial_height]
    
    for k, h in enumerate(height_sorted[1:]) :
        points_test = []
        for component in list_component :
                if component.height >= h :
                    points_test.extend(component.points)
        poly_test = vm.Polygon2D.points_convex_hull(points_test)
        
        if h > initial_height + minimum_stepfloor or k==len(height_sorted)-2:
            if poly_test.Area() < delta_surface*initial_area or k==len(height_sorted)-2: 
                list_polyfloor.append(poly_test)
                list_height_polyfloor.append(h)
                initial_height, initial_area = h, poly_test.Area()
            else :
                list_height_polyfloor[-1] = h
                if len(list_polyfloor) == 1 :
                    height_bot_belt = h
                    list_height_polyfloor[-1] += height_belt
        else :
            vol_seuil = minimum_stepfloor*delta_surface*list_polyfloor[-1].Area()
            delta_h = h-list_height_polyfloor[-1]
            vol = delta_h*(list_polyfloor[-1].Area()-poly_test.Area())
            if vol < vol_seuil and delta_h > 0.9*minimum_stepfloor:
                list_polyfloor.append(poly_test)
                list_height_polyfloor.append(h)
                initial_height, initial_area = h, poly_test.Area()
            else :
                list_height_polyfloor[-1] = h
                if len(list_polyfloor) == 1 :
                    height_bot_belt = h
                    list_height_polyfloor[-1] += height_belt
    
    while list_height_polyfloor[1]-list_height_polyfloor[0]<height_belt :
        del list_height_polyfloor[1]
        del list_polyfloor[0]
                
    # First Floor
    radius = {0: 0.1, 1: 0.1, 2: 0.1}
    for k in range(0, len(poly.points)-3) :
        radius[3+k] = 0.1 + 0.3 * random.random()
    
    contour = primitives2D.ClosedRoundedLineSegments2D(poly.points, radius,
                                                       adapt_radius=True)
    offset_radius = []
    for prim in contour.primitives :
        if prim.__class__ is vm.core.Arc2D :
            offset_radius.append(prim.radius)
            
    inner_contour = contour.Offset(-min(offset_radius))
    outer_contour = inner_contour.Offset(-thickness_min)
    
    sides = primitives3D.ExtrudedProfile(origin, x, y, outer_contour, [inner_contour],
                                         height_bot_belt * extrusion_vector, alpha=alpha, name='sides')
    
    bottom = primitives3D.ExtrudedProfile(origin, x, y, outer_contour, [],
                                          -thickness_min * extrusion_vector, alpha=alpha, name='bottom')
    
    # Screw
    screw_holes_rl = inner_contour.Offset(-(thickness+screw_holes_clearance + 0.5 * screw_holes_diameter))
    screw_holes = []
    l = screw_holes_rl.Length()
    for i in range(n_screws):
        s = i * l/n_screws
        p = screw_holes_rl.PointAtCurvilinearAbscissa(s)
        screw_holes.append(vm.Circle2D(p, screw_holes_diameter*0.5))  
    
    # Belt
    belt_outer_contour = inner_contour.Offset(-(2*screw_holes_clearance + screw_holes_diameter + thickness))
    
    belt = primitives3D.ExtrudedProfile(origin + basis_plane.normal*height_bot_belt, x, y,
                                        belt_outer_contour, [inner_contour]+screw_holes,
                                        height_belt * extrusion_vector, alpha=alpha, name='belt')
    
    # Hat floors
    list_component_hat, list_contour_floor = [], []    
    list_inner, list_outer = [], []
    primitives_floor = []
    for enum, polyfloor in enumerate(list_polyfloor[1:]) :
        primitive_of_floor = []
        radius_floor = {0: 0.1, 1: 0.1}
        for k in range(0, len(polyfloor.points)-2) :
            radius_floor[2+k] = 0.1 + 0.3 * random.random()
            
        contour_floor = primitives2D.ClosedRoundedLineSegments2D(polyfloor.points,
                                                                  radius_floor,
                                                                  adapt_radius=True)
        
        offset_radius_floor = []
        for prim in contour_floor.primitives :
            if prim.__class__ is vm.core.Arc2D :
                offset_radius_floor.append(prim.radius)
                
        inner_contour_floor = contour_floor.Offset(-min(offset_radius_floor))
        outer_contour_floor = inner_contour_floor.Offset(-thickness_min)
        
        list_inner.append(inner_contour_floor)
        list_outer.append(outer_contour_floor)
        
        height_origin = list_height_polyfloor[enum]
        
        if enum == 0 :
            belt_top = primitives3D.ExtrudedProfile(origin + basis_plane.normal*height_origin, x, y,
                                                    belt_outer_contour,
                                                    [inner_contour_floor]+screw_holes,
                                                    height_belt * basis_plane.normal, alpha=alpha, name='belt_top')
            list_component_hat.append(belt_top)
            height_origin += height_belt
            
        
        sides_floor = primitives3D.ExtrudedProfile(origin + extrusion_vector*height_origin, x, y,
                                                    outer_contour_floor, [inner_contour_floor],
                                                    (list_height_polyfloor[enum+1]-height_origin)*basis_plane.normal, alpha=alpha, name='sides_floor')
        
        list_component_hat.append(sides_floor)
        primitive_of_floor.append(sides_floor)
        
        if with_borders :
            r = 0.15
            arc = vm.Arc2D(vm.Point2D((0,0)), vm.Point2D(((-1+math.sqrt(2)/2)*r,r*math.sqrt(2)/2)), vm.Point2D((-r,r)))
            
            contour_sweep = outer_contour_floor.Offset(r)
            
            contour = vm.Contour2D([arc])
            wire_sweep = vm.Wire3D([p.To3D(origin + extrusion_vector*list_height_polyfloor[enum+1], x, y) for p in contour_sweep.primitives])
            sweep = primitives3D.Sweep(contour, wire_sweep, alpha=alpha, name = 'congé')
            
            list_component_hat.append(sweep)
            list_contour_floor.append(contour_sweep)
        else :
            r = 0
        
        if enum == len(list_polyfloor)-2: 
            if with_borders :
                list_component_hat.pop()
                r_top = thickness_min
                arc_top = vm.Arc2D(vm.Point2D((0,0)), vm.Point2D(((-1+math.sqrt(2)/2)*r_top,r_top*math.sqrt(2)/2)), vm.Point2D((-r_top,r_top)))
                contour_top = vm.Contour2D([arc_top])
                contour_sweep_top = outer_contour_floor.Offset(r_top)
            
                new_wire_sweep = vm.Wire3D([p.To3D(origin + extrusion_vector*(list_height_polyfloor[enum+1]), x, y) for p in contour_sweep_top.primitives])
                sweep = primitives3D.Sweep(contour_top, new_wire_sweep, alpha=alpha, name = 'congé')
                list_component_hat.append(sweep)
            
                top = primitives3D.ExtrudedProfile(origin + extrusion_vector*list_height_polyfloor[enum+1], x, y,
                                                   contour_sweep_top, [], 
                                                   thickness_min * basis_plane.normal, alpha=alpha, name='top')
            else :
                
                top = primitives3D.ExtrudedProfile(origin + extrusion_vector*list_height_polyfloor[enum+1], x, y,
                                                   outer_contour_floor, [], 
                                                   thickness_min * basis_plane.normal, alpha=alpha, name='top')
                
            
            list_component_hat.append(top)
            primitive_of_floor.append(top)
        primitives_floor.append(primitive_of_floor)
        
    # fig, ax = plt.subplots()
    # ax.set_aspect('equal')
    # belt_top.outer_contour2d.MPLPlot(ax=ax)
    # [pt.MPLPlot(ax=ax, color='r') for pt in belt_top.outer_contour2d.points]
    # [inner.MPLPlot(ax=ax) for inner in belt_top.inner_contours2d]
    # for inner in belt_top.inner_contours2d :
    #     [pt.MPLPlot(ax=ax, color='g') for pt in inner.points]
    
    if with_borders :
        for k, contour_floor in enumerate(list_contour_floor) :
            if k < len(list_contour_floor)-1 :
                rooftop = primitives3D.ExtrudedProfile(origin + extrusion_vector*list_height_polyfloor[k+1], x, y,
                                                        contour_floor, [list_inner[k+1]],
                                                        r*basis_plane.normal, alpha=alpha, name='rooftop')
            else :
                rooftop = primitives3D.ExtrudedProfile(origin + extrusion_vector*list_height_polyfloor[k+1], x, y,
                                                        contour_floor, [],
                                                        r*basis_plane.normal, alpha=alpha, name='rooftop')
            list_component_hat.append(rooftop)
    
    else :
        for k, outer_contour in enumerate(list_outer) :
            if k < len(list_outer)-1 :
                rooftop = primitives3D.ExtrudedProfile(origin + extrusion_vector*list_height_polyfloor[k+1], x, y,
                                                       outer_contour, [list_inner[k+1]],
                                                       thickness_min*basis_plane.normal, alpha=alpha, name='rooftop')
            else :
                rooftop = primitives3D.ExtrudedProfile(origin + extrusion_vector*list_height_polyfloor[k+1], x, y,
                                                       outer_contour, [],
                                                       thickness_min*basis_plane.normal, alpha=alpha, name='rooftop')
            list_component_hat.append(rooftop)
    
    # m = vm.VolumeModel(all_solid+[bottom])
    # m.babylonjs(debug=True)  
        
    # m = vm.VolumeModel(all_solid+[sides,bottom, belt])
    # m.babylonjs(debug=True) 
    
    # m = vm.VolumeModel(all_solid+list_component_hat)
    # m.babylonjs(debug=True)
    
    # m = vm.VolumeModel(all_solid+[sides,bottom, belt]+list_component_hat)
    # m.babylonjs()  
    primitives = all_solid+[bottom, sides, belt]+list_component_hat
    
nb_primitives = len(primitives)
away_frame = vm.Frame3D(vm.Point3D((100,100,100)), vm.X3D, vm.Y3D, vm.Z3D)
steps = [[vm.OXYZ]*nb_components+[away_frame]*(nb_primitives-nb_components),
         [vm.OXYZ]*(nb_components+1)+[away_frame]*(nb_primitives-nb_components-1), 
         [vm.OXYZ]*(nb_components+3)+[away_frame]*(nb_primitives-nb_components-3),
         [vm.OXYZ]*(nb_components+3+len(primitives_floor[0]))+[away_frame]*(nb_primitives-nb_components-3-len(primitives_floor[0])),
         [vm.OXYZ]*(nb_components+3+len(primitives_floor[0])+len(primitives_floor[1]))+[away_frame]*(nb_primitives-nb_components-3-len(primitives_floor[0])-len(primitives_floor[1])),
         [vm.OXYZ]*nb_primitives]
volmod = vm.MovingVolumeModel(primitives, steps)
volmod.babylonjs()
