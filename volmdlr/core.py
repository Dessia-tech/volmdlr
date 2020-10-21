#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


#import bezier

# from packaging import version
import warnings
import math
import numpy as npy
npy.seterr(divide='raise')

# from geomdl import NURBS
from geomdl import BSpline
from geomdl import utilities
import matplotlib.pyplot as plt
import  mpl_toolkits
from matplotlib.patches import Arc, FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d import proj3d
# from matplotlib import __version__ as _mpl_version

import networkx as nx

from .core_compiled import (Vector2D, Vector3D, Point2D, Point3D,
                   O2D, X2D, Y2D, OXY,
                   Basis2D, Basis3D, Frame2D, Frame3D,
                   O3D, X3D, Y3D, Z3D,
                   LineSegment2DPointDistance,
                   PolygonPointBelongs, Matrix22
                   )

from scipy.linalg import solve

import volmdlr.geometry as geometry
from volmdlr import plot_data
# from volmdlr import triangulation as tri
import triangle # doc : https://rufat.be/triangle/

import dessia_common as dc
# from typing import TypeVar, List, Tuple

from jinja2 import Environment, PackageLoader, select_autoescape

import webbrowser
import os

import tempfile
import subprocess

import random

import scipy as scp
import scipy.optimize 
# import cma


def standardize_knot_vector(knot_vector):
    u0 = knot_vector[0]
    u1 = knot_vector[-1]
    standard_u_knots = []
    if u0 != 0 or u1 != 1:
        x = 1/(u1-u0)
        y = u0/(u0-u1)
        for u in knot_vector:
            standard_u_knots.append(u*x+y)
        return standard_u_knots
    else:
        return knot_vector


def find_and_replace(string, find, replace):
    """
    Finds a string in a string and replace it
    """
    index = string.find(find)
    if index != -1:
        try:
            # verifie si il reste pas des chiffre apres
            int(string[index+len(find)])
        except (ValueError, IndexError):
            # on remplace
            return string[:index]+replace+string[index+len(find):]
        else:
            return string[:index]+find_and_replace(string[index+len(find)], find, replace)
    return string

def step_split_arguments(function_arg):
    """
    Split the arguments of a function that doesn't start with '(' but end with ')'
    ex: IN: '#123,#124,#125)'
       OUT: ['#123', '#124', '#125']
    """
    if len(function_arg) > 0 and function_arg[-1] != ')':
            function_arg += ')'
    arguments = []
    argument = ""
    parenthesis = 1
    for char in function_arg:
        if char == "(":
            parenthesis += 1

        if char != "," or parenthesis > 1:
            argument += char
        else:
            arguments.append(argument)
            argument = ""

        if char == ")":
            parenthesis -= 1
            if parenthesis == 0:
                arguments.append(argument[:-1])
                argument = ""
                break
    return arguments

def set_to_list(step_set):
    char_list = step_set.split(',')
    char_list[0] = char_list[0][1:]
    char_list[-1] = char_list[-1][:-1]
    return [elem for elem in char_list]

def delete_node_and_predecessors(graph, node):
    predecessors = list(graph.predecessors(node))
    graph.remove_node(node)
    for predecessor in predecessors:
        delete_node_and_predecessors(graph, predecessor)

def delete_node_and_successors(graph, node):
    successors = list(graph.successors(node))
    graph.remove_node(node)
    for successor in successors:
        delete_node_and_successors(graph, successor)

def clockwise_angle(vector1, vector2):
    """
    Return the clockwise angle in radians between vector1 and vector2.
    """
    # print('vector1',vector1)
    # print('vector2',vector2)
    vector0 = Vector2D((0, 0))
    if vector0 in (vector1, vector2):
        return 0
    
    dot = vector1.Dot(vector2)
    norm_vec_1 = vector1.Norm()
    norm_vec_2 = vector2.Norm()
    sol = dot/(norm_vec_1*norm_vec_2)
    cross = vector1.Cross(vector2)
    if math.isclose(sol, 1, abs_tol=1e-6) :
        inner_angle = 0
    elif math.isclose(sol, -1, abs_tol=1e-6) :
        inner_angle = math.pi
    else :
        inner_angle = math.acos(sol)
    # print('cross', cross)
    
    if cross < 0:
        return inner_angle

    return 2*math.pi-inner_angle

def vectors3d_angle(vector1, vector2):
    dot = vector1.Dot(vector2)
    theta = math.acos(dot/(vector1.Norm()*vector2.Norm()))
    # u1 = dot/(vector1.Norm()*vector2.Norm())
    # cross = vector1.Cross(vector2)
    # theta2 = math.asin(cross.Norm()/(vector1.Norm()*vector2.Norm()))
    # u2 = cross.Norm()/(vector1.Norm()*vector2.Norm())
    # theta = sin_cos_angle(u1, u2)
    
    return theta

def sin_cos_angle(u1, u2):
    #You have to give cos(theta)=u1, sin(theta)=u2
    #return an angle between 0 and 2pi
    if u1<-1 :
        u1 = -1
    elif u1>1 :
        u1 = 1
    if u2<-1 :
        u2 = -1
    elif u2>1 :
        u2 = 1
    
    if u1 > 0 :
        if u2 >= 0 :
            theta = math.acos(u1)
        else : 
            theta = 2*math.pi + math.asin(u2)
    else : 
        if u2 >= 0 :
            theta = math.acos(u1) 
        else :
            theta = 2*math.pi - math.acos(u1)
    return theta

def delete_double_pos(points, triangles):
    
    points_to_indexes = {}
    
    for index, point in enumerate(points):
        if point not in points_to_indexes:
            points_to_indexes[point] = [index]
        else:
            points_to_indexes[point].append(index)
            
    new_points = []
    index_to_modified_index = {}
    for i, (point, indexes) in enumerate(points_to_indexes.items()):
        new_points.append(point)
        index_to_modified_index[indexes[0]] = i
    
    index_to_new_index = {}
    
    for indexes in points_to_indexes.values():
        for index in indexes[1:]:
            index_to_new_index[index] = indexes[0]
    
    new_triangles = []
    # print('triangles', triangles)
    for face_triangles in triangles:
        if face_triangles is None : 
            continue
        # print('face_triangles', face_triangles)
        new_face_triangles = []
        for triangle_ in face_triangles:
            # print('triangle', triangle)
            new_triangle = []
            for index in triangle_:
                if index in index_to_new_index:
                    modified_index = index_to_modified_index[index_to_new_index[index]]
                else:
                    modified_index = index_to_modified_index[index]
                new_triangle.append(modified_index)
            new_face_triangles.append(tuple(new_triangle))
        new_triangles.append(new_face_triangles)
        
    return new_points, new_triangles

# def h_triangle(vec1, vec2) :
#     #Al-Kashi Formula
#     l1, l2 = vec1.Norm(), vec2.Norm() #two length, other than the basis
#     l3 = math.sqrt(l2**2 + l1**2 - 2*vec1.Dot(vec2))
#     pos_l3 = (l2**2 - l1**2 +l3**2)/(2*l3)
#     if pos_l3**2 > l2**2 :
#         return None
#     else :
#         h = math.sqrt(l2**2 - pos_l3**2)
#         return h

# def h2_triangle(pt1,pt2,pt3) :
#     #pt2 and pt3 form the side where we want the height
#     #Héron Formula
#     a, b, c = (pt1-pt2).Norm(), (pt1-pt3).Norm(), (pt2-pt3).Norm()
#     p = (a+b+c)/2 #half perimeter
#     area = math.sqrt(p*(p-a)*(p-b)*(p-c))
#     h = area*2/c
#     return h

def determinant(vec1, vec2, vec3) :
    a = npy.array((vec1.vector, vec2.vector, vec3.vector))
    return npy.linalg.det(a)

def delete_double_point(list_point) :
    points = []
    for pt in list_point:
        if pt not in points :
            points.append(pt)
        else :
            continue
    return points

def max_pos(list_of_float) :
    pos_max, max_float = 0, list_of_float[0]
    for pos, fl in enumerate(list_of_float):
        if pos == 0 :
            continue
        else :
            if fl > max_float :
                max_float = fl
                pos_max = pos
    return max_float, pos_max

def min_pos(list_of_float) :
    pos_min, min_float = 0, list_of_float[0]
    for pos, fl in enumerate(list_of_float):
        if pos == 0 :
            continue
        else :
            if fl < min_float :
                min_float = fl
                pos_min = pos
    return min_float, pos_min

def check_singularity(all_points):
    plus_pi, moins_pi = [], []
    for enum, pt in enumerate(all_points) :
        if pt.vector[0]>math.pi*1.01 :
            plus_pi.append(enum)
        elif pt.vector[0]<math.pi*0.99 :
            moins_pi.append(enum)
            
    if len(moins_pi) <= 2 and len(all_points)>4 :
        for pos in moins_pi :
            new_pt = all_points[pos].copy() + Point2D((2*math.pi, 0))
            if new_pt.vector[0]> 2*math.pi :
                new_pt.vector[0] = 2*math.pi
            all_points[pos] = new_pt   
    elif len(plus_pi) <=2 and len(all_points)>4 :
        for pos in plus_pi :
            new_pt = all_points[pos].copy() - Point2D((2*math.pi, 0))
            if new_pt.vector[0] < 0 :
                new_pt.vector[0] = 0
            all_points[pos] = new_pt
    if 3*len(moins_pi) <= len(plus_pi) and len(all_points)>4 :
        for pos in moins_pi :
            new_pt = all_points[pos].copy() + Point2D((2*math.pi, 0))
            if new_pt.vector[0]> 2*math.pi :
                new_pt.vector[0] = 2*math.pi
            all_points[pos] = new_pt
    elif 3*len(plus_pi) <= len(moins_pi) and len(all_points)>4 :
        for pos in plus_pi :
            new_pt = all_points[pos].copy() - Point2D((2*math.pi, 0))
            if new_pt.vector[0] < 0 :
                new_pt.vector[0] = 0
            all_points[pos] = new_pt
    
    return all_points

def posangle_arc(start, end, radius, frame=None) :
    if frame is None :
        p1_new, p2_new = start, end
    else :
        p1_new, p2_new = frame.NewCoordinates(start), frame.NewCoordinates(end)
    #Angle pour le p1
    u1, u2 = p1_new.vector[0]/radius, p1_new.vector[1]/radius
    theta1 = sin_cos_angle(u1, u2)
    #Angle pour le p2
    u3, u4 = p2_new.vector[0]/radius, p2_new.vector[1]/radius
    theta2 = sin_cos_angle(u3, u4)
    
    if math.isclose(theta1, theta2, abs_tol=1e-6) :
        if math.isclose(theta2, 0, abs_tol=1e-6) :
            theta2 += 2*math.pi
        elif math.isclose(theta1, 2*math.pi, abs_tol=1e-6) :
            theta1 -= 2*math.pi
    
    return theta1, theta2

def offset_angle(trigo, angle_start, angle_end) :
    if trigo :
        offset = angle_start
    else :
        offset = angle_end
    if angle_start > angle_end :
        angle = angle_start-angle_end
    else :
        angle = angle_end-angle_start
    return offset, angle

class Primitive2D(dc.DessiaObject):
    def __init__(self, name=''):
        self.name = name
    
        dc.DessiaObject.__init__(self, name=name)
    
    def generated_toroidalface(self, arcgen) :
        rcenter = arcgen.radius
        rcircle = self.radius
        
        center = arcgen.center
        normal = arcgen.normal
        normal.Normalize()
        
        center1 = arcgen.points[0]
        x = Vector3D((center1 - center).vector)
        x.Normalize()
        frame3d = Frame3D(center, x, normal.Cross(x), normal)
        toroidalsurface3d = ToroidalSurface3D(frame3d, rcenter*1000, rcircle*1000)

        theta = arcgen.angle
        phi = self.angle
        
        pt1, pt2, pt3, pt4 = Point2D((0, 0)), Point2D((0, phi)), Point2D((theta, phi)), Point2D((theta, 0))
        seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
        edges = [seg1, seg2, seg3, seg4]
        contours2d =  [Contour2D(edges)]
        points = [theta, phi]
        
        return ToroidalFace3D(contours2d, toroidalsurface3d, points)
    
    def generated_planeface(self, contours3D) : 
        
        return PlaneFace3D(contours3D)

    
class CompositePrimitive2D(Primitive2D):
    """
    A collection of simple primitives
    """
    def __init__(self, primitives, name=''):
        Primitive2D.__init__(self, name)
        self.primitives = primitives
        self.UpdateBasisPrimitives()

    def UpdateBasisPrimitives(self):
        basis_primitives = []
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)

        self.basis_primitives = basis_primitives


    def Rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.Rotation(center, angle, copy=True)\
                                   for p in self.primitives])
        else:
            for p in self.primitives:
                p.Rotation(center, angle, copy=False)
            self.UpdateBasisPrimitives()

    def Translation(self, offset, copy=True):
        if copy:
            return self.__class__([p.Translation(offset, copy=True)\
                                   for p in self.primitives])
        else:
            for p in self.primitives:
                p.Translation(offset, copy=False)
            self.UpdateBasisPrimitives()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return self.__class__([p.frame_mapping(frame, side, copy=True)\
                                   for p in self.primitives])
        else:
            for p in self.primitives:
                p.frame_mapping(frame, side, copy=False)
            self.UpdateBasisPrimitives()

    def To3D(self, plane_origin, x, y, name=None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return CompositePrimitive3D(primitives3D, name)

    def MPLPlot(self, ax=None, color='k', arrow=False, width=None, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        for element in self.primitives:
            if element.__class__.__name__ == 'LineSegment2D':
                element.MPLPlot(ax, color, arrow, width, plot_points=plot_points)
            else:

                element.MPLPlot(ax, color=color)

        ax.margins(0.1)
        plt.show()

        return ax

    def plot_data(self, name, fill=None, color='black', stroke_width=0.2, opacity=1):
        plot_data = {}
        plot_data['fill'] = fill
        plot_data['name'] = name
        plot_data['type'] = 'wire'
        plot_data['plot_data'] = []
        for item in self.primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                        stroke_width=stroke_width,
                                                        opacity=opacity))
        return plot_data


class Wire2D(CompositePrimitive2D):
    """
    A collection of simple primitives, following each other making a wire
    """
    def __init__(self, primitives, name=''):
        CompositePrimitive2D.__init__(self, primitives, name)

    # TODO: method to check if it is a wire

    def Length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        return ValueError

    def plot_data(self, name='', fill=None, color='black', stroke_width=1, opacity=1):
        plot_data = {}
        plot_data['name'] = name
        plot_data['type'] = 'wire'
        plot_data['plot_data'] = []
        for item in self.primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                        stroke_width=stroke_width,
                                                        opacity=opacity))
        return plot_data
    
    def line_intersection(self, line):
        intersection_points = []
        for primitive in self.primitives:
            pts = primitive.line_intersection(line)
            if pts is not None:
                if type(pts) is list:
                    intersection_points.extend(pts)
                else:
                    intersection_points.append(pts)
        if intersection_points:
            return intersection_points
        else:
            return None


class Contour2D(Wire2D):
    """
    A collection of 2D primitives forming a closed wire2D
    TODO : CenterOfMass and SecondMomentArea should be changed accordingly to
    Area considering the triangle drawn by the arcs
    """
    _non_serializable_attributes  = ['internal_arcs', 'external_arcs' ,'polygon', 'straight_line_contour_polygon']
    def __init__(self, primitives, name=''):
        Wire2D.__init__(self, primitives, name)
        self._utd_analysis = False
        self.tessel_points = self.clean_points()

    def _primitives_analysis(self):
        """
        An internal arc is an arc that has his interior point inside the polygon
        """

        arcs = []
        internal_arcs = []
        external_arcs = []
        points_polygon = []
        points_straight_line_contour = []
        for primitive in self.primitives:
            if primitive.__class__.__name__ == 'LineSegment2D':
                points_polygon.append(primitive.points[0])
                points_straight_line_contour.extend(primitive.points)
            elif primitive.__class__.__name__ == 'Arc2D':
                points_polygon.append(primitive.start)
                points_polygon.append(primitive.center)
                
                # points_polygon.append(primitive.end)
                arcs.append(primitive)
            elif primitive.__class__.__name__ == 'Circle2D':
                raise ValueError('Circle2D primitives should not be inserted in a contour, as a circle is already a contour. Use directcly the circle')
                # return None
            elif primitive.__class__.__name__ == 'OpenedRoundedLineSegments2D':
                for prim in primitive.primitives:
                    if prim.__class__.__name__ == 'LineSegment2D':
                        points_polygon.extend(prim.points)
                        points_straight_line_contour.extend(prim.points)
                    elif prim.__class__.__name__ == 'Arc2D':
        #                points_polygon.append(primitive.center)
                        points_polygon.append(prim.start)
                        points_polygon.append(prim.end)
                        arcs.append(prim)
            else:
                raise NotImplementedError('primitive of type {} is not handled'.format(primitive))

        # points_polygon = list(set(points_polygon))
        polygon = Polygon2D(points_polygon)
        points_straight_line_contour = list(set(points_straight_line_contour))
        straight_line_contour_polygon = Polygon2D(points_straight_line_contour)
        
        for arc in arcs:
            if polygon.PointBelongs(arc.interior):
                internal_arcs.append(arc)
            else:
                external_arcs.append(arc)

        return internal_arcs, external_arcs, polygon, straight_line_contour_polygon
    
    def _get_internal_arcs(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon, self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._internal_arcs
            
    internal_arcs = property(_get_internal_arcs)
    
    def _get_external_arcs(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon, self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._external_arcs
            
    external_arcs = property(_get_external_arcs)
    
    def _get_polygon(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon, self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._polygon
            
    polygon = property(_get_polygon)
    
    def _get_straight_line_contour_polygon(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon, self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._straight_line_contour_polygon
            
    straight_line_contour_polygon = property(_get_straight_line_contour_polygon)
    
    
    def point_belongs(self, point):
        for arc in self.internal_arcs:
            if arc.point_belongs(point):
                return False
        if self.polygon.PointBelongs(point):
            return True
        for arc in self.external_arcs:
            if arc.point_belongs(point):
                return True
        return False

    def point_distance(self, point):
        min_distance = self.primitives[0].point_distance(point)
        for primitive in self.primitives[1:]:
            distance = primitive.point_distance(point)
            if distance < min_distance:
                min_distance = distance
        return min_distance

    def bounding_points(self):
        points = self.straight_line_contour_polygon.points[:]
        for arc in self.internal_arcs + self.external_arcs:
            points.extend(arc.tessellation_points())
        xmin = min([p[0] for p in points])
        xmax = max([p[0] for p in points])
        ymin = min([p[1] for p in points])
        ymax = max([p[1] for p in points])
        return (Point2D((xmin, ymin)), Point2D((xmax, ymax)))

    def To3D(self, plane_origin, x, y, name=None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return Contour3D(edges=primitives3D, name=name)

    def Area(self):
        if len(self.primitives) == 1:
            return self.primitives[0].Area()

        A = self.polygon.Area()

        for arc in self.internal_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A - arc.Area() + triangle.Area()
        for arc in self.external_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A + arc.Area() - triangle.Area()

        return A

    def CenterOfMass(self):
        if len(self.primitives) == 1:
            return self.primitives[0].CenterOfMass()

        area = self.polygon.Area()
        if area > 0.:
            c = area*self.polygon.CenterOfMass()
        else:
            c = O2D

        for arc in self.internal_arcs:
            arc_area = arc.Area()
            c -= arc_area*arc.CenterOfMass()
            area -= arc_area
        for arc in self.external_arcs:
            arc_area = arc.Area()
            c += arc_area*arc.CenterOfMass()
            area += arc_area
        if area != 0:
            return c/area
        else:
            return False

    def SecondMomentArea(self, point):
        if len(self.primitives) == 1:
            return self.primitives[0].SecondMomentArea(point)

        A = self.polygon.SecondMomentArea(point)
        for arc in self.internal_arcs:
            A -= arc.SecondMomentArea(point)
        for arc in self.external_arcs:
            A += arc.SecondMomentArea(point)

        return A
    
    def plot_data(self, name='', fill=None, marker=None, color='black', 
                  stroke_width=1, dash=False, opacity=1):

        plot_data = {}
        plot_data['fill'] = fill
        plot_data['name'] = name
        plot_data['type'] = 'contour'
        plot_data['plot_data'] = []
        for item in self.primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                         stroke_width=stroke_width,
                                                         opacity=opacity))
        return plot_data

    def copy(self) :
        primitives_copy = []
        for primitive in self.primitives :
            primitives_copy.append(primitive.copy())
        return Contour2D(primitives_copy)

    def average_center_point(self):
        nb = len(self.tessel_points )
        x = npy.sum([p[0] for p in self.tessel_points]) / nb
        y = npy.sum([p[1] for p in self.tessel_points]) / nb
        return Point2D((x,y))
    
    def clean_points(self):
        """
        This method is copy from Contour3D, if changes are done there or here,
        please change both method
        Be aware about primitives = 2D, edges = 3D
        """
        if hasattr(self.primitives[0], 'endpoints'):
            points = self.primitives[0].endpoints[:]
        else:
            points = self.primitives[0].tessellation_points()
        for primitive in self.primitives[1:]:
            if hasattr(primitive, 'endpoints'):
                points_to_add = primitive.endpoints[:]
            else :
                points_to_add = primitive.tessellation_points()
            if points[0] == points[-1]: # Dans le cas où le (dernier) edge relie deux fois le même point
                points.extend(points_to_add[::-1])
            
            elif points_to_add[0] == points[-1]:
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[-1]:
                points.extend(points_to_add[-2::-1])
            elif points_to_add[0] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[-2::-1])
            else: 
                d1, d2 = (points_to_add[0]-points[0]).Norm(), (points_to_add[0]-points[-1]).Norm()
                d3, d4 = (points_to_add[-1]-points[0]).Norm(), (points_to_add[-1]-points[-1]).Norm()
                if math.isclose(d2, 0, abs_tol=1e-3):
                    points.extend(points_to_add[1:])
                elif math.isclose(d4, 0, abs_tol=1e-3):
                    points.extend(points_to_add[-2::-1])
                elif math.isclose(d1, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[1:])
                elif math.isclose(d3, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[-2::-1])
                    
        if len(points) > 1:
            if points[0] == points[-1]:
                points.pop()
        return points
    

class Mesh2D:
    def __init__(self, contours, points_densities, default_density):
        self.contours = contours
        self.points_densities = points_densities
        self.default_density = default_density

    def GeoScript(self, filepath=''):
        s = ''
        ipt = 1# point index
        ipr = 1# primitive index
        points_index = {}
        #assigning an index to point
        for contour in self.contours:
            for primitive in contour.primitives:
                for point in primitive.geo_points:
                    try:
                        points_index[point]
                    except KeyError:
                        points_index[point]=ipt
                        try:
                            d=self.points_densities[point]
                        except KeyError:
                            d=self.default_density
                        s += 'Point({})={{{},{},0.,{}}};\n'.format(ipt,*point.vector,d)
                        ipt += 1
        contours_indices = []
        for contour in self.contours:
            contour_iprs = []
            for primitive in contour.primitives:
                spr,ipr2=primitive.GeoScript(ipr,[points_index[p] for p in primitive.geo_points])
                s+=spr
                contour_iprs.extend(range(ipr,ipr2))
                ipr=ipr2
            s+='Line Loop({}) = {{{}}};\n'.format(ipr,str(contour_iprs)[1:-1])
            contours_indices.append(ipr)
            ipr+=1
        s+='Plane Surface({}) = {{{}}};\n'.format(ipr,str(contours_indices)[1:-1])
        # Saving to file if required
        if filepath!='':
            with open(filepath, 'w') as file:
                file.write(s)
        return s

class Line(dc.DessiaObject):
    def __neg__(self):
        return self.__class__(self.points[::-1])

    def DirectionVector(self, unit=False):
        u = self.points[1] - self.points[0]
        if unit:
            u.Normalize()
        return u

    def NormalVector(self, unit=False):
        return self.DirectionVector(unit).NormalVector()

class Line2D(Primitive2D, Line):
    """
    Define an infinite line given by two points.
    """
    def __init__(self, point1, point2,*, name=''):
        Primitive2D.__init__(self, name=name)
        self.points=[point1, point2]
        
        self.endpoints=[point1, point2]
        
        self.point1 = point1
        self.point2 = point2

    def To3D(self, plane_origin, x1, x2):
        p3D = [p.To3D(plane_origin, x1, x2) for p in self.points]
        return Line2D(*p3D, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Line2D(*[p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Line2D(*[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def point_distance(self, point, return_other_point=False):
        """
        Computes the distance of a point to line
        """
        p1, p2 = self.points
        u = p2 - p1
        t = (point-p1).Dot(u) / u.Norm()**2
        projection = p1 + t * u # Projection falls on the segment
        if return_other_point:
            return (point-projection).Norm(), projection
        return (point-projection).Norm()

    def PointProjection(self, point, curvilinear_abscissa=False):
        p1, p2 = self.points
        u = p2 - p1
        t = (point-p1).Dot(u) / u.Norm()**2
        projection = p1 + t * u
        if curvilinear_abscissa:
            return projection,t
        return projection

    def MPLPlot(self, ax=None, color='k', dashed=True):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = ax.figure
        
        p1, p2 = self.points

        # Axline disappeared in matplotlib 3.2.0 but was in 3.2.0.rc ...
        # TODO: if it comes back implement it...

        # if version.parse(_mpl_version) >= version.parse('3.2'):
        #     if dashed:
        #         ax.axline(*p1, *p2, dashes=[30, 5, 10, 5])
        #     else:
        #         ax.axline(*p1, *p2)
        # else:
        u = p2 - p1
        p3 = p1 - 3*u
        p4 = p2 + 4*u
        if dashed:
            ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color, dashes=[30, 5, 10, 5])
        else:
            ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color)


        return ax
    
    def plot_data(self, marker=None, color='black', stroke_width=1,
                 dash=False, opacity=1, arrow=False):
        p1, p2 = self.points
        u = p2 - p1
        p3 = p1 - 3*u
        p4 = p2 + 4*u
        return {'type' : 'line',
                'data' : [p3[0], p3[1],
                          p4[0], p4[1]],
                'color' : color,
                'marker' : marker,
                'size' : stroke_width,
                'dash' : dash,
                'opacity' : opacity,
                'arrow': arrow
                }

    def CreateTangentCircle(self, point, other_line):
        """
        Computes the two circles that are tangent to 2 lines and intersect
        a point located on one of the two lines.
        """

        # point will be called I(x_I, y_I)
        # self will be (AB)
        # line will be (CD)

        if math.isclose(self.point_distance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            A = Vector2D((self.points[0][0], self.points[0][1]))
            B = Vector2D((self.points[1][0], self.points[1][1]))
            C = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            D = Vector2D((other_line.points[1][0], other_line.points[1][1]))

        elif math.isclose(other_line.point_distance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            C = Vector2D((self.points[0][0], self.points[0][1]))
            D = Vector2D((self.points[1][0], self.points[1][1]))
            A = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            B = Vector2D((other_line.points[1][0], other_line.points[1][1]))
        else:
            raise AttributeError("The point isn't on any of the two lines")

        # CHANGEMENT DE REPAIRE
        new_u = Vector2D((B-A))
        new_u.Normalize()
        new_v = new_u.NormalVector(unit=True)
        new_basis = Frame2D(I, new_u, new_v)

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
        elif math.isclose(self.DirectionVector(unit=True).Dot(other_line.NormalVector(unit=True)), 0, abs_tol=1e-06):

            segments_distance = abs(new_C[1] - new_A[1])
            r = segments_distance / 2
            new_circle_center = Point2D((0, npy.sign(new_C[1] - new_A[1])*r))
            circle_center = new_basis.OldCoordinates(new_circle_center)
            circle = Circle2D(circle_center, r)

            return circle, None

# =============================================================================
# LES SEGMENTS SONT PERPENDICULAIRES
#   => 2 SOLUTIONS
# =============================================================================
        elif math.isclose(self.DirectionVector(unit=True).Dot(other_line.DirectionVector(unit=True)), 0, abs_tol=1e-06):

            line_AB = Line2D(Point2D(new_A), Point2D(new_B))
            line_CD = Line2D(Point2D(new_C), Point2D(new_D))
            new_pt_K = Point2D.LinesIntersection(line_AB ,line_CD)

            r = abs(new_pt_K[0])
            new_circle_center1 = Point2D((0, r))
            new_circle_center2 = Point2D((0, -r))
            circle_center1 = new_basis.OldCoordinates(new_circle_center1)
            circle_center2 = new_basis.OldCoordinates(new_circle_center2)
            circle1 = Circle2D(circle_center1, r)
            circle2 = Circle2D(circle_center2, r)

            return circle1, circle2

# =============================================================================
# LES SEGMENTS SONT QUELCONQUES
#   => 2 SOLUTIONS
# =============================================================================
        else:

            line_AB = Line2D(Point2D(new_A), Point2D(new_B))
            line_CD = Line2D(Point2D(new_C), Point2D(new_D))
            new_pt_K = Point2D.LinesIntersection(line_AB ,line_CD)
            pt_K = Point2D(new_basis.OldCoordinates(new_pt_K))

            if pt_K == I:
                return None, None

            # CHANGEMENT DE REPERE:
            new_u2 = Vector2D(pt_K-I)
            new_u2.Normalize()
            new_v2 = new_u2.NormalVector(unit=True)
            new_basis2 = Frame2D(I, new_u2, new_v2)

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

            new_circle_center1 = Point2D((0, -r1))
            new_circle_center2 = Point2D((0, r2))

            circle_center1 = new_basis2.OldCoordinates(new_circle_center1)
            circle_center2 = new_basis2.OldCoordinates(new_circle_center2)

            if new_basis.NewCoordinates(circle_center1)[1] > 0:
                circle1 = Circle2D(circle_center1, r1)
                circle2 = Circle2D(circle_center2, r2)
            else:
                circle1 = Circle2D(circle_center2, r2)
                circle2 = Circle2D(circle_center1, r1)

            return circle1, circle2

class BSplineCurve2D(Primitive2D) :
    def __init__(self, degree, control_points, knot_multiplicities, knots, weights=None, periodic=False, name=''):
        Primitive2D.__init__(self, name=name)
        self.control_points = control_points
        self.degree = degree
        knots = standardize_knot_vector(knots)
        self.knots = knots
        self.knot_multiplicities = knot_multiplicities
        self.weights = weights
        self.periodic = periodic
        self.name = name

        curve = BSpline.Curve()
        curve.degree = degree
        if weights is None:
            P = [(control_points[i][0], control_points[i][1]) for i in range(len(control_points))]
            curve.ctrlpts = P
        else:
            Pw = [(control_points[i][0]*weights[i], control_points[i][1]*weights[i], weights[i]) for i in range(len(control_points))]
            curve.ctrlptsw = Pw
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot]*knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.1
        curve_points = curve.evalpts

        self.curve = curve
        self.points = [Point2D((p[0], p[1])) for p in curve_points]

    def Length(self):
        #Approximately
        length = 0
        for k in range(0,len(self.points)-1) :
            length += (self.points[k] - self.points[k+1]).Norm()
        return length
    
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        #copy paste from wire3D
        length = 0.
        primitives = []
        for k in range(0, len(self.points)-1) :
            primitives.append(LineSegment2D(self.points[k], self.points[k+1]))
        for primitive in primitives:
            primitive_length = primitive.Length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError
    
    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        x=[p.vector[0] for p in self.points]
        y=[p.vector[1] for p in self.points]
        ax.plot(x,y, 'o-k')
        return fig, ax

    def To3D(self, plane_origin, x1, x2):
        control_points3D=[p.To3D(plane_origin,x1,x2) for p in self.control_points]
        return BSplineCurve3D(self.degree, control_points3D, self.knot_multiplicities, self.knots, self.weights, self.periodic, self.name)
    
    def tessellation_points(self) :
        return self.points

class LineSegment2D(Line2D):
    """
    Define a line segment limited by two points
    """
    def __init__(self,point1, point2, *,name=''):
        Line2D.__init__(self, point1, point2, name = name)
       
        
    def to_dict(self):
        # improve the object structure ?
        dict_ = {}
        dict_['name'] = self.name
        dict_['point1'] = self.points[0].to_dict()
        dict_['point2'] = self.points[1].to_dict()
        dict_['object_class'] = 'volmdlr.core.LineSegment2D'
        return dict_
    
    @classmethod
    def dict_to_object(cls, dict_):
        return cls(point1 = Point2D.dict_to_object(dict_['point1']), 
                   point2 = Point2D.dict_to_object(dict_['point2']), name = dict_['name'])

    def _get_geo_points(self):
        return self.points

    geo_points=property(_get_geo_points)

    def Length(self):
        return self.points[1].point_distance(self.points[0])

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + self.DirectionVector(unit=True) * curvilinear_abscissa

    def point_distance(self, point, return_other_point=False):
        """
        Computes the distance of a point to segment of line
        """
        if self.points[0] == self.points[1]:
            if return_other_point:
                return 0, Point2D(point)
            return 0
        distance, point = LineSegment2DPointDistance([p.vector for p in self.points], point.vector)
        if return_other_point:
            return distance, Point2D(point)
        return distance

    def PointProjection(self, point, curvilinear_abscissa=False):
        point, curv_abs = Line2D.PointProjection(self, point, True)
        if curv_abs <= 0.:
            point = self.points[0]
            curv_abs = 0.
        elif curv_abs >= 1.:
            point = self.points[1]
            curv_abs = 1.

        if curvilinear_abscissa:
            return point, curv_abs
        else:
            return point
        
    def PointProjection2(self, point, curvilinear_abscissa=False):
        """
        If the projection falls outside the LineSegment2D, returns None.
        """
        point, curv_abs = Line2D.PointProjection(self, point, True)
        if curv_abs < 0 or curv_abs > 1:
            if curvilinear_abscissa:
                return None, curv_abs
            else:
                return None
        if curvilinear_abscissa:
            return point, curv_abs
        else:
            return point
        
    # def line_intersection(self, line):
    #     point = Point2D.LinesIntersection(self, line)
    #     if point is not None:
    #         point_projection = self.PointProjection2(point)
    #         if line.__class__ is LineSegment2D:
    #             point_projection2 = line.PointProjection2(point)
    #             if point_projection is None or point_projection2 is None:
    #                 return None
    #         return point_projection
    #     else:
    #         return None
        
    def line_intersection(self, line):
        point = Point2D.LinesIntersection(self, line)
        if point is not None:
            point_projection1 = self.PointProjection2(point)
            if point_projection1 is None:
                return None
            
            if line.__class__.__name__ == 'LineSegment2D':
                point_projection2 = line.PointProjection2(point)
                if point_projection2 is None:
                    return None
                
            return point_projection1
        else:
            return None

    def MPLPlot(self, ax=None, color='k', arrow=False, width=None, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            # ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        p1, p2 = self.points
        if arrow:
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, style='o-')
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color)
                
            length = ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5
            if width is None:
                width = length / 1000.
                head_length = length/20.
                head_width = head_length/2.
            else:
                head_width = 2*width
                head_length = head_width
            ax.arrow(p1[0], p1[1], (p2[0] - p1[0])/length*(length - head_length),
                     (p2[1] - p1[1])/length*(length - head_length),
                     head_width = head_width, fc = 'b', linewidth = 0,
                     head_length = head_length, width = width, alpha = 0.3)
        else:
            if width is None:
                width=1
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, marker='o', linewidth=width)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=width)
        return ax

    def To3D(self, plane_origin, x1, x2):
        p3D=[p.To3D(plane_origin,x1,x2) for p in self.points]
        return LineSegment3D(*p3D,self.name)
    
    def reverse(self):
        return LineSegment2D(self.points[1].copy(), self.points[0].copy())
    
    def to_line(self):
        return Line2D(*self.points)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return LineSegment2D(*[p.Rotation(center,angle,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center,angle,copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment2D(*[p.Translation(offset,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset,copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return LineSegment2D(*[frame.OldCoordinates(p) for p in self.points])
            else:
                self.points = [frame.OldCoordinates(p) for p in self.points]
        if side == 'new':
            if copy:
                return LineSegment2D(*[frame.NewCoordinates(p) for p in self.points])
            else:
                self.points = [frame.NewCoordinates(p) for p in self.points]

    def GeoScript(self, primitive_index, points_indices):
        s='Line({}) = {{{}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1

    def plot_data(self, marker=None, color='black', stroke_width=1,
                 dash=False, opacity=1, arrow=False):
        return {'type' : 'line',
                'data' : [self.points[0].vector[0], self.points[0].vector[1],
                          self.points[1].vector[0], self.points[1].vector[1]],
                'color' : color,
                'marker' : marker,
                'size' : stroke_width,
                'dash' : dash,
                'opacity' : opacity,
                'arrow': arrow
                }

    def CreateTangentCircle(self, point, other_line):
        circle1, circle2 = Line2D.CreateTangentCircle(other_line, point, self)
        if circle1 is not None:
            point_J1, curv_abs1 = Line2D.PointProjection(self, circle1.center, True)
            if curv_abs1 < 0. or curv_abs1 > 1.:
                circle1 = None
        if circle2 is not None:
            point_J2, curv_abs2 = Line2D.PointProjection(self, circle2.center, True)
            if curv_abs2 < 0. or curv_abs2 > 1.:
                circle2 = None
        return circle1, circle2


class Arc2D(Primitive2D):
    """
    angle: the angle measure always >= 0
    """
    def __init__(self, start, interior, end, name=''):
        Primitive2D.__init__(self, name)
        self.interior = interior
        self.start = start
        self.end = end
        xi, yi = interior.vector
        xe, ye = end.vector
        xs, ys = start.vector
        try : 
            A = Matrix22(2*(xs-xi), 2*(ys-yi),
                          2*(xs-xe), 2*(ys-ye))
            b = - Vector2D((xi**2 + yi**2 - xs**2 - ys**2,
                            xe**2 + ye**2 - xs**2 - ys**2))
            inv_A = A.inverse()
            x = inv_A.vector_multiplication(b)
            self.center = Point2D(x.vector)
        except ValueError : 
            A = npy.array([[2*(xs-xi), 2*(ys-yi)],
                           [2*(xs-xe), 2*(ys-ye)]])
            b = - npy.array([xi**2 + yi**2 - xs**2 - ys**2,
                               xe**2 + ye**2 - xs**2 - ys**2])
            self.center = Point2D(solve(A,b))

        r1 = self.start - self.center
        r2 = self.end - self.center
        ri = self.interior - self.center

        self.radius = r1.Norm()
        angle1 = math.atan2(r1.vector[1], r1.vector[0])
        anglei = math.atan2(ri.vector[1], ri.vector[0])
        angle2 = math.atan2(r2.vector[1], r2.vector[0])
        
        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei+2*math.pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + 2*math.pi
            
        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2+2*math.pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + 2*math.pi     
            
        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle1 = angle1
            self.angle2 = angle2
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle1 = angle2
            self.angle2 = angle1
            self.angle = clockwise_path
            
    
    def _get_points(self):
        return [self.start,self.interior,self.end]

    points=property(_get_points)


    def _get_geo_points(self):
        return [self.start,self.interior,self.end]

    geo_points=property(_get_geo_points)

    # def tessellation_points(self, resolution_for_circle=40):
    #     number_points_tesselation = math.ceil(resolution_for_circle*abs(self.angle)/2/math.pi)
    #     if number_points_tesselation == 1:
    #         number_points_tesselation += 1
            
    #     points = []
    #     if not self.is_trigo:
    #         delta_angle = -abs(self.angle1-self.angle2)/(number_points_tesselation-1)
    #         delta_angle = -self.angle/(number_points_tesselation-1)
    #     else:
    #         delta_angle =  abs(self.angle1-self.angle2)/(number_points_tesselation-1)
    #         delta_angle = self.angle/(number_points_tesselation-1)
    #     points.append(self.start)
    #     for i in range(number_points_tesselation-2):
    #         point_to_add = points[-1].Rotation(self.center, delta_angle)
    #         points.append(point_to_add)
    #     points.append(self.end)
    #     return points
    
    def tessellation_points(self, resolution_for_circle=40):
        # TODO: change this to a simple rotation?
        number_points_tesselation = math.ceil(resolution_for_circle*abs(self.angle)/2/math.pi)
        number_points_tesselation = max(number_points_tesselation, 5)
        
        # vector_start = Vector2D((self.start - self.center).vector)
        # vector_end = Vector2D((self.end - self.center).vector)
        # angle = clockwise_angle(vector_start, vector_end)
        # if not self.is_trigo:
        #     angle = - angle
        # delta_angle = angle/(number_points_tesselation-1)
        # points = []
        # points.append(self.start)
        # for i in range(number_points_tesselation-2):
        #     point_to_add = points[-1].Rotation(self.center, delta_angle)
        #     points.append(point_to_add)
        # points.append(self.end)
        l = self.Length()
        return [self.PointAtCurvilinearAbscissa(i/(number_points_tesselation-1)*l) for i in range(number_points_tesselation)]
        # return points
    
    def point_belongs(self, point):
        """
        Computes if the point belongs to the pizza slice drawn by the arc and its center
        """
        circle = Circle2D(self.center, self.radius)
        if not circle.point_belongs(point):
            return False
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end  = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = clockwise_angle(vector_start, vector_end)
        point_angle = clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return True

    def point_distance(self, point):
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end  = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = clockwise_angle(vector_start, vector_end)
        point_angle = clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return abs(LineSegment2D(point, self.center).Length()-self.radius)
        else:
            return min(LineSegment2D(point, self.start).Length(), LineSegment2D(point, self.end).Length())

    def line_intersection(self, line):
        circle = Circle2D(self.center, self.radius)
        circle_intersection_points = circle.line_intersection(line)
        
        if circle_intersection_points is None:
            return None
        
        intersection_points = []
        for pt in circle_intersection_points:
            if self.point_belongs(pt):
                intersection_points.append(pt)
        return intersection_points

    def Length(self):
        return self.radius * abs(self.angle)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        if self.is_trigo:
            return self.start.Rotation(self.center, curvilinear_abscissa/self.radius)
            # return self.start.Rotation(self.center, curvilinear_abscissa*self.angle)
        else:
            return self.start.Rotation(self.center, -curvilinear_abscissa/self.radius)
            # return self.start.Rotation(self.center, -curvilinear_abscissa*self.angle)

    def MiddlePoint(self):
        l = self.Length()
        return self.PointAtCurvilinearAbscissa(0.5*l)

    def GeoScript(self, primitive_index, points_indices):
        s='Circle({}) = {{{}, {}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1

    def Area(self):
        if self.angle2<self.angle1:
            angle=self.angle2+2*math.pi-self.angle1
        else:
            angle=self.angle2-self.angle1
        return self.radius**2*angle/2

    def CenterOfMass(self):
#        u=self.middle.vector-self.center.vector
        u = self.MiddlePoint() - self.center
        u.Normalize()
        alpha = abs(self.angle)
        return self.center + 4/(3*alpha)*self.radius*math.sin(alpha*0.5)*u

    def MPLPlot(self, ax=None, color='k'):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        pc = self.center.vector
#        ax.plot([pc[0]], [pc[1]], 'or')
#        ax.plot([self.interior[0]], [self.interior[1]], 'ob')
        ax.add_patch(Arc(pc, 2*self.radius, 2*self.radius, angle=0,
                    theta1=self.angle1*0.5/math.pi*360,
                    theta2=self.angle2*0.5/math.pi*360,
                    color=color))
        
        return ax

    def To3D(self,plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)

        return Arc3D(ps, pi, pe, name=self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Arc2D(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.interior,self.end]])
        else:
            self.__init__(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.interior,self.end]])
#            self.start.Rotation(center,angle,copy=False)
#            self.interior.Rotation(center,angle,copy=False)
#            self.end.Rotation(center,angle,copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Arc2D(*[p.Translation(offset,copy=True) for p in [self.start,self.interior,self.end]])
        else:
            self.__init__(*[p.Translation(offset,copy=True) for p in [self.start,self.interior,self.end]])
#            self.start.Translation(offset,copy=False)
#            self.interior.Translation(offset,copy=False)
#            self.end.Translation(offset,copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return Arc2D(*[p.frame_mapping(frame, side, copy=True) for p in [self.start,self.interior,self.end]])
        else:
            self.__init__(*[p.frame_mapping(frame, side, copy=True) for p in [self.start,self.interior,self.end]])
#            self.start.frame_mapping(frame, side,copy=False)
#            self.interior.frame_mapping(frame, side,copy=False)
#            self.end.frame_mapping(frame, side,copy=False)


    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        if self.angle2<self.angle1:
            angle2 = self.angle2+2*math.pi

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        Ix = self.radius**4/8*(angle2-angle1+0.5*(math.sin(2*angle1)-math.sin(2*angle2)))
        Iy = self.radius**4/8*(angle2-angle1+0.5*(math.sin(2*angle2)-math.sin(2*angle1)))
        Ixy = self.radius**4/8*(math.cos(angle1)**2-math.cos(angle2)**2)
        Ic = npy.array([[Ix, Ixy], [Ixy, Iy]])
        return geometry.Huygens2D(Ic, self.Area(), self.center, point)

    def Discretise(self, num=10):
        list_node = []
        if (self.angle1 < 0) and (self.angle2 > 0):
            delta_angle = -self.angle1 + self.angle2
        elif (self.angle1 > 0) and (self.angle2 < 0):
            delta_angle =  (2*npy.pi + self.angle2) - self.angle1
        else:
            delta_angle = self.angle2 - self.angle1
        for angle in npy.arange(self.angle1, self.angle1 + delta_angle, delta_angle/(num*1.)):
            list_node.append(Point2D(self.center + self.radius*Vector2D((npy.cos(angle), npy.sin(angle)))))
        list_node.append(Point2D(self.center + self.radius*Vector2D((npy.cos(self.angle1 + delta_angle), npy.sin(self.angle1 + delta_angle)))))
        if list_node[0] == self.start:
            return list_node
        else:
            return list_node[::-1]

    def plot_data(self, marker=None, color='black', stroke_width=1, dash=False, opacity=1):
        list_node = self.Discretise()
        data = []
        for nd in list_node:
            data.append({'x': nd.vector[0], 'y': nd.vector[1]})
        return {'type' : 'arc',
                    'cx' : self.center.vector[0],
                    'cy' : self.center.vector[1],
                    'data' : data,
                    'r' : self.radius,
                    'color' : color,
                    'opacity' : opacity,
                    'size' : stroke_width,
                    'dash' : None,
                    'marker' : marker,
                    'angle1' : self.angle1,
                    'angle2' : self.angle2, }

    def copy(self) :
        return Arc2D(self.start.copy(), self.interior.copy(), self.end.copy())

class ArcEllipse2D(Primitive2D) :
    """
    An arc is defined by a starting point, an end point and an interior point
    
    """
    def __init__(self, start, interior, end, center, major_dir, name='', extra=None):
        self.start = start
        self.interior = interior
        self.end = end
        self.center = center
        self.extra = extra
        self.major_dir = major_dir
        self.minor_dir = self.major_dir.deterministic_unit_normal_vector()
        
        frame = Frame2D(self.center, self.major_dir, self.minor_dir)
        start_new, end_new = frame.NewCoordinates(self.start), frame.NewCoordinates(self.end)
        interior_new, center_new = frame.NewCoordinates(self.interior), frame.NewCoordinates(self.center)
        
        #### from : https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a-center-and-3-arbitrary-points-on-it-are-given
        def theta_A_B(s,i,e,c): #theta=angle d'inclinaison ellipse par rapport à horizontal(sens horaire),A=demi grd axe, B=demi petit axe
            xs, ys, xi, yi, xe, ye = s[0]-c[0], s[1]-c[1], i[0]-c[0], i[1]-c[1], e[0]-c[0], e[1]-c[1]
            A = npy.array(([xs**2, ys**2, 2*xs*ys],
                          [xi**2, yi**2, 2*xi*yi],
                          [xe**2, ye**2, 2*xe*ye]))
            invA = npy.linalg.inv(A)
            One = npy.array(([1],
                            [1],
                            [1]))
            C = npy.dot(invA,One) #matrice colonne de taille 3
            theta = 0
            c1 = C[0]+C[1]
            c2 = (C[1]-C[0])/math.cos(2*theta)
            gdaxe = math.sqrt((2/(c1-c2)))
            ptax = math.sqrt((2/(c1+c2)))
            return theta, gdaxe, ptax
        
        if start==end :
            extra_new = frame.NewCoordinates(self.extra)
            theta, A, B = theta_A_B(start_new,extra_new,interior_new,center_new)
        else : 
            theta, A, B = theta_A_B(start_new,interior_new,end_new,center_new)
        
        self.Gradius = A
        self.Sradius = B
        self.theta = theta
        
        #Angle pour start
        u1, u2 = start_new.vector[0]/self.Gradius, start_new.vector[1]/self.Sradius
        angle1 = sin_cos_angle(u1, u2)
        #Angle pour end
        u3, u4 = end_new.vector[0]/self.Gradius, end_new.vector[1]/self.Sradius
        angle2 = sin_cos_angle(u3, u4)
        #Angle pour interior
        u5, u6 = interior_new.vector[0]/self.Gradius, interior_new.vector[1]/self.Sradius
        anglei = sin_cos_angle(u5, u6)
        
        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei+2*math.pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + 2*math.pi
            
        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2+2*math.pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + 2*math.pi     
            
        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path 
        
        if self.start == self.end or self.angle == 0:
            self.angle = 2*math.pi
        
        if self.is_trigo : #sens trigo
            self.offset_angle = angle1
        else :
            self.offset_angle = angle2
        
        Primitive2D.__init__(self, name=name)

    def _get_points(self):
        return self.tessellation_points()
    
    points=property(_get_points)

    def tessellation_points(self, resolution_for_ellipse=40):
        number_points_tesselation = math.ceil(resolution_for_ellipse*abs(0.5*self.angle/math.pi))
        
        frame2d = Frame2D(self.center, self.major_dir, self.minor_dir)
        
        tessellation_points_2D = [(Point2D((self.Gradius*math.cos(self.offset_angle + self.angle*i/(number_points_tesselation)), self.Sradius*math.sin(self.offset_angle + self.angle*i/(number_points_tesselation))))) for i in range(number_points_tesselation+1)]
        
        global_points = []
        for pt in tessellation_points_2D:
            global_points.append(frame2d.OldCoordinates(pt))
        
        return global_points

    def To3D(self,plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)
        pc = self.center.To3D(plane_origin, x, y)
        if self.extra is None :
            pextra = None
        else :
            pextra = self.extra.To3D(plane_origin, x, y)
        if ps == pe :
            p3 = pextra
        else :
            p3 = pe
        plane = Plane3D.from_3_points(ps, pi, p3)
        n = plane.normal
        major_dir = self.major_dir.To3D(plane_origin, x, y)
        major_dir.Normalize()
        
        return ArcEllipse3D(ps, pi, pe, pc, major_dir, normal=n, name=self.name, extra=pextra)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = ax.figure

        self.interior.MPLPlot(ax=ax, color='m')
        self.start.MPLPlot(ax=ax, color='r')
        self.end.MPLPlot(ax=ax, color='b')
        self.center.MPLPlot(ax=ax, color='y')

        x = []
        y = []
        for px, py in self.tessellation_points():
            x.append(px)
            y.append(py)

        plt.plot(x, y, 'k')
        return fig, ax

class Circle2D(Contour2D):
    _non_serializable_attributes  = ['internal_arcs', 'external_arcs',
                                     'polygon', 'straight_line_contour_polygon',
                                     'primitives', 'basis_primitives']
    def __init__(self,center,radius,name=''):
        self.center = center
        self.radius = radius
        self.angle = 2*math.pi
        self.utd_geo_points = False
        
        self.points = self.tessellation_points()
        
        Contour2D.__init__(self, [self], name=name) # !!! this is dangerous
    
    def __hash__ (self):
        return int(round(1e6*(self.center.vector[0] + self.center.vector[1] + self.radius)))

    def __eq__(self, other_circle):
        return math.isclose(self.center.vector[0], other_circle.center.vector[0], abs_tol=1e-06) \
           and math.isclose(self.center.vector[1], other_circle.center.vector[1], abs_tol=1e-06) \
           and math.isclose(self.radius, other_circle.radius, abs_tol=1e-06)
        
    def _get_geo_points(self):
        if not self.utd_geo_points:
            self._geo_start = self.center+self.radius*Point2D((1,0))
            self.utd_geo_points = True
        return [self._geo_start, self.center, self._geo_start]

    geo_points = property(_get_geo_points)

    def tessellation_points(self, resolution=40):
        return [self.center + self.radius*math.cos(teta)*Vector2D((1,0)) + self.radius*math.sin(teta)*Vector2D((0,1)) \
                for teta in npy.linspace(0, 2*math.pi, resolution+1)][:-1]

    def point_belongs(self, point):
        epsilon = 1e-6
        return point.point_distance(self.center) <= self.radius + epsilon
    
    def line_intersection(self, line):

        Q = Vector2D(self.center.vector)
        if line.points[0].vector == self.center.vector:
            P1 = Vector2D(line.points[1].vector)
            V = Vector2D((line.points[0] - line.points[1]).vector)
        else:
            P1 = Vector2D(line.points[0].vector)
            V = Vector2D((line.points[1] - line.points[0]).vector)
        
        a = V.Dot(V)
        b = 2 * V.Dot(P1 - Q)
        c = P1.Dot(P1) + Q.Dot(Q) - 2 * P1.Dot(Q) - self.radius**2
        
        disc = b**2 - 4 * a * c
        if disc < 0:
            return None

        sqrt_disc = math.sqrt(disc)
        t1 = (-b + sqrt_disc) / (2 * a)
        t2 = (-b - sqrt_disc) / (2 * a)
        if line.__class__ is Line2D:
            if t1 == t2:
                return [Point2D((P1+t1*V).vector)]
            else:
                return [Point2D((P1+t1*V).vector), Point2D((P1+t2*V).vector)]
        else:
            if not 0 <= t1 <= 1 and not 0 <= t2 <= 1:
                return None
            elif 0 <= t1 <= 1 and not 0 <= t2 <= 1:
                return [Point2D((P1+t1*V).vector)]
            elif not 0 <= t1 <= 1 and 0 <= t2 <= 1:
                return [Point2D((P1+t2*V).vector)]
            else:
                return [Point2D((P1+t1*V).vector), Point2D((P1+t2*V).vector)]

    def Length(self):
        return 2* math.pi * self.radius

    def GeoScript(self, primitive_index, points_indices):
        s = 'Circle({}) = {{{}, {}, {}}};\n'.format(primitive_index,*points_indices)
        return s, primitive_index+1

    def MPLPlot(self, ax=None, linestyle='-', color='k', linewidth=1):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        pc = self.center.vector
        if self.radius > 0:
            ax.add_patch(Arc(pc,
                             2*self.radius,
                             2*self.radius,
                             angle=0,
                             theta1=0,
                             theta2=360,
                             color=color,
                             linestyle=linestyle,
                             linewidth=linewidth))
        return ax

    def To3D(self, plane_origin, x, y):
        normal = x.Cross(y)
        pc = self.center.To3D(plane_origin, x, y)
        return Circle3D(pc, self.radius, normal, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Circle2D(self.center.Rotation(center,angle,copy=True),self.radius)
        else:
            self.center.Rotation(center,angle,copy=False)
            self.utd_geo_points=False

    def Translation(self,offset,copy=True):
        if copy:
            return Circle2D(self.center.Translation(offset,copy=True),self.radius)
        else:
            self.center.Translation(offset,copy=False)
            self.utd_geo_points=False

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return Circle2D(frame.OldCoordinates(self.center), self.radius)
            else:
                self.center = frame.OldCoordinates(self.center)
        if side == 'new':
            if copy:
                return Circle2D(frame.NewCoordinates(self.center), self.radius)
            else:
                self.points = frame.NewCoordinates(self.center)

    def Area(self):
        return math.pi*self.radius**2

    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        I = math.pi*self.radius**4/4
        Ic = npy.array([[I,0],[0,I]])
        return geometry.Huygens2D(Ic,self.Area(),self.center,point)

    def CenterOfMass(self):
        return self.center

    def point_symmetric(self, point):
        center = 2*point - self.center
        return Circle2D(center, self.radius)

    def plot_data(self, marker=None, color='black', stroke_width=1, opacity=1, fill=None):
        return {'type' : 'circle',
                'cx' : self.center.vector[0],
                'cy' : self.center.vector[1],
                'r' : self.radius,
                'color' : color,
                'opacity' : opacity,
                'size' : stroke_width,
                'dash' : None,
                'fill' : fill}
    
    def copy(self) :
        return Circle2D(self.center.copy(), self.radius)
        

class Polygon2D(Contour2D):
    # TODO: inherit from contour?
    def __init__(self,points, name=''):
        self.points = points
        # primitives = []
        # for p1,p2 in zip(points,points[1:]+[points[0]]):
        #     primitives.append(LineSegment2D(p1,p2))

        # TODO: remove this?
        self.line_segments = self._LineSegments()

        Contour2D.__init__(self, self.line_segments, name)

    def copy(self):
        points = [p.copy() for p in self.points]
        return Polygon2D(points, self.name)

    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def __eq__(self, other_):
        equal = True
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point == other_point)
        return equal

    def Area(self):

        x=[point.vector[0]for point in self.points]
        y=[point.vector[1]for point in self.points]

        return 0.5*npy.abs(npy.dot(x,npy.roll(y,1))-npy.dot(y,npy.roll(x,1)))

    def CenterOfMass(self):

        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]


        xi_xi1 = x+npy.roll(x,-1)
        yi_yi1 = y+npy.roll(y,-1)
        xi_yi1 = npy.multiply(x,npy.roll(y,-1))
        xi1_yi = npy.multiply(npy.roll(x,-1),y)

        a = 0.5*npy.sum(xi_yi1-xi1_yi)# signed area!
#        a=self.Area()
        if not math.isclose(a, 0, abs_tol=1e-08):
            cx = npy.sum(npy.multiply(xi_xi1,(xi_yi1-xi1_yi)))/6./a
            cy = npy.sum(npy.multiply(yi_yi1,(xi_yi1-xi1_yi)))/6./a
            return Point2D((cx, cy))

        else:
            raise NotImplementedError

    def PointBelongs(self, point):
        """
        Ray casting algorithm copied from internet...
        """
        return PolygonPointBelongs(point.vector, [p.vector for p in self.points])

    def SecondMomentArea(self, point):
        Ix, Iy, Ixy = 0, 0, 0
        for pi, pj in zip(self.points,self.points[1:]+[self.points[0]]):
            xi, yi = (pi-point).vector
            xj, yj = (pj-point).vector
            Ix += (yi**2 + yi*yj + yj**2)*(xi*yj - xj*yi)
            Iy += (xi**2 + xi*xj + xj**2)*(xi*yj - xj*yi)
            Ixy += (xi*yj + 2*xi*yi + 2*xj*yj + xj*yi)*(xi*yj - xj*yi)
        if Ix < 0:
            Ix =- Ix
            Iy =- Iy
            Ixy =- Ixy
        return npy.array([[Ix/12., Ixy/24.], [Ixy/24., Iy/12.]])

    def _LineSegments(self):
        lines=[]
        for p1, p2 in zip(self.points, self.points[1:]+[self.points[0]]):
            lines.append(LineSegment2D(p1,p2))
        return lines

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Polygon2D([p.Rotation(center,angle,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center,angle,copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Polygon2D([p.Translation(offset,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset,copy=False)

    def PointBorderDistance(self, point, return_other_point=False):
        """
        Compute the distance to the border distance of polygon
        Output is always positive, even if the point belongs to the polygon
        """
        d_min, other_point_min = self.line_segments[0].point_distance(point, return_other_point=True)
        for line in self.line_segments[1:]:
            d, other_point = line.point_distance(point, return_other_point=True)
            if d < d_min:
                d_min = d
                other_point_min = other_point
        if return_other_point:
            return d_min, other_point_min
        return d_min

    def SelfIntersect(self):
        epsilon = 0
        # BENTLEY-OTTMANN ALGORITHM
        # Sort the points along ascending x for the Sweep Line method
        sorted_index = sorted(range(len(self.points)), key=lambda p: (self.points[p][0], self.points[p][1]))
        nb = len(sorted_index)
        segments = []
        deleted = []

        while len(sorted_index) != 0: # While all the points haven't been swept
            # Stock the segments between 2 consecutive edges
            # Ex: for the ABCDE polygon, if Sweep Line is on C, the segments
            #   will be (C,B) and (C,D)
            if sorted_index[0]-1 < 0:
                segments.append((sorted_index[0], nb-1))
            else:
                segments.append((sorted_index[0], sorted_index[0]-1))
            if sorted_index[0] >= len(self.points)-1:
                segments.append((sorted_index[0], 0))
            else:
                segments.append((sorted_index[0], sorted_index[0]+1))


            # Once two edges linked by a segment have been swept, delete the
            # segment from the list
            to_del = []
            for index in deleted:
                if abs(index-sorted_index[0]) == 1 or abs(index-sorted_index[0]) == nb-1:
                    to_del.append((index, sorted_index[0]))
                    to_del.append((sorted_index[0], index))

            # Keep track of which edges have been swept
            deleted.append(sorted_index[0])
            sorted_index.pop(0)

            # Delete the segments that have just been swept
            index_to_del = []
            for i, segment in enumerate(segments):
                for seg_to_del in to_del:
                    if segment == seg_to_del:
                        index_to_del.append(i)
            for index in index_to_del[::-1]:
                segments.pop(index)

            # Checks if two segments are intersecting each other, returns True
            # if yes, otherwise the algorithm continues at WHILE
            for segment1 in segments:
                for segment2 in segments:
                    if segment1[0] != segment2[0] and segment1[1] != segment2[1] and segment1[0] != segment2[1] and segment1[1] != segment2[0]:

                        line1 = LineSegment2D(Point2D(self.points[segment1[0]]), Point2D(self.points[segment1[1]]))
                        line2 = LineSegment2D(Point2D(self.points[segment2[0]]), Point2D(self.points[segment2[1]]))

                        p, a, b = Point2D.LinesIntersection(line1, line2, True)

                        if p is not None:
                            if a >= 0+epsilon and a <= 1-epsilon and b >= 0+epsilon and b <= 1-epsilon:
                                return True, line1, line2

        return False, None, None


#    def Dict(self):
#        d = {'points': [point.Dict() for point in self.points], 'name':self.name}
#        return d
#
#    @classmethod
#    def DictToObject(cls, dict_):
#        return cls([Point2D.DictToObject(p) for p in dict_['points']], name=dict_['name'])



   
        
   
    def plot_data(self, marker=None, color='black', stroke_width=1, opacity=1):
        data = []
        for nd in self.points:
            data.append({'x': nd.vector[0], 'y': nd.vector[1]})
        return {'type' : 'wire',
                    'data' : data,
                    'color' : color,
                    'size' : stroke_width,
                    'dash' : None,
                    'marker' : marker,
                    'opacity' : opacity}

    @classmethod
    def points_convex_hull(cls, points):
        ymax, pos_ymax = max_pos([pt.vector[1] for pt in points])
        point_start = points[pos_ymax]
        hull, thetac = [point_start], 0 #thetac is the current theta
        
        barycenter = points[0]
        for pt in points[1:] :
            barycenter += pt
        barycenter = barycenter/(len(points))
        #second point of hull
        theta = []
        remaining_points = points
        del remaining_points[pos_ymax]
        
        vec1 = point_start - barycenter
        for pt in remaining_points :
            vec2 = pt - point_start
            theta_i = -clockwise_angle(vec1, vec2)
            theta.append(theta_i)
        
        min_theta, posmin_theta = min_pos(theta)
        thetac += min_theta
        next_point = remaining_points[posmin_theta]
        hull.append(next_point)
        del remaining_points[posmin_theta]
        #Adding first point to close the loop at the end
        remaining_points.append(hull[0])
        
        while next_point != point_start :
            vec1 = next_point - barycenter
            theta = []
            for pt in remaining_points :
                vec2 = pt - next_point
                theta_i = -clockwise_angle(vec1, vec2)
                theta.append(theta_i)
                
            min_theta, posmin_theta = min_pos(theta)
            thetac += min_theta
            next_point = remaining_points[posmin_theta]
            hull.append(next_point)
            del remaining_points[posmin_theta]
        
        hull.pop()
        
        return cls(hull)
        

class Primitive3D(dc.DessiaObject):
    def __init__(self, basis_primitives=None, name=''):
        self.name = name
        self.primitives = basis_primitives # une liste
        if basis_primitives is None:
            self.primitives = []
        
        dc.DessiaObject.__init__(self, name=name)
        
    def volmdlr_primitives(self):
        return [self]

class Plane3D(Primitive3D):
    def __init__(self, origin, vector1, vector2, name=''):
        self.origin = Point3D(origin.vector)
        vector1 = Vector3D(vector1.vector)
        vector1.Normalize()
        vector2 = Vector3D(vector2.vector)
        vector2.Normalize()
        self.vectors = [vector1, vector2]
        self.name = name
        self.normal = self.vectors[0].Cross(self.vectors[1])
        self.normal.Normalize()

    def __hash__(self):
        return sum([hash(v) for v in self.vectors]) + hash(self.origin)

    def __eq__(self, other_):
        equal = (self.origin == other_.origin
                 and self.vectors[0] == other_.vectors[0]
                 and self.vectors[1] == other_.vectors[1])
        return equal

    def to_dict(self):
        # improve the object structure ?
        dict_ = dc.DessiaObject.base_dict(self)
        dict_['vector1'] = self.vectors[0].to_dict()
        dict_['vector2'] = self.vectors[1].to_dict()
        dict_['origin'] = self.origin.to_dict()
        dict_['name'] = self.name
        dict_['object_class'] = 'volmdlr.core.Plane3D'
        return dict_

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        origin = frame3d.origin
        vector1 = frame3d.u
        vector2 = frame3d.v
        # TRANSFORMER EN 3D TOUS LES OBJETS LIES AU PLAN
        return cls(origin, vector1, vector2, arguments[0][1:-1])

    @classmethod
    def from_3_points(cls, point1, point2, point3):
        vector1 = point2 - point1
        vector2 = point3 - point1
        vector1.Normalize()
        vector2.Normalize()
        normal = vector1.Cross(vector2)
        normal.Normalize()
        vector = normal.Cross(vector1)
        return cls(point1.copy(), vector1.copy(), vector.copy())

    @classmethod
    def from_normal(cls, point, normal):
        v1 = normal.deterministic_unit_normal_vector()
        v2 = v1.Cross(normal)
        return cls(point, v1, v2)

    @classmethod
    def from_points(cls, points):
        if len(points) < 3:
            raise ValueError
        elif len(points) == 3:
            return cls.from_3_points(Point3D(points[0].vector), Vector3D(points[1].vector), Vector3D(points[2].vector))
        else:
            points = [p.copy() for p in points]
            indexes_to_del = []
            for i, point in enumerate(points[1:]):
                if point == points[0]:
                    indexes_to_del.append(i)
            for index in indexes_to_del[::-1]:
                del points[index+1]

            origin = Point3D(points[0].vector)
            vector1 = Vector3D(points[1]-origin)
            vector1.Normalize()
            vector2_min = Vector3D(points[2]-origin)
            vector2_min.Normalize()
            dot_min = abs(vector1.Dot(vector2_min))
            for point in points[3:]:
                vector2 = Vector3D(point-origin)
                vector2.Normalize()
                dot = abs(vector1.Dot(vector2))
                if dot < dot_min:
                    vector2_min = vector2
                    dot_min = dot
            return cls.from_3_points(origin, vector1+origin, vector2_min+origin)

    def point_on_plane(self, point):
        if math.isclose(self.normal.Dot(point-self.origin), 0, abs_tol=1e-6):
            return True
        return False

    def line_intersection(self, line):
        u = line.points[1] - line.points[0]
        w = line.points[0] - self.origin
        if math.isclose(self.normal.Dot(u), 0, abs_tol=1e-08):
            return None
        intersection_abscissea = - self.normal.Dot(w) / self.normal.Dot(u)
        return line.points[0] + intersection_abscissea * u

    def linesegment_intersection(self, linesegment, abscissea=False):
        u = linesegment.points[1] - linesegment.points[0]
        w = linesegment.points[0] - self.origin
        normalDotu = self.normal.Dot(u)
        if math.isclose(normalDotu, 0, abs_tol=1e-08):
            if abscissea:
                return None, None
            return None
        intersection_abscissea = - self.normal.Dot(w) / normalDotu
        if intersection_abscissea < 0 or intersection_abscissea > 1:
            if abscissea:
                return None, None
            return None
        if abscissea:
            return linesegment.points[0] + intersection_abscissea * u, intersection_abscissea
        return linesegment.points[0] + intersection_abscissea * u

    def equation_coefficients(self):
        """
        returns the a,b,c,d coefficient from equation ax+by+cz+d = 0
        """
        a, b, c = self.normal.vector
        d = -self.origin.Dot(self.normal)
        return (a, b, c, d)

    def plane_intersection(self, other_plane):
        line_direction = self.normal.Cross(other_plane.normal)

        if line_direction.Norm() < 1e-6:
            return None

        a1, b1, c1, d1 = self.equation_coefficients()
        a2, b2, c2, d2 = other_plane.equation_coefficients()

        if a1*b2-a2*b1 != 0.:
            x0 = (b1*d2-b2*d1)/(a1*b2-a2*b1)
            y0 = (a2*d1-a1*d2)/(a1*b2-a2*b1)
            point1 = Point3D((x0, y0, 0))
        else:
            y0 = (b2*d2-c2*d1)/(b1*c2-c1*b2)
            z0 = (c1*d1-b1*d2)/(b1*c2-c1*b2)
            point1 = Point3D((0, y0, z0))

        point2 = point1 + line_direction
        return Line3D(point1, point2)

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_origin = self.origin.Rotation(center, axis, angle, True)
            new_vector1 = self.vectors[0].Rotation(center, axis, angle, True)
            new_vector2 = self.vectors[1].Rotation(center, axis, angle, True)
            return Plane3D(new_origin, new_vector1, new_vector2, self.name)
        else:
            self.origin.Rotation(center, axis, angle, True)
            self.vectors[0].Rotation(center, axis, angle, True)
            self.vectors[1].Rotation(center, axis, angle, True)

    def Translation(self, offset, copy=True):
        if copy:
            new_origin = self.origin.Translation(offset, True)
            return Plane3D(new_origin, self.vectors[0], self.vectors[1], self.name)
        else:
            self.origin.Translation(offset, False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            new_origin = frame.OldCoordinates(self.origin)
            new_vector1 = frame.Basis().OldCoordinates(self.vectors[0])
            new_vector2 = frame.Basis().OldCoordinates(self.vectors[1])
            if copy:
                return Plane3D(new_origin, new_vector1, new_vector2, self.name)
            else:
                self.origin = new_origin
                self.vectors = [new_vector1, new_vector2]
                self.normal = frame.Basis().OldCoordinates(self.normal)
                self.normal.Normalize()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.origin)
            new_vector1 = frame.Basis().NewCoordinates(self.vectors[0])
            new_vector2 = frame.Basis().NewCoordinates(self.vectors[1])
            if copy:
                return Plane3D(new_origin, new_vector1, new_vector2, self.name)
            else:
                self.origin = new_origin
                self.vectors = [new_vector1, new_vector2]
                self.normal = frame.Basis().NewCoordinates(self.normal)
                self.normal.Normalize()

    def copy(self):
        new_origin = self.origin.copy()
        new_vector1 = self.vectors[0].copy()
        new_vector2 = self.vectors[1].copy()
        return Plane3D(new_origin, new_vector1, new_vector2, self.name)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        self.origin.MPLPlot(ax)
        self.vectors[0].MPLPlot(ax, starting_point=self.origin, color='r')
        self.vectors[1].MPLPlot(ax, starting_point=self.origin, color='g')
        return ax

    def Babylon(self):
        s = 'var myPlane = BABYLON.MeshBuilder.CreatePlane("myPlane", {width: 0.5, height: 0.5, sideOrientation: BABYLON.Mesh.DOUBLESIDE}, scene);\n'
        s += 'myPlane.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n'.format(self.origin[0], self.origin[1], self.origin[2])

        s += 'var axis1 = new BABYLON.Vector3({}, {}, {});\n'.format(self.vectors[0][0], self.vectors[0][1], self.vectors[0][2])
        s += 'var axis2 = new BABYLON.Vector3({}, {}, {});\n'.format(self.vectors[1][0], self.vectors[1][1], self.vectors[1][2])
        s += 'var axis3 = new BABYLON.Vector3({}, {}, {});\n'.format(self.normal[0], self.normal[1], self.normal[2])
        s += 'var orientation = BABYLON.Vector3.RotationFromAxis(axis1, axis2, axis3);\n'
        s += 'myPlane.rotation = orientation;\n'

        s += 'var planemat = new BABYLON.StandardMaterial("planemat", scene);\n'
        s += 'planemat.alpha = 0.4;\n'
        s += 'myPlane.material = planemat;\n'

        return s

PLANE3D_OXY = Plane3D(O3D, X3D, Y3D)
PLANE3D_OYZ = Plane3D(O3D, Y3D, Z3D)
PLANE3D_OZX = Plane3D(O3D, Z3D, X3D)



#    @classmethod
#    def DictToObject(cls, dict_):
#        vectors = [Vector3D.DictToObject(vector_dict) for vector_dict in dict_['vectors']]
#        return cls(*vectors)



XYZ = Basis3D(X3D, Y3D, Z3D)
YZX = Basis3D(Y3D, Z3D, X3D)
ZXY = Basis3D(Z3D, X3D, Y3D)




OXYZ = Frame3D(O3D, X3D, Y3D, Z3D)
OYZX = Frame3D(O3D, Y3D, Z3D, X3D)
OZXY = Frame3D(O3D, Z3D, X3D, Y3D)

class CylindricalSurface3D(Primitive3D):
    """
    :param frame: Cylinder's frame to position it 
    :type frame: Frame3D
    :param radius: Cylinder's radius
    :type radius: float
    """
    
    def __init__(self, frame, radius, name=''): 
        self.frame = frame
        self.radius = radius
        self.name = name
        V=frame.v
        V.Normalize()
        W=frame.w
        W.Normalize()
        self.plane = Plane3D(frame.origin,V,W)
        
    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, -frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = Frame3D(frame3d.origin, U, V, W)
        radius = float(arguments[2])/1000
        return cls(frame_direct, radius, arguments[0][1:-1])
    
    def frame_mapping(self, frame, side, copy=True) :
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return CylindricalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return CylindricalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

class ToroidalSurface3D(Primitive3D):
    """
    :param frame: Tore's frame to position it 
    :type frame: Frame3D
    :param rcenter: Tore's radius
    :type rcenter: float
    :param rcircle: Circle to revolute radius
    :type rcircle: float
    """
    def __init__(self, frame, rcenter, rcircle, name=''): 
        self.frame = frame
        self.rcenter = rcenter
        self.rcircle = rcircle
        self.name = name
        V=frame.v
        V.Normalize()
        W=frame.w
        W.Normalize()
        self.plane = Plane3D(frame.origin,V,W)
        
    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, -frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = Frame3D(frame3d.origin, U, V, W)
        rcenter = float(arguments[2])/1000
        rcircle = float(arguments[3])/1000
        return cls(frame_direct, rcenter, rcircle, arguments[0][1:-1])
    
    def frame_mapping(self, frame, side, copy=True) :
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ToroidalSurface3D(new_frame,self.rcenter,self.rcircle,name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ToroidalSurface3D(new_frame,self.rcenter,self.rcircle,name=self.name)
            else:
                self.frame = new_frame

class ConicalSurface3D(Primitive3D):
    """
    :param frame: Cone's frame to position it 
    :type frame: Frame3D
    :param r: Cone's bottom radius
    :type r: float
    :param semi_angle: Cone's semi-angle
    :type semi_angle: float
    """
    def __init__(self, frame, radius, semi_angle, name=''): 
        self.frame = frame
        self.radius = radius
        self.semi_angle = semi_angle
        self.name = name
        V=frame.v
        V.Normalize()
        W=frame.w
        W.Normalize()
        self.plane = Plane3D(frame.origin,V,W)
        
    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = Frame3D(frame3d.origin, U, V, W)
        radius = float(arguments[2])/1000
        semi_angle = float(arguments[3])
        return cls(frame_direct, radius, semi_angle, arguments[0][1:-1])
    
    def frame_mapping(self, frame, side, copy=True) :
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ConicalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ConicalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

class SphericalSurface3D(Primitive3D):
    """
    :param frame: Sphere's frame to position it 
    :type frame: Frame3D
    :param radius: Sphere's radius
    :type radius: float
    """
    
    def __init__(self, frame, radius, name=''): 
        self.frame = frame
        self.radius = radius
        self.name = name
        V=frame.v
        V.Normalize()
        W=frame.w
        W.Normalize()
        self.plane = Plane3D(frame.origin,V,W)
        
    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = Frame3D(frame3d.origin, U, V, W)
        radius = float(arguments[2])/1000
        return cls(frame_direct, radius, arguments[0][1:-1])

# class EllipseSurface3D(Primitive3D):
    
#     def __init__(self, frame, R, r, name=''): 
#         self.frame = frame
#         self.R = R
#         self.r = r
#         self.name = name
#         V=frame.v #U représente la normale
#         V.Normalize()
#         W=frame.w
#         W.Normalize()
#         self.plane = Plane3D(frame.origin,V,W)
        
#     @classmethod
#     def from_step(cls, arguments, object_dict):
#         frame3d = object_dict[arguments[1]]
#         R=arguments[2]
#         r=arguments[3]
#         return cls(frame3d, R, r, arguments[0][1:-1])
        
class Line3D(Primitive3D, Line):
    _non_eq_attributes = ['name', 'basis_primitives', 'bounding_box']

    """
    Define an infinite line passing through the 2 points
    """
    def __init__(self, point1, point2, name=''):
        Primitive3D.__init__(self, basis_primitives=[point1, point2], name=name)
        self.points = [point1, point2]
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return sum([hash(p) for p in self.points]) + hash(self.bounding_box)

    def to_dict(self):
        # improve the object structure ?
        dict_ = {}
        dict_['name'] = self.name
        dict_['point1'] = self.points[0].to_dict()
        dict_['point2'] = self.points[1].to_dict()
        return dict_

    def _bounding_box(self):
        points = self.points

        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])

        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + (self.points[1]-self.points[0]) * curvilinear_abscissa

    def MPLPlot(self, ax=None, color='k', dashed=True):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = ax.figure

        # Line segment
        x = [p.vector[0] for p in self.points]
        y = [p.vector[1] for p in self.points]
        z = [p.vector[2] for p in self.points]
        ax.plot(x,y,z, 'ok')

        # Drawing 3 times length of segment on each side
        u = self.points[1] - self.points[0]
        x1, y1, z1 = (self.points[0] - 3*u).vector
        x2, y2, z2 = (self.points[1] + 3*u).vector
        if dashed:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color, dashes=[30, 5, 10, 5])
        else:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color)            
        return ax

    def PlaneProjection2D(self,center, x, y):
        return Line2D(self.points[0].PlaneProjection2D(center, x, y),
                      self.points[1].PlaneProjection2D(center, x, y))

    def MinimumDistancePoints(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.points[1] - self.points[0]
        v = other_line.points[1] - other_line.points[0]
        w = self.points[0] - other_line.points[0]
        a = u.Dot(u)
        b = u.Dot(v)
        c = v.Dot(v)
        d = u.Dot(w)
        e = v.Dot(w)

        s = (b*e -c*d) / (a*c - b**2)
        t = (a*e -b*d) / (a*c - b**2)
        p1 = self.points[0] + s*u
        p2 = other_line.points[0] + t*v
        return p1, p2

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            return Line3D(*[p.Rotation(center, axis, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Line3D(*[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return Line3D(*[frame.OldCoordinates(p) for p in self.points])
            else:
                for p in self.points:
                    self.points = [frame.OldCoordinates(p) for p in self.points]
        if side == 'new':
            if copy:
                return Line3D(*[frame.NewCoordinates(p) for p in self.points])
            else:
                for p in self.points:
                    self.points = [frame.NewCoordinates(p) for p in self.points]

    def copy(self):
        return Line3D(*[p.copy() for p in self.points])

    @classmethod
    def from_step(cls, arguments, object_dict):
        point1 = object_dict[arguments[1]]
        direction = object_dict[arguments[2]]
        point2 = point1 + direction
        return cls(point1, point2, arguments[0][1:-1])

    def Intersection(self,line2):
        
        x1 = self.points[0].vector[0]
        y1 = self.points[0].vector[1]
        z1 = self.points[0].vector[2]
        x2 = self.points[1].vector[0]
        y2 = self.points[1].vector[1]
        z2 = self.points[1].vector[2]
        x3 = line2.points[0].vector[0]
        y3 = line2.points[0].vector[1]
        z3 = line2.points[0].vector[2]
        x4 = line2.points[1].vector[0]
        y4 = line2.points[1].vector[1]
        z4 = line2.points[1].vector[2]
        
        if x3 == 0 and x4 ==0 and y4-y3 == 0 :
            x5, y5, z5 = x3, y3, z3
            x6, y6, z6 = x4, y4, z4
            x3, y3, z3 = x1, y1, z1
            x4, y4, z4 = x2, y2, z2
            x1, y1, z1 = x5, y5, z5
            x2, y2, z2 = x6, y6, z6
        
        elif y3 == 0 and y4 ==0 and x4-x3 == 0 :
            x5, y5, z5 = x3, y3, z3
            x6, y6, z6 = x4, y4, z4
            x3, y3, z3 = x1, y1, z1
            x4, y4, z4 = x2, y2, z2
            x1, y1, z1 = x5, y5, z5
            x2, y2, z2 = x6, y6, z6
            
            
        res, list_t1 = [], []
        
        #2 unknown 3eq with t1 et t2 unknown
        
        if (x2-x1+y1-y2) != 0 and (y4-y3) != 0:
            t1 = (x3-x1 + (x4-x3)*(y1-y3)/(y4-y3))/(x2-x1+y1-y2)
            t2 = (y1-y3 + (y2-y1)*t1)/(y4-y3)
            res1 = z1 + (z2-z1)*t1
            res2 = z3 + (z4-z3)*t2
            list_t1.append(t1)
            res.append([res1, res2])
        
        if (z2-z1+y1-y2) != 0 and (y4-y3) != 0:
            t1 = (z3-z1 + (z4-z3)*(y1-y3)/(y4-y3))/(z2-z1+y1-y2)
            t2 = (y1-y3 + (y2-y1)*t1)/(y4-y3)
            res1 = x1 + (x2-x1)*t1
            res2 = x3 + (x4-x3)*t2
            list_t1.append(t1)
            res.append([res1, res2])
        
        if (z2-z1+x1-x2) != 0 and (x4-x3) != 0: 
            t1 = (z3-z1 + (z4-z3)*(x1-x3)/(x4-x3))/(z2-z1+x1-x2)
            t2 = (x1-x3 + (x2-x1)*t1)/(x4-x3)
            res1 = y1 + (y2-y1)*t1
            res2 = y3 + (y4-y3)*t2
            list_t1.append(t1)
            res.append([res1, res2])
        
        if len(res)==0 :
            return None
        
        for pair, t1 in zip(res, list_t1) :
            res1, res2 = pair[0], pair[1]
            if math.isclose(res1, res2, abs_tol=1e-7) : #if there is an intersection point
                return Point3D([x1+(x2-x1)*t1, y1+(y2-y1)*t1, z1+(z2-z1)*t1])
        
        return None

class BSplineCurve3D(Primitive3D):
    def __init__(self, degree, control_points, knot_multiplicities, knots, weights=None, periodic=False, name=''):
        Primitive3D.__init__(self, basis_primitives=control_points, name=name)
        self.control_points = control_points
        self.degree = degree
        knots = standardize_knot_vector(knots)
        self.knots = knots
        self.knot_multiplicities = knot_multiplicities
        self.weights = weights
        self.periodic = periodic
        self.name = name

        curve = BSpline.Curve()
        curve.degree = degree
        if weights is None:
            P = [(control_points[i][0], control_points[i][1], control_points[i][2]) for i in range(len(control_points))]
            curve.ctrlpts = P
        else:
            Pw = [(control_points[i][0]*weights[i], control_points[i][1]*weights[i], control_points[i][2]*weights[i], weights[i]) for i in range(len(control_points))]
            curve.ctrlptsw = Pw
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot]*knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.1
        curve_points = curve.evalpts

        self.curve = curve
        self.points = [Point3D((p[0], p[1], p[2])) for p in curve_points]

    def Length(self):
        #Approximately
        length = 0
        for k in range(0,len(self.points)-1) :
            length += (self.points[k] - self.points[k+1]).Norm()
        return length
    
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        #copy paste from wire3D
        length = 0.
        primitives = []
        for k in range(0, len(self.points)-1) :
            primitives.append(LineSegment3D(self.points[k], self.points[k+1]))
        for primitive in primitives:
            primitive_length = primitive.Length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        points = '['
        for i in range(len(self.control_points)):
            point = 'fc.Vector({},{},{}),'.format(self.control_points[i][0],self.control_points[i][1],self.control_points[i][2])
            points += point
        points = points[:-1]
        points += ']'
        # !!! : A QUOI SERT LE DERNIER ARG DE BSplineCurve (False)?
        # LA MULTIPLICITE EN 3e ARG ET LES KNOTS EN 2e ARG ?
        return '{} = Part.BSplineCurve({},{},{},{},{},{},{})\n'.format(name,points,self.knot_multiplicities,self.knots,self.periodic,self.degree,self.weights,False)

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        degree = int(arguments[1])
        points = [object_dict[int(i[1:])] for i in arguments[2]]
        # curve_form = arguments[3]
        if arguments[4] == '.F.':
            closed_curve = False
        elif arguments[4] == '.T.':
            closed_curve = True
        else:
            raise ValueError
        # self_intersect = arguments[5]
        knot_multiplicities = [int(i) for i in arguments[6][1:-1].split(",")]
        knots = [float(i) for i in arguments[7][1:-1].split(",")]
        # knot_spec = arguments[8]
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot]*knot_multiplicities[i])

        if 9 in range(len(arguments)):
            weight_data = [float(i) for i in arguments[9][1:-1].split(",")]
        else:
            weight_data = None

        # FORCING CLOSED_CURVE = FALSE:
        closed_curve = False
        return cls(degree, points, knot_multiplicities, knots, weight_data, closed_curve, name)

    def point_distance(self, pt1):
        distances = []
        for point in self.points:
#            vmpt = Point3D((point[1], point[2], point[3]))
            distances.append(pt1.point_distance(point))
        return min(distances)

    def Rotation(self, center, axis, angle, copy=True):
        new_control_points = [p.Rotation(center, axis, angle, True) for p in self.control_points]
        new_BSplineCurve3D = BSplineCurve3D(self.degree, new_control_points, self.knot_multiplicities, self.knots, self.weights, self.periodic, self.name)
        if copy:
            return new_BSplineCurve3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineCurve3D.curve
            self.points = new_BSplineCurve3D.points

    def Translation(self, offset, copy=True):
        new_control_points = [p.Translation(offset, True) for p in self.control_points]
        new_BSplineCurve3D = BSplineCurve3D(self.degree, new_control_points, self.knot_multiplicities, self.knots, self.weights, self.periodic, self.name)
        if copy:
            return new_BSplineCurve3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineCurve3D.curve
            self.points = new_BSplineCurve3D.points
    #Copy paste du LineSegment3D
    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        x=[p.vector[0] for p in self.points]
        y=[p.vector[1] for p in self.points]
        z=[p.vector[2] for p in self.points]
        ax.plot(x,y,z, 'o-k')
        return ax

    def To2D(self, plane_origin, x1, x2):
        control_points2D=[p.To2D(plane_origin,x1,x2) for p in self.control_points]
        return BSplineCurve2D(self.degree, control_points2D, self.knot_multiplicities, self.knots, self.weights, self.periodic, self.name)

    def tessellation_points(self) :
        return self.points

class Arc3D(Primitive3D):
    """
    An arc is defined by a starting point, an end point and an interior point
    
    """
    def __init__(self, start, interior, end, normal=None, name='', other_vec=None):
        """
        TODO : vérifier la position du centre et valeur du rayon pas seulement pour le cas particulier (s=e et i à 180°)    
    
        """
        self.start = start
        self.interior = interior
        self.end = end
        self.other_vec = other_vec
        self.setup_arc(start, interior, end, normal=normal, name=name)

    def setup_arc(self, start, interior, end, normal=None, name='') :
        u1 = (self.interior - self.start)
        u2 = (self.interior - self.end)
        try:
            u1.Normalize()
            u2.Normalize()
        except ZeroDivisionError:
            raise ValueError('Start, end and interior points  of an arc must be distincts')
        
        if normal is None:
            n = u2.Cross(u1)
            n.Normalize()
            self.normal = n
        else:
            normal.Normalize()
            self.normal = normal
        
        if u1 == u2:
            u2=self.normal.Cross(u1) 
            u2.Normalize()
            
        v1 = self.normal.Cross(u1)# v1 is normal, equal u2
        v2 = self.normal.Cross(u2)#equal -u1
        
        p11 = 0.5 * (start + interior)# Mid point of segment s,m
        p12 = p11 + v1
        p21 = 0.5 * (end + interior)# Mid point of segment s,m
        p22 = p21 + v2

        l1 = Line3D(p11, p12)
        l2 = Line3D(p21, p22)
        
        c1, _ = l1.MinimumDistancePoints(l2)
        
        self.center = c1 
        self.radius = (self.center - self.start).Norm()

        # Determining angle

        if self.other_vec is None :
            vec1 = (self.start - self.center)
        else :
            vec1 = self.other_vec
        vec1.Normalize()
        vec2 = self.normal.Cross(vec1)
            
        # r1 = (self.start).To2D(self.center, u1, v1)
        # r2 = (self.end).To2D(self.center, u1, v1)
        # ri = (self.interior).To2D(self.center, u1, v1)
        
        r1 = (self.start).To2D(self.center, vec1, vec2)
        r2 = (self.end).To2D(self.center, vec1, vec2)
        ri = (self.interior).To2D(self.center, vec1, vec2)
        
        angle1 = math.atan2(r1.vector[1], r1.vector[0])
        anglei = math.atan2(ri.vector[1], ri.vector[0])
        angle2 = math.atan2(r2.vector[1], r2.vector[0])
        
        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei+2*math.pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + 2*math.pi
            
        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2+2*math.pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + 2*math.pi     
            
        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path                
        Primitive3D.__init__(self, basis_primitives=self.tessellation_points(), name=name)
        
    # @classmethod
    # def from_step(cls, arguments, object_dict):
    #     edges = []
    #     for edge in arguments[1]:
    #         edges.append(object_dict[int(edge[1:])])

    def _get_points(self):
        return self.tessellation_points()
        # return [self.start,self.interior,self.end]

    points=property(_get_points)

    def tessellation_points(self, resolution_for_circle=40):
#        number_points_tesselation = math.ceil(resolution_for_circle*abs(0.5*self.angle/math.pi))
        number_points_tesselation = resolution_for_circle
        l = self.Length()
        tessellation_points_3D = [self.PointAtCurvilinearAbscissa(l*i/(number_points_tesselation)) for i in range(number_points_tesselation+1)]
        return tessellation_points_3D

    def Length(self):
        return self.radius * abs(self.angle)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.start.Rotation(self.center, self.normal, curvilinear_abscissa/self.radius, copy=True)

    def Rotation(self, rot_center, axis, angle, copy=True):
        if copy:
            new_start = self.start.Rotation(rot_center, axis, angle, True)
            new_interior = self.interior.Rotation(rot_center, axis, angle, True)
            new_end = self.end.Rotation(rot_center, axis, angle, True)
            return Arc3D(new_start, new_interior, new_end, name=self.name)
        else:
            self.center.Rotation(rot_center, axis, angle, False)
            self.start.Rotation(rot_center, axis, angle, False)
            self.interior.Rotation(rot_center, axis, angle, False)
            self.end.Rotation(rot_center, axis, angle, False)
            [p.Rotation(rot_center, axis, angle, False) for p in self.primitives]

    def Translation(self, offset, copy=True):
        if copy:
            new_start = self.start.Translation(offset, True)
            new_interior = self.interior.Translation(offset, True)
            new_end = self.end.Translation(offset, True)
            return Arc3D(new_start, new_interior, new_end, name=self.name)
        else:
            self.center.Translation(offset, False)
            self.start.Translation(offset, False)
            self.interior.Translation(offset, False)
            self.end.Translation(offset, False)
            [p.Translation(offset, False) for p in self.primitives]

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]], color='b')
        ax.plot([self.start[0]], [self.start[1]], [self.start[2]], c='r')
        ax.plot([self.end[0]], [self.end[1]], [self.end[2]], c='r')
        ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]], c='g')
        x = []
        y = []
        z = []
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)

        ax.plot(x, y, z, 'k')
        return ax

    def MPLPlot2D(self, center=O3D, x3d=X3D, y3D=Y3D, ax=None, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        # TODO: Enhance this plot
        l = self.Length()
        x = []
        y = []
        for i in range(30):
            p = self.PointAtCurvilinearAbscissa(i/(29.)*l)
            xi, yi = p.PlaneProjection2D(center, X3D, Y3D)
            x.append(xi)
            y.append(yi)
        ax.plot(x, y, color=color)

        return ax

    def FreeCADExport(self, name, ndigits=6):
        xs, ys, zs = round(1000*self.start, ndigits).vector
        xi, yi, zi = round(1000*self.interior, ndigits).vector
        xe, ye, ze = round(1000*self.end, ndigits).vector
        return '{} = Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,xs,ys,zs,xi,yi,zi,xe,ye,ze)

    def copy(self):
        return Arc3D(self.start.copy(), self.interior.copy(), self.end.copy())

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.OldCoordinates(self.start.copy())
            new_interior = frame.OldCoordinates(self.interior.copy())
            new_end = frame.OldCoordinates(self.end.copy())
            if copy:
                return Arc3D(new_start, new_interior, new_end, normal=None, name=self.name)
            else:
                self.start, self.interior, self.end = new_start, new_interior, new_end
                self.setup_arc(self.start, self.interior, self.end)
                
        if side == 'new':
            new_start = frame.NewCoordinates(self.start.copy())
            new_interior = frame.NewCoordinates(self.interior.copy())
            new_end = frame.NewCoordinates(self.end.copy())
            if copy:
                return Arc3D(new_start, new_interior, new_end, normal=None, name=self.name)
            else:
                self.start, self.interior, self.end = new_start, new_interior, new_end
                self.setup_arc(self.start, self.interior, self.end)
                
    def To2D(self,plane_origin, x, y):
        ps = self.start.To2D(plane_origin, x, y)
        pi = self.interior.To2D(plane_origin, x, y)
        pe = self.end.To2D(plane_origin, x, y)
        if ps==pe :
            pc = self.center.To2D(plane_origin, x, y)
            return Circle2D(pc, self.radius, name=self.name)
        else :
            return Arc2D(ps, pi, pe, name=self.name)
 
    def minimum_distance_points_arc(self, other_arc) :
    
        u1 = self.start - self.center
        u1.Normalize()
        u2 = self.normal.Cross(u1)
        
        w = other_arc.center - self.center
        
        u3 = other_arc.start - other_arc.center
        u3.Normalize()
        u4 = other_arc.normal.Cross(u3)
        
        r1, r2 = self.radius, other_arc.radius

        a, b, c, d = u1.Dot(u1), u1.Dot(u2), u1.Dot(u3), u1.Dot(u4) 
        e, f, g = u2.Dot(u2), u2.Dot(u3), u2.Dot(u4)
        h, i = u3.Dot(u3), u3.Dot(u4)
        j = u4.Dot(u4)
        k, l, m, n, o = w.Dot(u1), w.Dot(u2), w.Dot(u3), w.Dot(u4), w.Dot(w)
        
        # x = (theta1, theta2)
        def distance_squared(x):
            return (a*((math.cos(x[0]))**2)*r1**2 + e*((math.sin(x[0]))**2)*r1**2
                    + o + h*((math.cos(x[1]))**2)*r2**2 + j*((math.sin(x[1]))**2)*r2**2
                    + b*math.sin(2*x[0])*r1**2 - 2*r1*math.cos(x[0])*k
                    - 2*r1*r2*math.cos(x[0])*math.cos(x[1])*c
                    - 2*r1*r2*math.cos(x[0])*math.sin(x[1])*d - 2*r1*math.sin(x[0])*l
                    - 2*r1*r2*math.sin(x[0])*math.cos(x[1])*f
                    - 2*r1*r2*math.sin(x[0])*math.sin(x[1])*g + 2*r2*math.cos(x[1])*m
                    + 2*r2*math.sin(x[1])*n + i*math.sin(2*x[1])*r2**2)
                    
        
        x01 = npy.array([self.angle/2, other_arc.angle/2])
        
        res1 = scp.optimize.least_squares(distance_squared, x01, bounds=[(0,0), (self.angle,other_arc.angle)])
            
        p1 = self.PointAtCurvilinearAbscissa(res1.x[0]*r1)
        p2 = other_arc.PointAtCurvilinearAbscissa(res1.x[1]*r2)
        
        return p1, p2
    
    def minimum_distance_points_line(self, other_line) :
    
        u = other_line.DirectionVector()
        k = self.start - self.center
        k.Normalize()
        w = self.center - other_line.points[0] 
        v = self.normal.Cross(k)
        
        r = self.radius

        a = u.Dot(u)
        b = u.Dot(v)
        c = u.Dot(k)
        d = v.Dot(v)
        e = v.Dot(k)
        f = k.Dot(k)
        g = w.Dot(u)
        h = w.Dot(v)
        i = w.Dot(k)
        j = w.Dot(w)
        
        # x = (s, theta)
        def distance_squared(x):
            return (a*x[0]**2 + j + d*((math.sin(x[1]))**2)*r**2 + f*((math.cos(x[1]))**2)*r**2
                    - 2*x[0]*g - 2*x[0]*r*math.sin(x[1])*b - 2*x[0]*r*math.cos(x[1])*c
                    + 2*r*math.sin(x[1])*h + 2*r*math.cos(x[1])*i
                    + math.sin(2*x[1])*e*r**2)
        x01 = npy.array([0.5, self.angle/2])
        x02 = npy.array([0.5, 0])
        x03 = npy.array([0.5, self.angle])
        
        res1 = scp.optimize.least_squares(distance_squared, x01, bounds=[(0,0), (1,self.angle)])
        res2 = scp.optimize.least_squares(distance_squared, x02, bounds=[(0,0), (1,self.angle)])
        res3 = scp.optimize.least_squares(distance_squared, x03, bounds=[(0,0), (1,self.angle)])
            
        p1 = other_line.PointAtCurvilinearAbscissa(res1.x[0]*other_line.Length())
        p2 = self.PointAtCurvilinearAbscissa(res1.x[1]*r)
        # d1 = p1.point_distance(p2)
        
        res = [res2, res3]
        for couple in res :
            ptest1 = other_line.PointAtCurvilinearAbscissa(couple.x[0]*other_line.Length())
            ptest2 = self.PointAtCurvilinearAbscissa(couple.x[1]*r)
            dtest = ptest1.point_distance(ptest2)
            if dtest < d :
                p1, p2 = ptest1, ptest2
        
        return p1, p2
               
    def minimum_distance(self, element, return_points=False) :
        if element.__class__ is Arc3D or element.__class__ is Circle3D :
            p1, p2 = self.minimum_distance_points_arc(element)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
            
        elif element.__class__ is LineSegment3D :
            pt1, pt2 = self.minimum_distance_points_line(element)
            if return_points : 
                return pt1.point_distance(pt2), pt1, pt2
            else :
                return pt1.point_distance(pt2)
        else :
            return NotImplementedError

class ArcEllipse3D(Primitive3D) :
    """
    An arc is defined by a starting point, an end point and an interior point
    
    """
    def __init__(self, start, interior, end, center, major_dir, normal=None, name='', extra=None):
        #Extra is an additionnal point if start=end because you need 3 points on the arcellipse to define it
        self.start = start
        self.interior = interior
        self.end = end
        self.center = center
        major_dir.Normalize()
        self.major_dir = major_dir #Vector for Gradius
        self.extra = extra
        
        u1 = (self.interior - self.start)
        u2 = (self.interior - self.end)
        u1.Normalize()
        u2.Normalize()
        
        if u1 == u2:
            u2 = (self.interior - self.extra)
            u2.Normalize()
        
        if normal is None:
            n = u2.Cross(u1)
            n.Normalize()
            self.normal = n
        else:
            n = normal
            n.Normalize()
            self.normal = normal
        
        self.minor_dir = self.normal.Cross(self.major_dir)
        
        frame = Frame3D(self.center, self.major_dir, self.minor_dir, self.normal)
        start_new, end_new = frame.NewCoordinates(self.start), frame.NewCoordinates(self.end)
        interior_new, center_new = frame.NewCoordinates(self.interior), frame.NewCoordinates(self.center)
        
        #### from : https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a-center-and-3-arbitrary-points-on-it-are-given
        def theta_A_B(s,i,e,c): #theta=angle d'inclinaison ellipse par rapport à horizontal(sens horaire),A=demi grd axe, B=demi petit axe
            xs, ys, xi, yi, xe, ye = s[0]-c[0], s[1]-c[1], i[0]-c[0], i[1]-c[1], e[0]-c[0], e[1]-c[1]
            A = npy.array(([xs**2, ys**2, 2*xs*ys],
                          [xi**2, yi**2, 2*xi*yi],
                          [xe**2, ye**2, 2*xe*ye]))
            invA = npy.linalg.inv(A)
            One = npy.array(([1],
                            [1],
                            [1]))
            C = npy.dot(invA,One) #matrice colonne de taille 3
            theta = 0.5*math.atan(2*C[2]/(C[1]-C[0]))
            c1 = C[0]+C[1]
            c2 = (C[1]-C[0])/math.cos(2*theta)
            gdaxe = math.sqrt((2/(c1-c2)))
            ptax = math.sqrt((2/(c1+c2)))
            return theta, gdaxe, ptax
        
        if start==end :
            extra_new = frame.NewCoordinates(self.extra)
            theta, A, B = theta_A_B(start_new,extra_new,interior_new,center_new)
        else : 
            theta, A, B = theta_A_B(start_new,interior_new,end_new,center_new)
        
        self.Gradius = A
        self.Sradius = B
        self.theta = theta
        
        #Angle pour start
        u1, u2 = start_new.vector[0]/self.Gradius, start_new.vector[1]/self.Sradius
        angle1 = sin_cos_angle(u1, u2)
        #Angle pour end
        u3, u4 = end_new.vector[0]/self.Gradius, end_new.vector[1]/self.Sradius
        angle2 = sin_cos_angle(u3, u4)
        #Angle pour interior
        u5, u6 = interior_new.vector[0]/self.Gradius, interior_new.vector[1]/self.Sradius
        anglei = sin_cos_angle(u5, u6)
        
        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei+2*math.pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + 2*math.pi
            
        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2+2*math.pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + 2*math.pi     
            
        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path 
            
        if self.start == self.end :
            self.angle = 2*math.pi

        if self.is_trigo : 
            self.offset_angle = angle1
        else :
            self.offset_angle = angle2
    
        Primitive3D.__init__(self, basis_primitives=self.tessellation_points(), name=name)
        
    
    def _get_points(self):
        return self.tessellation_points()
    
    points=property(_get_points)
    
    def tessellation_points(self, resolution_for_ellipse=40):
        number_points_tesselation = math.ceil(resolution_for_ellipse*abs(0.5*self.angle/math.pi))
        
        plane3d = Plane3D(self.center, self.major_dir, self.minor_dir, self.normal)
        frame3d = Frame3D(self.center, plane3d.vectors[0], plane3d.vectors[1], plane3d.normal)
        
        tessellation_points_3D = [Point3D((self.Gradius*math.cos(self.offset_angle + self.angle*i/(number_points_tesselation)), self.Sradius*math.sin(self.offset_angle + self.angle*i/(number_points_tesselation)), 0)) for i in range(number_points_tesselation+1)]
        
        global_points = []
        for pt in tessellation_points_3D:
            global_points.append(frame3d.OldCoordinates(pt))
        
        return global_points
    
    def To2D(self,plane_origin, x, y):
        ps = self.start.To2D(plane_origin, x, y)
        pi = self.interior.To2D(plane_origin, x, y)
        pe = self.end.To2D(plane_origin, x, y)
        center = self.center.To2D(plane_origin, x, y)
        
        if self.extra is None :
            pextra = None
        else :
            pextra = self.extra.To2D(plane_origin, x, y)
        
        maj_dir2d = self.major_dir.To2D(plane_origin, x, y)
        maj_dir2d.Normalize()
        return ArcEllipse2D(ps, pi, pe, center, maj_dir2d, name=self.name, extra=pextra)
    
    def Length(self):
        return self.angle*math.sqrt((self.Gradius**2+self.Sradius**2)/2)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]], color='b')
        ax.plot([self.start[0]], [self.start[1]], [self.start[2]], c='r')
        ax.plot([self.end[0]], [self.end[1]], [self.end[2]], c='r')
        ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]], c='g')
        x = []
        y = []
        z = []
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)

        ax.plot(x, y, z, 'k')
        return ax

    def MPLPlot2D(self, x3d, y3D, ax, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        # TODO: Enhance this plot
        l = self.Length()
        x = []
        y = []
        for i in range(30):
            p = self.PointAtCurvilinearAbscissa(i/(29.)*l)
            xi, yi = p.PlaneProjection2D(X3D, Y3D)
            x.append(xi)
            y.append(yi)
        ax.plot(x, y, color=color)

        return ax

    def FreeCADExport(self, name, ndigits=6):
        xs, ys, zs = round(1000*self.start, ndigits).vector
        xi, yi, zi = round(1000*self.interior, ndigits).vector
        xe, ye, ze = round(1000*self.end, ndigits).vector
        return '{} = Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,xs,ys,zs,xi,yi,zi,xe,ye,ze)
 
class BSplineSurface3D(Primitive3D):
    def __init__(self, degree_u, degree_v, control_points, nb_u, nb_v, u_multiplicities, v_multiplicities, u_knots, v_knots, weights=None, name=''):
        Primitive3D.__init__(self, basis_primitives=control_points, name=name)
        self.control_points = control_points
        self.degree_u = degree_u
        self.degree_v = degree_v
        self.nb_u = nb_u
        self.nb_v = nb_v
        u_knots = standardize_knot_vector(u_knots)
        v_knots = standardize_knot_vector(v_knots)
        self.u_knots = u_knots
        self.v_knots = v_knots
        self.u_multiplicities = u_multiplicities
        self.v_multiplicities = v_multiplicities
        self.weights = weights


        self.control_points_table = []
        points_row = []
        i = 1
        for pt in control_points:
            points_row.append(pt)
            if i == nb_v:
                self.control_points_table.append(points_row)
                points_row = []
                i = 1
            else:
                i += 1
        surface = BSpline.Surface()
        surface.degree_u = degree_u
        surface.degree_v = degree_v
        if weights is None:
            P = [(control_points[i][0], control_points[i][1],
                  control_points[i][2]) for i in range(len(control_points))]
            surface.set_ctrlpts(P, nb_u, nb_v)
        else:
            Pw = [(control_points[i][0]*weights[i],
                   control_points[i][1]*weights[i], control_points[i][2]*weights[i],
                   weights[i]) for i in range(len(control_points))]
            surface.set_ctrlpts(Pw, nb_u, nb_v)
        knot_vector_u = []
        for i, u_knot in enumerate(u_knots):
            knot_vector_u.extend([u_knot]*u_multiplicities[i])
        knot_vector_v = []
        for i, v_knot in enumerate(v_knots):
            knot_vector_v.extend([v_knot]*v_multiplicities[i])
        surface.knotvector_u = knot_vector_u
        surface.knotvector_v = knot_vector_v
        surface.delta = 0.05
        surface_points = surface.evalpts

        self.surface = surface
        self.points = [Point3D((p[0], p[1], p[2])) for p in surface_points]
        
    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        script = ""
        points = '['
        for i, pts_row in enumerate(self.control_points_table):
            pts = '['
            for j, pt in enumerate(pts_row):
                point = 'fc.Vector({},{},{}),'.format(pt[0],pt[1],pt[2])
                pts += point
            pts = pts[:-1] + '],'
            points += pts
        points = points[:-1] + ']'

        script += '{} = Part.BSplineSurface()\n'.format(name)
        if self.weights is None:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={})\n'.format(name,points,self.u_multiplicities,self.v_multiplicities,self.degree_u,self.degree_v,self.u_knots,self.v_knots)
        else:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={},weights={})\n'.format(name,points,self.u_multiplicities,self.v_multiplicities,self.degree_u,self.degree_v,self.u_knots,self.v_knots,self.weights)

        return script
    
    def Rotation(self, center, axis, angle, copy=True):
        new_control_points = [p.Rotation(center, axis, angle, True) for p in self.control_points]
        new_BSplineSurface3D = BSplineSurface3D(self.degree_u, self.degree_v, new_control_points, self.nb_u, self.nb_v, self.u_multiplicities, self.v_multiplicities, self.u_knots, self.v_knots, self.weights, self.name)
        if copy:
            return new_BSplineSurface3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineSurface3D.curve
            self.points = new_BSplineSurface3D.points

    def Translation(self, offset, copy=True):
        new_control_points = [p.Translation(offset, True) for p in self.control_points]
        new_BSplineSurface3D = BSplineSurface3D(self.degree_u, self.degree_v, new_control_points, self.nb_u, self.nb_v, self.u_multiplicities, self.v_multiplicities, self.u_knots, self.v_knots, self.weights, self.name)
        if copy:
            return new_BSplineSurface3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineSurface3D.curve
            self.points = new_BSplineSurface3D.points

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        
        degree_u = int(arguments[1])
        degree_v = int(arguments[2])
        points_sets = arguments[3][1:-1].split("),")
        points_sets = [elem+")" for elem in points_sets[:-1]]+[points_sets[-1]]
        control_points = []
        for points_set in points_sets:
            points = [object_dict[int(i[1:])] for i in points_set[1:-1].split(",")]
            nb_v = len(points)
            control_points.extend(points)
        nb_u = int(len(control_points) / nb_v)
        surface_form = arguments[4]
        if arguments[5] == '.F.':
            u_closed = False
        elif arguments[5] == '.T.':
            u_closed = True
        else:
            raise ValueError
        if arguments[6] == '.F.':
            v_closed = False
        elif arguments[6] == '.T.':
            v_closed = True
        else:
            raise ValueError
        self_intersect = arguments[7]
        u_multiplicities = [int(i) for i in arguments[8][1:-1].split(",")]
        v_multiplicities = [int(i) for i in arguments[9][1:-1].split(",")]
        u_knots = [float(i) for i in arguments[10][1:-1].split(",")]
        v_knots = [float(i) for i in arguments[11][1:-1].split(",")]
        knot_spec = arguments[12]
    
        if 13 in range(len(arguments)):
            weight_data = [float(i) for i in arguments[13][1:-1].replace("(", "").replace(")", "").split(",")]
        else:
            weight_data = None
    
        return cls(degree_u, degree_v, control_points, nb_u, nb_v, u_multiplicities, v_multiplicities, u_knots, v_knots, weight_data, name)
    
    
class BSplineExtrusion(Primitive3D):

    def __init__(self, obj, vectorextru, name=''):
        self.obj = obj
        vectorextru.Normalize()
        self.vectorextru = vectorextru
        if obj.__class__ is Ellipse3D :
            self.points = obj.tessel_points
            # self.surface = obj.points
        else :
            self.points = obj.points
    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        if object_dict[arguments[1]].__class__ is Ellipse3D:
            ell = object_dict[arguments[1]]
            vectextru = -object_dict[arguments[2]]
            return cls(ell, vectextru, name)
        
        elif object_dict[arguments[1]].__class__ is BSplineCurve3D:
            bsplinecurve = object_dict[arguments[1]]
            vectextru = object_dict[arguments[2]]
            return cls(bsplinecurve, vectextru, name)
        else : 
            # surface = BSpline.Surface()
            # surface.degree_u = degree_u
            # surface.degree_v = degree_v
            # if weights is None:
            #     P = [(control_points[i][0], control_points[i][1], control_points[i][2]) for i in range(len(control_points))]
            #     surface.set_ctrlpts(P, nb_u, nb_v)
            # else:
            #     Pw = [(control_points[i][0]*weights[i], control_points[i][1]*weights[i], control_points[i][2]*weights[i], weights[i]) for i in range(len(control_points))]
            #     surface.set_ctrlpts(Pw, nb_u, nb_v)
            # knot_vector_u = []
            # for i, u_knot in enumerate(u_knots):
            #     knot_vector_u.extend([u_knot]*u_multiplicities[i])
            # knot_vector_v = []
            # for i, v_knot in enumerate(v_knots):
            #     knot_vector_v.extend([v_knot]*v_multiplicities[i])
            # surface.knotvector_u = knot_vector_u
            # surface.knotvector_v = knot_vector_v
            # surface.delta = 0.05
            # surface_points = surface.evalpts
            
            # self.surface = surface
            # self.points = [Point3D((p[0], p[1], p[2])) for p in surface_points]
            raise NotImplementedError  ## a adapter pour les bpsline
            

    # @classmethod
    # def from_step(cls, arguments, object_dict):
    #     name = arguments[0][1:-1]
        
        # degree_u = int(arguments[1])
        # degree_v = int(arguments[2])
        # points_sets = arguments[3][1:-1].split("),")
        # points_sets = [elem+")" for elem in points_sets[:-1]]+[points_sets[-1]]
        # control_points = []
        # for points_set in points_sets:
        #     points = [object_dict[int(i[1:])] for i in points_set[1:-1].split(",")]
        #     nb_v = len(points)
        #     control_points.extend(points)
        # nb_u = int(len(control_points) / nb_v)
        # surface_form = arguments[4]
        # if arguments[5] == '.F.':
        #     u_closed = False
        # elif arguments[5] == '.T.':
        #     u_closed = True
        # else:
        #     raise ValueError
        # if arguments[6] == '.F.':
        #     v_closed = False
        # elif arguments[6] == '.T.':
        #     v_closed = True
        # else:
        #     raise ValueError
        # self_intersect = arguments[7]
        # u_multiplicities = [int(i) for i in arguments[8][1:-1].split(",")]
        # v_multiplicities = [int(i) for i in arguments[9][1:-1].split(",")]
        # u_knots = [float(i) for i in arguments[10][1:-1].split(",")]
        # v_knots = [float(i) for i in arguments[11][1:-1].split(",")]
        # knot_spec = arguments[12]
    
        # if 13 in range(len(arguments)):
        #     weight_data = [float(i) for i in arguments[13][1:-1].replace("(", "").replace(")", "").split(",")]
        # else:
        #     weight_data = None
    
        # return cls(degree_u, degree_v, control_points, nb_u, nb_v, u_multiplicities, v_multiplicities, u_knots, v_knots, weight_data, name)
    
    
    
class CompositePrimitive3D(Primitive3D):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes  = ['basis_primitives']
    _non_eq_attributes = ['name', 'basis_primitives']
    _non_hash_attributes = []
    """
    A collection of simple primitives3D
    """
    def __init__(self, primitives, name=''):
        self.primitives = primitives
        basis_primitives=[]
        for primitive in primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.primitives)
            else:
                basis_primitives.append(primitive)

        Primitive3D.__init__(self, basis_primitives=basis_primitives, name=name)

    # def __eq__(self, other_):
    #     equal = True
    #     for primitive, other_primitive in zip(self.primitives, other_.primitives):
    #         equal = (equal and primitive == other_primitive)
    #     return equal

    def UpdateBasisPrimitives(self):
        # TODO: This is a copy/paste from CompositePrimitive2D, in the future make a Common abstract class
        basis_primitives=[]
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.primitives)
            else:
                basis_primitives.append(primitive)

        self.primitives = basis_primitives

    def MPLPlot(self, ax = None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        for primitive in self.edges:
            primitive.MPLPlot(ax)

        # ax.set_aspect('equal')

        return ax


class Wire3D(CompositePrimitive3D):
    """
    A collection of simple primitives, following each other making a wire
    """
    def __init__(self, primitives, name=''):
        CompositePrimitive3D.__init__(self, primitives, name)

    def Length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.Length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError


    # TODO: method to check if it is a wire
    def FreeCADExport(self, ip):
        name='primitive'+str(ip)

        s = 'E = []\n'
        for ip, primitive in enumerate(self.primitives):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'E.append(Part.Edge(L{}))\n'.format(ip)
        s += '{} = Part.Wire(E[:])\n'.format(name)

        return s

    def frame_mapping(self, frame, side, copy=True):
        new_wire = []
        if side == 'new':
            if copy :
                for primitive in self.primitives :
                    new_wire.append(primitive.frame_mapping(frame, side, copy))
                return Wire3D(new_wire)
            else :
                for primitive in self.primitives :
                    primitive.frame_mapping(frame, side, copy=False)
                
        if side == 'old':
            if copy :
                for primitive in self.primitives :
                    new_wire.append(primitive.frame_mapping(frame, side, copy))
                return Wire3D(new_wire)
            else :
                for primitive in self.primitives :
                    primitive.frame_mapping(frame, side, copy=False)
                
        
    def minimum_distance(self, wire2) :
        distance = []
        for element in self.primitives :
            for element2 in wire2.primitives :
                # fig = plt.figure()
                # ax = fig.add_subplot(111, projection='3d')
                # element.MPLPlot(ax=ax)
                # element2.MPLPlot(ax=ax)
                distance.append(element.minimum_distance(element2))
                
        return min(distance)
    
    def copy(self) :
        primitives_copy = []
        for primitive in self.primitives :
            primitives_copy.append(primitive.copy())
        return Wire3D(primitives_copy)

class Edge3D(Primitive3D):
    def __init__(self, edge_start, edge_end, name=''):
        Primitive3D.__init__(self, basis_primitives=[edge_start, edge_end], name=name)
        self.points = [edge_start, edge_end]

    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def __eq__(self, other_):
        equal = True
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point == other_point)
        return equal

    def to_dict(self):
        # improve the object structure ?
        dict_ = dc.DessiaObject.base_dict(self)
        dict_['edge_start'] = self.points[0]
        dict_['edge_end'] = self.points[1]
        return dict_

    @classmethod
    def from_step(cls, arguments, object_dict):
        if object_dict[arguments[3]].__class__ is Line3D:
            return LineSegment3D(object_dict[arguments[1]], object_dict[arguments[2]], arguments[0][1:-1])
        
        elif object_dict[arguments[3]].__class__ is Circle3D:
            # We supposed that STEP file is reading on trigo way
            center = object_dict[arguments[3]].center
            normal = object_dict[arguments[3]].normal
            normal.Normalize()
            radius = object_dict[arguments[3]].radius
            p1 = object_dict[arguments[1]]
            p2 = object_dict[arguments[2]]
            other_vec = object_dict[arguments[3]].other_vec
            if other_vec is None :
                other_vec = p1 - center
            other_vec.Normalize()
            frame = Frame3D(center, other_vec, normal.Cross(other_vec), normal)
            if p1 == p2:
                angle = math.pi
            else: 
                # p1_new, p2_new = frame.NewCoordinates(p1), frame.NewCoordinates(p2)
                # #Angle for p1
                # u1, u2 = p1_new.vector[0]/radius, p1_new.vector[1]/radius
                # theta1 = sin_cos_angle(u1, u2)
                # #Angle for p2
                # u3, u4 = p2_new.vector[0]/radius, p2_new.vector[1]/radius
                # theta2 = sin_cos_angle(u3, u4)
                theta1, theta2 = posangle_arc(p1, p2, radius, frame)
                if theta1 > theta2 : #sens trigo
                    angle = math.pi + (theta1 + theta2)/2
                else :
                    angle = (theta1 + theta2)/2
            p_3 = Point3D((radius*math.cos(angle), radius*math.sin(angle),0))
            p3 = frame.OldCoordinates(p_3)
            if p1==p3 or p2==p3 :
                p_3 = Point3D((radius*math.cos(0), radius*math.sin(0),0))
                p3 = frame.OldCoordinates(p_3)
            arc = Arc3D(p1, p3, p2, normal, arguments[0][1:-1], other_vec)
            if math.isclose(arc.radius, 0, abs_tol=1e-9):
                if p1 == p2 :
                    p_3 = Point3D((radius*math.cos(0), radius*math.sin(0),0))
                    p3 = frame.OldCoordinates(p_3)
                    arc = Arc3D(p1, p3, p2, normal, arguments[0][1:-1], other_vec)
            return arc

        elif object_dict[arguments[3]].__class__ is Ellipse3D:
            majorax = object_dict[arguments[3]].major_axis
            minorax = object_dict[arguments[3]].minor_axis
            center = object_dict[arguments[3]].center
            normal = object_dict[arguments[3]].normal
            normal.Normalize()
            majordir = object_dict[arguments[3]].major_dir
            majordir.Normalize()
            minordir = normal.Cross(majordir)
            minordir.Normalize()
            frame = Frame3D(center, majordir, minordir, normal)
            p1 = object_dict[arguments[1]] #on part du principe que p1 suivant majordir
            p2 = object_dict[arguments[2]]
            if p1 == p2: 
                angle = 5*math.pi/4
                xtra = Point3D((majorax*math.cos(math.pi/2), minorax*math.sin(math.pi/2),0))
                extra = frame.OldCoordinates(xtra) 
            else :
                extra = None
                ## Positionnement des points dans leur frame
                p1_new, p2_new = frame.NewCoordinates(p1), frame.NewCoordinates(p2)
                #Angle pour le p1
                u1, u2 = p1_new.vector[0]/majorax, p1_new.vector[1]/minorax
                theta1 = sin_cos_angle(u1, u2)
                #Angle pour le p2
                u3, u4 = p2_new.vector[0]/majorax, p2_new.vector[1]/minorax
                theta2 = sin_cos_angle(u3, u4)
                
                if theta1 > theta2 : #sens trigo
                    angle = math.pi + (theta1 + theta2)/2
                else :
                    angle = (theta1 + theta2)/2
                
            p_3 = Point3D((majorax*math.cos(angle), minorax*math.sin(angle),0))
            p3 = frame.OldCoordinates(p_3)
            
            arcellipse = ArcEllipse3D(p1, p3, p2, center, majordir, normal, arguments[0][1:-1], extra)
            
            return arcellipse
        
        elif object_dict[arguments[3]].__class__ is BSplineCurve3D:
            # print(object_dict[arguments[1]], object_dict[arguments[2]])
            # BSplineCurve3D à couper à gauche et à droite avec les points ci dessus ?
            return object_dict[arguments[3]]

        else:
            print(object_dict[arguments[3]])
            raise NotImplementedError
        
    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_edge_start = self.points[0].Rotation(center, axis, angle, copy=True)
            new_edge_end = self.points[1].Rotation(center, axis, angle, copy=True)
            return Edge3D(new_edge_start, new_edge_end)
        else:
            self.points[0].Rotation(center, axis, angle, copy=False)
            self.points[1].Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_edge_start = self.points[0].Translation(offset, copy=True)
            new_edge_end = self.points[1].Translation(offset, copy=True)
            return Edge3D(new_edge_start, new_edge_end)
        else:
            self.points[0].Translation(offset, copy=False)
            self.points[1].Translation(offset, copy=False)


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_edge_start = self.points[0].frame_mapping(frame, side, copy=True)
            new_edge_end = self.points[1].frame_mapping(frame, side, copy=True)
            return Edge3D(new_edge_start, new_edge_end)
        else:
            # self.points[0].frame_mapping(frame, side, copy=False)
            # self.points[1].frame_mapping(frame, side, copy=False)
            self.points[0] = self.points[0].frame_mapping(frame, side, copy=True)
            self.points[1] = self.points[1].frame_mapping(frame, side, copy=True)

    def copy(self):
        new_edge_start = self.points[0].copy()
        new_edge_end = self.points[1].copy()
        return Edge3D(new_edge_start, new_edge_end)


class LineSegment3D(Edge3D):
    """
    Define a line segment limited by two points
    """
    def __init__(self, point1, point2, name=''):
        Edge3D.__init__(self, point1, point2, name='')
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return hash(self.points[0]) + hash(self.points[1])

    def to_dict(self):
        # improve the object structure ?
        dict_ = dc.DessiaObject.base_dict(self)
        dict_['point1'] = self.points[0].to_dict()
        dict_['point2'] = self.points[1].to_dict()
        dict_['object_class'] = 'volmdlr.core.LineSegment3D'
        return dict_
    
    @classmethod
    def dict_to_object(cls, dict_):
        return cls(point1 = Point3D.dict_to_object(dict_['point1']), 
                   point2 = Point3D.dict_to_object(dict_['point2']), name = dict_['name'])
                       

    def _bounding_box(self):
        points = self.points

        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])

        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def DirectionVector(self, unit=False):
        u = self.points[1] - self.points[0]
        if unit:
            u.Normalize()
        return u


    def Length(self):
        return self.points[1].point_distance(self.points[0])
    
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + curvilinear_abscissa*(self.points[1] - self.points[0]) / self.Length()

    def middle_point(self):
        l = self.Length()
        return self.PointAtCurvilinearAbscissa(0.5*l)

    def PlaneProjection2D(self, center, x, y):
        return LineSegment2D(self.points[0].PlaneProjection2D(center, x, y),
                             self.points[1].PlaneProjection2D(center, x, y))

    def Intersection(self, segment2):
        x1 = self.points[0].vector[0]
        y1 = self.points[0].vector[1]
        z1 = self.points[0].vector[2]
        x2 = self.points[1].vector[0]
        y2 = self.points[1].vector[1]
        z2 = self.points[1].vector[2]
        x3 = segment2.points[0].vector[0]
        y3 = segment2.points[0].vector[1]
        z3 = segment2.points[0].vector[2]
        x4 = segment2.points[1].vector[0]
        y4 = segment2.points[1].vector[1]
        z4 = segment2.points[1].vector[2]
        
        if x3 == 0 and x4 ==0 and y4-y3 == 0 :
            x5, y5, z5 = x3, y3, z3
            x6, y6, z6 = x4, y4, z4
            x3, y3, z3 = x1, y1, z1
            x4, y4, z4 = x2, y2, z2
            x1, y1, z1 = x5, y5, z5
            x2, y2, z2 = x6, y6, z6
        
        elif y3 == 0 and y4 ==0 and x4-x3 == 0 :
            x5, y5, z5 = x3, y3, z3
            x6, y6, z6 = x4, y4, z4
            x3, y3, z3 = x1, y1, z1
            x4, y4, z4 = x2, y2, z2
            x1, y1, z1 = x5, y5, z5
            x2, y2, z2 = x6, y6, z6
        
        res, list_t1 = [], []

        #2 unknown 3eq with t1 et t2 unknown
        if (x2-x1+y1-y2) != 0 and (y4-y3) != 0:
            t1 = (x3-x1 + (x4-x3)*(y1-y3)/(y4-y3))/(x2-x1+y1-y2)
            t2 = (y1-y3 + (y2-y1)*t1)/(y4-y3)
            res1 = z1 + (z2-z1)*t1
            res2 = z3 + (z4-z3)*t2
            list_t1.append(t1)
            res.append([res1, res2])
        
        if (z2-z1+y1-y2) != 0 and (y4-y3) != 0:
            t1 = (z3-z1 + (z4-z3)*(y1-y3)/(y4-y3))/(z2-z1+y1-y2)
            t2 = (y1-y3 + (y2-y1)*t1)/(y4-y3)
            res1 = x1 + (x2-x1)*t1
            res2 = x3 + (x4-x3)*t2
            list_t1.append(t1)
            res.append([res1, res2])
        
        if (z2-z1+x1-x2) !=0 and (x4-x3) != 0 :
            t1 = (z3-z1 + (z4-z3)*(x1-x3)/(x4-x3))/(z2-z1+x1-x2)
            t2 = (x1-x3 + (x2-x1)*t1)/(x4-x3)
            res1 = y1 + (y2-y1)*t1
            res2 = y3 + (y4-y3)*t2
            list_t1.append(t1)
            res.append([res1, res2])
        
        if len(res)==0 :
            return None
        
        for pair, t1 in zip(res, list_t1) :
            res1, res2 = pair[0], pair[1]
            if math.isclose(res1, res2, abs_tol=1e-7) : #if there is an intersection point
                if t1>=0 or t1<=1:
                    return Point3D([x1+(x2-x1)*t1, y1+(y2-y1)*t1, z1+(z2-z1)*t1])
        
        return None

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            return LineSegment3D(*[p.Rotation(center, axis, angle, copy=True) for p in self.points])
        else:
            Edge3D.Rotation(self, center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def __contains__(self, point):
        point1, point2 = self.points[0], self.points[1]
        axis = Vector3D(point2 - point1)
        test = point.Rotation(point1, axis, math.pi)
        if test==point :
            return True
        else : 
            return False
    
    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment3D(*[p.Translation(offset, copy=True) for p in self.points])
        else:
            Edge3D.Translation(self, offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return LineSegment3D(*[frame.OldCoordinates(p) for p in self.points])
            else:
                Edge3D.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()
        if side == 'new':
            if copy:
                return LineSegment3D(*[frame.NewCoordinates(p) for p in self.points])
            else:
                Edge3D.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()

    def copy(self):
        return LineSegment3D(*[p.copy() for p in self.points])

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        x=[p.vector[0] for p in self.points]
        y=[p.vector[1] for p in self.points]
        z=[p.vector[2] for p in self.points]
        ax.plot(x,y,z, 'o-k')
        return ax

    def MPLPlot2D(self, x_3D, y_3D, ax=None, color='k', width=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        edge2D =  self.PlaneProjection2D(O3D, x_3D, y_3D)
        edge2D.MPLPlot(ax=ax, color=color, width=width)
        return ax

    def plot_data(self, x_3D, y_3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        edge2D =  self.PlaneProjection2D(O3D, x_3D, y_3D)
        return edge2D.plot_data(marker, color, stroke_width,
                         dash, opacity, arrow)

    def FreeCADExport(self, name, ndigits=6):
        name = 'primitive'+str(name)
        x1, y1, z1 = round(1000*self.points[0], ndigits).vector
        x2, y2, z2 = round(1000*self.points[1], ndigits).vector
        print('name', name)
        return '{} = Part.LineSegment(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,x1,y1,z1,x2,y2,z2)

    def to_line(self):
        return Line3D(*self.points)

    def babylon_script(self, color=(1, 1, 1), name='line',  type_='line', parent=None):
        if type_ == 'line' or type_ == 'dashed':
            s = 'var myPoints = [];\n'
            s += 'var point1 = new BABYLON.Vector3({},{},{});\n'.format(*self.points[0])
            s += 'myPoints.push(point1);\n'
            s += 'var point2 = new BABYLON.Vector3({},{},{});\n'.format(*self.points[1])
            s += 'myPoints.push(point2);\n'
            if type_ == 'line':
                s += 'var {} = BABYLON.MeshBuilder.CreateLines("lines", {{points: myPoints}}, scene);\n'.format(name)
            elif type_ == 'dashed':
                s += 'var {} = BABYLON.MeshBuilder.CreateDashedLines("lines", {{points: myPoints, dashNb:20}}, scene);'.format(name)
            s += '{}.color = new BABYLON.Color3{};\n'.format(name, tuple(color))
        elif type_ == 'tube':
            radius = 0.03*self.points[0].point_distance(self.points[1])
            s = 'var points = [new BABYLON.Vector3({},{},{}), new BABYLON.Vector3({},{},{})];\n'.format(*self.points[0], *self.points[1])
            s += 'var {} = BABYLON.MeshBuilder.CreateTube("frame_U", {{path: points, radius: {}}}, {});'.format(name, radius, parent)
#            s += 'line.material = red_material;\n'

        else:
            raise NotImplementedError

        if parent is not None:
            s += '{}.parent = {};\n'.format(name, parent)

        return s

    def To2D(self, plane_origin, x1, x2):
        p2D=[p.To2D(plane_origin,x1,x2) for p in self.points]
        return LineSegment2D(*p2D,name=self.name)

    def reverse(self):
        return LineSegment3D(self.points[1].copy(), self.points[0].copy())

    def MinimumDistancePoints(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.points[1] - self.points[0]
        v = other_line.points[1] - other_line.points[0]
        w = self.points[0] - other_line.points[0]
        a = u.Dot(u)
        b = u.Dot(v)
        c = v.Dot(v)
        d = u.Dot(w)
        e = v.Dot(w)
        if (a*c - b**2)!=0 :
            s = (b*e -c*d) / (a*c - b**2)
            t = (a*e -b*d) / (a*c - b**2)
            p1 = self.points[0] + s*u
            p2 = other_line.points[0] + t*v
            return p1, p2, s, t
        else :
            return None, None, -1, -1

    def Matrix_distance(self, other_line) :
        u = self.DirectionVector()
        v = other_line.DirectionVector()
        w = other_line.points[0] - self.points[0]
        
        a = u.Dot(u)
        b = -u.Dot(v)
        d = v.Dot(v)
        
        e = w.Dot(u)
        f = -w.Dot(v)
        
        A = npy.array([[a, b],
                       [b, d]])
        B = npy.array([e, f])
        
        res = scp.optimize.lsq_linear(A, B, bounds=(0,1))
        p1 = self.PointAtCurvilinearAbscissa(res.x[0]*self.Length())
        p2 = other_line.PointAtCurvilinearAbscissa(res.x[1]*other_line.Length())
        return p1, p2
    
    def parallele_distance(self, LS2):
        ptA, ptB, ptC = self.points[0], self.points[1], LS2.points[0]
        u = Vector3D((ptA - ptB).vector)
        u.Normalize()
        plane1 = Plane3D.from_3_points(ptA, ptB, ptC)
        v = u.Cross(plane1.normal) #distance vector
        #ptA = k*u + c*v + ptC
        res = (ptA - ptC).vector
        x, y, z = res[0], res[1], res[2]
        u1, u2, u3 = u.vector[0], u.vector[1], u.vector[2]
        v1, v2, v3 = v.vector[0], v.vector[1], v.vector[2]
        
        if (u1*v2-v1*u2)!=0 and u1!=0: 
            c = (y*u1-x*u2)/(u1*v2-v1*u2)
            k = (x-c*v1)/u1
            if math.isclose(k*u3+c*v3, z, abs_tol=1e-7) :
                return k
        elif (u1*v3-v1*u3)!=0 and u1!=0: 
            c = (z*u1-x*u3)/(u1*v3-v1*u3)
            k = (x-c*v1)/u1
            if math.isclose(k*u2+c*v2, y, abs_tol=1e-7) :
                return k
        elif (v1*u2-v2*u1)!=0 and u2!=0:
            c = (u2*x-y*u1)/(v1*u2-v2*u1)
            k = (y-c*v2)/u2
            if math.isclose(k*u3+c*v3, z, abs_tol=1e-7) :
                return k
        elif (v3*u2-v2*u3)!=0 and u2!=0:
            c = (u2*z-y*u3)/(v3*u2-v2*u3)
            k = (y-c*v2)/u2
            if math.isclose(k*u1+c*v1, x, abs_tol=1e-7) :
                return k
        elif (u1*v3-v1*u3)!=0 and u3!=0: 
            c = (z*u1-x*u3)/(u1*v3-v1*u3)
            k = (z-c*v3)/u3
            if math.isclose(k*u2+c*v2, y, abs_tol=1e-7) :
                return k
        elif (u2*v3-v2*u3)!=0 and u3!=0: 
            c = (z*u2-y*u3)/(u2*v3-v2*u3)
            k = (z-c*v3)/u3
            if math.isclose(k*u1+c*v1, x, abs_tol=1e-7) :
                return k
        else :
            return NotImplementedError
        
    def minimum_distance(self, element, return_points=False):
        if element.__class__ is Arc3D or element.__class__ is Circle3D:
            pt1, pt2 = element.minimum_distance_points_line(self)
            if return_points : 
                return pt1.point_distance(pt2), pt1, pt2
            else :
                return pt1.point_distance(pt2)
        
        elif element.__class__ is LineSegment3D :
            p1, p2 = self.Matrix_distance(element)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)

        else :
            return NotImplementedError
            
        

class Contour3D(Wire3D):
    _non_serializable_attributes = ['points']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['points', 'name']
    _generic_eq = True
    """
    A collection of 3D primitives forming a closed wire3D
    """
    def __init__(self, edges, point_inside_contour=None, name=''):
        # TODO: docstring in english
        """
        Faire un choix : soit edges c'est un CompositePrimitives3D
        ou alors un ensemble de primitives
        ou alors un ensemble de basis_primtives (qui sont des points pour le moment)
        """

        self.name = name
        self.point_inside_contour = point_inside_contour

        edges_primitives = []
        for edge in edges:
            if edge.__class__ == CompositePrimitive3D:
                edges_primitives.extend(edge.primitives)
            else:
                edges_primitives.append(edge)
        self.edges = edges_primitives

        if self.edges[0].__class__.__name__ == 'Contour3D':
            raise ValueError

        self.tessel_points = self.clean_points()

    def __hash__(self):
        return sum([hash(e) for e in self.edges]) + sum([hash(p) for p in self.tessel_points])

    def __eq__(self, other_):
        equal = True
        for edge, other_edge in zip(self.edges, other_.edges):
            equal = (equal and edge == other_edge)
        # for point, other_point in zip(self.points, other_.points):
        #     equal = (equal and point == other_point)
        #     print('contour', equal, point.vector, other_point.vector)
        return equal

    @classmethod
    def from_step(cls, arguments, object_dict):
        edges = []
        for edge in arguments[1]:
            # print(arguments[1])
            edges.append(object_dict[int(edge[1:])])
            # print(edges)
#        points = edges[0].points[:]
#        for i, edge in enumerate(edges[1:-1]):
#            if edge.points[0] in points[-2:]:
#                points.append(edge.points[1])
#            elif edge.points[1] in points[-2:]:
#                points.append(edge.points[0])
#            else:
#                raise NotImplementedError
#        contour_points = [p.copy() for p in points]
        ################################################################ print('arg', arguments)

        return cls(edges, point_inside_contour=None, name=arguments[0][1:-1])

    def clean_points(self):
        """
        TODO : verifier si le dernier point est toujours le meme que le premier point
        lors d'un import step par exemple
        """
        # print('!' , self.edges)
        if hasattr(self.edges[0], 'points'):
            points = self.edges[0].points[:]
        else:
            points = self.edges[0].tessellation_points()
        for edge in self.edges[1:]:
            if hasattr(edge, 'points'):
                points_to_add = edge.points[:]
            else :
                points_to_add = edge.tessellation_points()
            if points[0] == points[-1]: # Dans le cas où le (dernier) edge relie deux fois le même point
                points.extend(points_to_add[::-1])
            
            elif points_to_add[0] == points[-1]:
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[-1]:
                points.extend(points_to_add[-2::-1])
            elif points_to_add[0] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[-2::-1])
            else:
                d1, d2 = (points_to_add[0]-points[0]).Norm(), (points_to_add[0]-points[-1]).Norm()
                d3, d4 = (points_to_add[-1]-points[0]).Norm(), (points_to_add[-1]-points[-1]).Norm()
                if math.isclose(d2, 0, abs_tol=1e-3):
                    points.extend(points_to_add[1:])
                elif math.isclose(d4, 0, abs_tol=1e-3):
                    points.extend(points_to_add[-2::-1])
                elif math.isclose(d1, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[1:])
                elif math.isclose(d3, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[-2::-1])
                
        if len(points) > 1:
            if points[0] == points[-1]:
                points.pop()
                
        return points

    def average_center_point(self):
        nb = len(self.tessel_points)
        x = npy.sum([p[0] for p in self.points]) / nb
        y = npy.sum([p[1] for p in self.points]) / nb
        z = npy.sum([p[2] for p in self.points]) / nb
        return Point3D((x,y,z))

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_edges = [edge.Rotation(center, axis, angle, copy=True) for edge in self.edges]
            # new_points = [p.Rotation(center, axis, copy=True) for p in self.points]
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.edges:
                edge.Rotation(center, axis, angle, copy=False)
            for point in self.tessel_points:
                point.Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_edges = [edge.Translation(offset, copy=True) for edge in self.edges]
            # new_points = [p.Translation(offset, copy=True) for p in self.points]
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.edges:
                edge.Translation(offset, copy=False)
            for point in self.tessel_points:
                point.Translation(offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_edges = [edge.frame_mapping(frame, side, copy=True) for edge in self.edges]
            # new_points = [p.frame_mapping(frame, side, copy=True) for p in self.points]
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.edges:
                edge.frame_mapping(frame, side, copy=False)
            for point in self.tessel_points:
                point.frame_mapping(frame, side, copy=False)

    def copy(self):
        new_edges = [edge.copy() for edge in self.edges]
        if self.point_inside_contour is not None:
            new_point_inside_contour = self.point_inside_contour.copy()
        else:
            new_point_inside_contour = None
        return Contour3D(new_edges, new_point_inside_contour, self.name)
    
    def Length(self):
        # TODO: this is duplicated code from Wire3D!
        length = 0.
        for edge in self.edges:
            length += edge.Length()
        return length
    
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        # TODO: this is duplicated code from Wire3D!
        length = 0.
        for primitive in self.edges:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError
        
    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None
            
        for edge in self.edges:
            edge.MPLPlot(ax=ax)
            
        return ax
            

class Circle3D(Contour3D):
    _non_serializable_attributes = ['point', 'edges', 'point_inside_contour']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    
    def __init__(self, center, radius, normal, name='', other_vec=None):
        self.center = center
        self.radius = radius
        self.normal = normal
        self.angle = 2*math.pi
        self.points = self.tessellation_points()
        self.other_vec = other_vec
        Contour3D.__init__(self, [self], name=name)
        
    def __hash__ (self):
        return int(round(1e6*(self.center.vector[0] + self.center.vector[1] + self.center.vector[2] + self.radius)))

    def __eq__(self, other_circle):
        return math.isclose(self.center.vector[0], other_circle.center.vector[0], abs_tol=1e-06) \
           and math.isclose(self.center.vector[1], other_circle.center.vector[1], abs_tol=1e-06) \
           and math.isclose(self.center.vector[2], other_circle.center.vector[2], abs_tol=1e-06) \
           and math.isclose(self.radius, other_circle.radius, abs_tol=1e-06)
        
#    def _get_points(self):
#        vr = Vector3D(npy.random.random(3))
#        vr.Normalize()
#        vn = vr.Cross(self.normal)
#        dir_radius = vn*self.radius
#        pt1 = self.center + dir_radius
#        pt2 = pt1.Rotation(self.center, self.normal, 2*math.pi/3.)
#        pt3 = pt2.Rotation(self.center, self.normal, 2*math.pi/3.)
#        return [pt1, pt2, pt3]
#
#    points=property(_get_points)

    def tessellation_points(self, resolution=20):
        plane = Plane3D.from_normal(self.center, self.normal)
        tessellation_points_3D = [self.center + self.radius*math.cos(teta)*plane.vectors[0] + self.radius*math.sin(teta)*plane.vectors[1] \
            for teta in npy.linspace(0, 2*math.pi, resolution+1)][:-1]
        return tessellation_points_3D
    
    def Length(self):
        return 2* math.pi * self.radius

    def FreeCADExport(self, name, ndigits=3):
#        name = 'primitive{}'.format(ip)
        xc,yc,zc = round(1000*self.center, ndigits)
        xn,yn,zn = round(self.normal, ndigits)
        return '{} = Part.Circle(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(name,xc,yc,zc,xn,yn,zn,1000*self.radius)

    def Rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.Rotation(rot_center, axis, angle, True)
        new_normal = self.normal.Rotation(rot_center, axis, angle, True)
        if copy:
            return Circle3D(new_center, self.radius, new_normal, self.name)
        else:
            self.center = new_center
            self.normal = new_normal

    def Translation(self, offset, copy=True):
        new_center = self.center.Translation(offset, True)
        new_normal = self.normal.Translation(offset, True)
        if copy:
            return Circle3D(new_center, self.radius, new_normal, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            
    def MPLPlot(self, ax=None, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        x = []
        y = []
        z = []
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color)
        return ax
            

    @classmethod
    def from_step(cls, arguments, object_dict):
        center = object_dict[arguments[1]].origin
        radius = float(arguments[2])/1000
        if object_dict[arguments[1]].u is not None:
            normal = object_dict[arguments[1]].u
            other_vec = object_dict[arguments[1]].v
            if other_vec is not None :
                other_vec.Normalize()
        else:
            normal = object_dict[arguments[1]].v ### ou w 
            other_vec = None
        normal.Normalize()
        
        return cls(center, radius, normal, arguments[0][1:-1], other_vec)
    
    def To2D(self, plane_origin, x, y):
        pc = self.center.To2D(plane_origin, x, y)
        return Circle2D(pc, self.radius, self.name)

class Ellipse3D(Contour3D):
    """
    :param major_axis: Largest radius of the ellipse
    :type major_axis: float
    :param minor_axis: Smallest radius of the ellipse
    :type minor_axis: float
    :param center: Ellipse's center
    :type center: Point3D
    :param normal: Ellipse's normal
    :type normal: Vector3D
    :param major_dir: Direction of the largest radius/major_axis
    :type major_dir: Vector3D
    """
    def __init__(self, major_axis, minor_axis, center, normal, major_dir, name=''):
        
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.center = center
        normal.Normalize()
        self.normal = normal
        major_dir.Normalize()
        self.major_dir = major_dir
        Contour3D.__init__(self, [self], name=name)
        
    def tessellation_points(self, resolution=20):
        # plane = Plane3D.from_normal(self.center, self.normal)
        tessellation_points_3D = [self.center + self.major_axis*math.cos(teta)*self.major_dir + self.minor_axis*math.sin(teta)*self.major_dir.Cross(self.normal) \
                                  for teta in npy.linspace(0, 2*math.pi, resolution+1)][:-1]
        return tessellation_points_3D

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        xc, yc, zc = npy.round(1000*self.center.vector, ndigits)
        major_vector = self.center + self.major_axis/2 * self.major_dir
        xmaj, ymaj, zmaj = npy.round(1000*major_vector.vector, ndigits)
        minor_vector = self.center + self.minor_axis/2 * self.normal.Cross(self.major_dir)
        xmin, ymin, zmin = npy.round(1000*minor_vector.vector, ndigits)
        return '{} = Part.Ellipse(fc.Vector({},{},{}), fc.Vector({},{},{}), fc.Vector({},{},{}))\n'.format(name,xmaj,ymaj,zmaj,xmin,ymin,zmin,xc,yc,zc)

    def Rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.Rotation(rot_center, axis, angle, True)
        new_normal = self.normal.Rotation(rot_center, axis, angle, True)
        new_major_dir = self.major_dir.Rotation(rot_center, axis, angle, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center, new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir

    def Translation(self, offset, copy=True):
        new_center = self.center.Translation(offset, True)
        new_normal = self.normal.Translation(offset, True)
        new_major_dir = self.major_dir.Translation(offset, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center, new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir
            
    def MPLPlot(self, ax=None, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        x = []
        y = []
        z = []
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color)
        return ax

    @classmethod
    def from_step(cls, arguments, object_dict):
        center = object_dict[arguments[1]].origin
        normal = object_dict[arguments[1]].u   #ancien w
        major_dir = object_dict[arguments[1]].v # ancien u
        major_axis = float(arguments[2])/1000
        minor_axis = float(arguments[3])/1000
        return cls(major_axis, minor_axis, center, normal, major_dir, arguments[0][1:-1])

class Face3D(Primitive3D):
    def __init__(self, contours, radius=None):
        self.contours3d = contours
        self.bounding_box = self._bounding_box()
        
    def _bounding_box(self):
        points = self.contours3d[0].tessel_points

        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])

        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)
    
    @classmethod
    def from_step(cls, arguments, object_dict):
        contours = []
        contours.append(object_dict[int(arguments[1][0][1:])])
        
        # print(object_dict[int(arguments[2])].__class__)
        # plane = Plane3D.from_points(contours[0].points)
        # contours[0].points, polygon2D = cls._repair_points_and_polygon2d(contours[0].points, plane)
        # points = [p.copy() for p in contours[0].points[:]]
        
        # if int(arguments[2]) not in list(object_dict.keys()):
            # print(int(arguments[2]))
#            return None
        if object_dict[int(arguments[2])].__class__  is Plane3D:
            # [print(e.points) for e in contours]
            # print(arguments)
            # object_dict[110].MPLPlot()
            
            # return PlaneFace3D(contours, plane=None, points=None, polygon2D=None, name=arguments[0][1:-1])
            return PlaneFace3D.from_contours3d(contours, name=arguments[0][1:-1])
    
        elif object_dict[int(arguments[2])].__class__  is CylindricalSurface3D:
            return CylindricalFace3D.from_contour3d(contours, object_dict[int(arguments[2])], name=arguments[0][1:-1])

        elif object_dict[int(arguments[2])].__class__  is BSplineExtrusion:
            return BSplineFace3D(contours, object_dict[int(arguments[2])], name=arguments[0][1:-1])
        
        elif object_dict[int(arguments[2])].__class__  is BSplineSurface3D:
            # print(object_dict[int(arguments[2])])
            return BSplineFace3D(contours, object_dict[int(arguments[2])], name=arguments[0][1:-1])
        
        elif object_dict[int(arguments[2])].__class__  is ToroidalSurface3D:
            return ToroidalFace3D.from_contour3d(contours, object_dict[int(arguments[2])], name=arguments[0][1:-1])
        
        elif object_dict[int(arguments[2])].__class__  is ConicalSurface3D:
            return ConicalFace3D.from_contour3d(contours, object_dict[int(arguments[2])], name=arguments[0][1:-1])
        
        elif object_dict[int(arguments[2])].__class__  is SphericalSurface3D:
            return SphericalFace3D.from_contour3d(contours, object_dict[int(arguments[2])], name=arguments[0][1:-1])
        
        else:
            print('arguments', arguments)
            print(object_dict[int(arguments[2])])
            raise NotImplementedError
    
    def delete_double(self, Le):
            Ls = []
            for i in Le:
                if i not in Ls:
                    Ls.append(i)
            return Ls
        
    def min_max(self, Le, pos):
        Ls = []
        for i in range (0,len(Le)) :
            Ls.append(Le[i][pos])
        return (min(Ls), max(Ls))
    
    def range_trigo(list_point) :
        points_set = delete_double_point(list_point)
        xmax, xmin = max(pt[0] for pt in points_set), min(pt[0] for pt in points_set)
        ymax, ymin = max(pt[1] for pt in points_set), min(pt[1] for pt in points_set) 
        center = Point2D(((xmax+xmin)/2, (ymax+ymin)/2))
        frame2d = Frame2D(center, X2D, Y2D)
        points_test = [frame2d.NewCoordinates(pt) for pt in points_set]
        
        points_2dint = []
        s = 0
        for k in range(0, len(points_test)) :
            closest = points_test[s]
            while closest is None : 
                s += 1
                closest = points_test[s]
            angle_min = math.atan2(closest.vector[1], closest.vector[0]) + math.pi
            pos = s
            for i in range (s+1, len(points_test)) :
                close_test = points_test[i]
                if close_test is None : 
                    continue
                else :
                    angle_test = math.atan2(close_test.vector[1], close_test.vector[0]) + math.pi
                    if angle_test < angle_min : #and dist_test <= dist_min:
                        angle_min = angle_test
                        closest = close_test
                        pos = i
            points_2dint.append(closest)
            points_test[pos] = None
        
        points_old = [frame2d.OldCoordinates(pt) for pt in points_2dint]
        return points_old
    
    def range_closest(list_point, r1=None, r2=None) :
        #use r1, r2 to compare h and r1*angle or r1*angle and r2*angle
        points_set = delete_double_point(list_point)
        if r1 is not None :
            for k in range(0, len(points_set)) :
                points_set[k].vector[0] = points_set[k].vector[0]*r1
        if r2 is not None :
            for k in range(0, len(points_set)) :
                points_set[k].vector[1] = points_set[k].vector[1]*r2
        
        points_2dint = [points_set[0]]
        s = 1
        for k in range(1, len(points_set)) :
            closest = points_set[s]
            while closest is None : 
                s += 1
                closest = points_set[s]
            dist_min = (points_2dint[-1] - closest).Norm()
            pos = s
            for i in range (s+1, len(points_set)) :
                close_test = points_set[i]
                if close_test is None : 
                    continue
                else :
                    dist_test = (points_2dint[-1] - close_test).Norm()
                    if dist_test <= dist_min:
                        dist_min = dist_test
                        closest = close_test
                        pos = i
            points_2dint.append(closest)
            points_set[pos] = None
        
        if r1 is not None :
            for k in range(0, len(points_2dint)) :
                points_2dint[k].vector[0] = points_2dint[k].vector[0]/r1
        if r2 is not None :
            for k in range(0, len(points_2dint)) :
                points_2dint[k].vector[1] = points_2dint[k].vector[1]/r2
        
        return points_2dint
    
    def cut_contours(self, contours2d, resolution) :
        placement2d = contours2d[0].tessel_points 
        segmentspt = []
        offset = 0.0001
        xmin, xmax = min([pt[0] for pt in placement2d]), max([pt[0] for pt in placement2d])
        ymin, ymax = min([pt[1] for pt in placement2d]), max([pt[1] for pt in placement2d])
        for k in range (0,len(placement2d)):
            if k == len(placement2d)-1 :
                pt1, pt2 = segmentspt[-1].points[1], segmentspt[0].points[0]
                segmentspt.append(LineSegment2D(pt1,pt2))
            else : 
                pt1, pt2 = placement2d[k], placement2d[k+1]
                segpoints = []
                if pt1[0] == xmin :
                    segpoints.append(pt1 - Point2D((offset,0)))
                elif pt1 [0] == xmax :
                    segpoints.append(pt1 + Point2D((offset,0)))
                else : 
                    segpoints.append(pt1)
                if pt2[0] == xmin :
                    segpoints.append(pt2 - Point2D((offset,0)))
                elif pt2 [0] == xmax :
                    segpoints.append(pt2 + Point2D((offset,0)))
                else : 
                    segpoints.append(pt2)
                segmentspt.append(LineSegment2D(segpoints[0],segpoints[1]))
        
        pas = (xmax-xmin)/resolution
        pointsxzmin = [Point2D([xmin+i*pas,ymin-1]) for i in range (0,resolution+1)]
        pointsxzmax = [Point2D([xmin+i*pas,ymax+1]) for i in range (0,resolution+1)]
        
        line = [LineSegment2D(ptxmax, ptxmin) for ptxmin,ptxmax in list(zip(pointsxzmin, pointsxzmax))]

        all_contours_points, contours_points = [], []
        register = []
        for l in line :
            nb_pt = 0
            for s in segmentspt :
                p = l.line_intersection(s)
                if p is None :
                    continue
                else :
                    nb_pt += 1
                    contours_points.append(p)
            register.append(nb_pt)
        
        seuil = 0
        for k in range(0, len(register)):
            list_pt = []
            for i in range(seuil, seuil+register[k]):
                list_pt.append(contours_points[i])
            seuil += register[k]
            if k == len(register)-1:
                for i in range(0, register[0]):
                    list_pt.append(contours_points[register[0]-1-i])
            else :
                for i in range(seuil+register[k+1]-1, seuil-1, -1):
                    list_pt.append(contours_points[i])
            all_contours_points.append(list_pt)

        all_contours_points.pop() 

        return all_contours_points

    def create_primitives(points) :
        primitives = []
        for k in range(0, len(points)) :
            if k == len(points)-1 :
                primitives.append(LineSegment2D(points[k], points[0]))
            else :
                primitives.append(LineSegment2D(points[k], points[k+1]))
        return primitives
    
    def LS2D_inprimitives(ls_toadd, primitives) :
        same = False
        for list_prim in primitives :
            for prim in list_prim :
                if ls_toadd.points[0] == prim.points[0] and ls_toadd.points[-1] == prim.points[-1] :
                    same = True
                elif ls_toadd.points[0] == prim.points[-1] and ls_toadd.points[-1] == prim.points[0] :
                    same = True
                else : 
                    continue
        return same

class PlaneFace3D(Face3D):
    """
    :param contours: The face's contour2D 
    :type contours: Contour2D
    :param plane: Plane used to place your face
    :type plane: Plane3D
    """
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes  = ['bounding_box', 'polygon2D']
    _non_eq_attributes = ['name', 'bounding_box']
    _non_hash_attributes = []

    def __init__(self, contours, plane, points=None, polygon2D=None, name=''):
        if contours[0].__class__ is Contour3D :
            raise ValueError('You must use Contour2D or use from_contours3d')
            
        self.name = name
        self.contours = contours
        self.plane = plane
        self.setup_planeface(self.contours, self.plane, points=points, polygon2D=polygon2D, name=self.name)

    def setup_planeface(self, ctrs2d, plane, points=None, polygon2D=None, name=''):
        if points is None or polygon2D is None:
            self.points = self.contours[0].tessel_points
            self.polygon2D = Polygon2D(self.points)
        else :
            self.points = points
            self.polygon2D = polygon2D
        ctr3d = ctrs2d[0].copy()
        Face3D.__init__(self, [ctr3d.To3D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1])])
            

    def __hash__(self):
        return hash(self.plane) + sum([hash(p) for p in self.points])

    def __eq__(self, other_):
        equal = (self.plane == other_.plane
                  and self.polygon2D == other_.polygon2D)
        for contour, other_contour in zip(self.contours, other_.contours):
            equal = (equal and contour == other_contour)
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point == other_point)
        return equal
    
# """
#     @classmethod
#     def from_step(cls, arguments, object_dict):
#         contours = []
#         contours.append(object_dict[int(arguments[1][0][1:])])
        
#         plane = Plane3D.from_points(contours[0].points)
#         contours[0].points, polygon2D = cls._repair_points_and_polygon2d(contours[0].points, plane)
#         points = [p.copy() for p in contours[0].points[:]]

#         return cls(contours, plane=plane, points=points, polygon2D=polygon2D, name=arguments[0][1:-1])
# """
    @classmethod
    def _repair_points_and_polygon2d(cls, points, plane):
        if points[0] == points[-1]:
            points = points[:-1]
        polygon_points = [p.To2D(plane.origin, plane.vectors[0], plane.vectors[1]) for p in points]
        repaired_points = [p.copy() for p in points]
        polygon2D = Polygon2D(polygon_points)
        if polygon2D.SelfIntersect()[0]:
            repaired_points = [repaired_points[1]]+[repaired_points[0]]+repaired_points[2:]
            polygon_points = [polygon_points[1]]+[polygon_points[0]]+polygon_points[2:]
            if polygon_points[0] == polygon_points[-1]:
                repaired_points = repaired_points[:-1]
                polygon_points = polygon_points[:-1]
            polygon2D = Polygon2D(polygon_points)
        return repaired_points, polygon2D

    @classmethod
    def from_contours3d(cls, contours3d, name=''):
        """
        :param contours3d: The face's contour3D
        :type contours3d: Contour3D
        """
        
        contour_points = [p.copy() for p in contours3d[0].tessel_points[:]]
        plane = Plane3D.from_points(contour_points)
        O, x, y = plane.origin, plane.vectors[0], plane.vectors[1]
        contours2d = []
        for contour in contours3d :
            prim = []
            for edge in contour.edges : 
                prim.append(edge.To2D(O,x,y))
            contours2d.append(Contour2D(prim))
        
        return cls(contours2d, plane, points=None, polygon2D=None, name=name)
    
    # @classmethod
    # def from_two_contours(cls, contour1, contour2) :
    #     #test, which is the smallest
    #     poly1, poly2 = contour1.polygon, contour2.polygon
    #     area1, area2 = poly1.Area(), poly2.Area()
    #     if area1 < area2 :
    #         basis_contour, other_contour = contour2, contour1
    #     else :
    #         basis_contour, other_contour = contour1, contour2
    #     #creation of polygon
    #     basis_points, other_points = basis_contour.tessel_points, other_contour.tessel_points
    #     basis_polygon, other_polygon = Polygon2D(basis_points), Polygon2D(other_points)
    #     #test to know if other_points in basis_polygon
    #     points_in_polygon = []
    #     for point in other_points :
    #         inside = basis_polygon.PointBelongs(point)
    #         if inside :
    #             points_in_polygon.append(point)
        
    #     fig, ax = plt.subplots()
    #     ax.set_aspect('equal')
    #     [pt.MPLPlot(ax=ax) for pt in basis_points]
    #     [pt.MPLPlot(ax=ax, color='b') for pt in other_points]
    #     [pt.MPLPlot(ax=ax, color='r') for pt in points_in_polygon]
        
        
    #     raise NotImplementedError
        
    #     return cls()
    
    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_contour = [subcontour.Rotation(center, axis, angle, copy=True) for subcontour in self.contour]
            new_plane = self.plane.Rotation(center, axis, angle, copy=True)
            new_points = [p.Rotation(center, axis, angle, copy=True) for p in self.points]
            return (new_contour, new_plane, new_points, self.polygon2D, self.name)
        else:
            for contour in self.contours:
                contour.Rotation(center, axis, angle, copy=False)
            for point in self.points:
                point.Rotation(center, axis, angle, copy=False)
            self.plane.Rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_contour = [subcontour.Translation(offset, copy=True) for subcontour in self.contours]
            new_plane = self.plane.Translation(offset, copy=True)
            new_points = [p.Translation(offset, copy=True) for p in self.points]
            return (new_contour, new_plane, new_points, self.polygon2D, self.name)
        else:
            for contour in self.contours:
                contour.Translation(offset, copy=False)
            for point in self.points:
                point.Translation(offset, copy=False)
            self.plane.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_plane = self.plane.frame_mapping(frame, side, copy=True)
            return PlaneFace3D(self.contours, new_plane, None, None, self.name)
        else:
            self.plane.frame_mapping(frame, side, copy=False)
            for contour in self.contours3d :
                contour.frame_mapping(frame, side, copy=False)
            self.setup_planeface(self.contours, self.plane, name=self.name)

    def copy(self):
        new_contours = [contour.copy() for contour in self.contours]
        new_plane = self.plane.copy()
        new_points = [p.copy() for p in self.points]
        return PlaneFace3D(new_contours, new_plane, new_points, self.polygon2D.copy(), self.name)

    def average_center_point(self):
        """
        excluding holes
        """
        points = self.points
        nb = len(points)
        x = npy.sum([p[0] for p in points]) / nb
        y = npy.sum([p[1] for p in points]) / nb
        z = npy.sum([p[2] for p in points]) / nb
        return Point3D((x,y,z))

    def triangulation(self):
        points_3D = []
        vertices = []
        segments = []
        holes = []
        total_len = 0
        for i, contour in enumerate(self.contours):
            
            points_2D = contour.tessel_points
            vertices.extend([tuple(p.vector) for p in points_2D])
            
            # #### Problem with BSplineCurve3D in 2D
            # print('contour.primitives',contour.primitives)
            # fig, ax = plt.subplots()
            # ax.set_aspect('equal')
            # [prim.MPLPlot(ax=ax) for prim in contour.primitives]
            # for prim in contour.primitives :
            #     prim.points[0].MPLPlot(ax=ax, color='r')
            #     prim.points[-1].MPLPlot(ax=ax, color='b')
            
            # fig, ax = plt.subplots()
            # ax.set_aspect('equal')
            # [pt.MPLPlot(ax=ax, color='g') for pt in points_2D]
            
            if len(vertices) != len(set(vertices)):
                raise ValueError('Clean_points problem')
            
            len_points = len(points_2D)
            segments += [[a+total_len, a+total_len+1] for a in range(len_points-1)]+[[len_points+total_len-1, 0+total_len]]
            total_len += len_points
            points_3D.extend([pt.To3D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1]) for pt in points_2D])
            if i > 0:
                polygon2D = Polygon2D(points_2D)
                mid_point_2D = contour.average_center_point()
                holes.append(mid_point_2D.vector)
                if not polygon2D.PointBelongs(mid_point_2D):
                    warnings.warn('average_center_point is not included inside its contour.')
        
        if holes:
            tri = {'vertices': vertices, 'segments': segments, 'holes': holes}
        else:
            tri = {'vertices': vertices, 'segments': segments}
        t = triangle.triangulate(tri, 'p')
        if 'triangles' in t:
            triangles = t['triangles'].tolist()
            return points_3D, [triangles]
        else:
            return None, [None]

    def distance_to_point(self, point, return_other_point=False):
        ## """
        ## Only works if the surface is planar
        ## TODO : this function does not take into account if Face has holes
        ## """
        # On projette le point sur la surface plane
        # Si le point est à l'intérieur de la face, on retourne la distance de projection
        # Si le point est à l'extérieur, on projette le point sur le plan
        # On calcule en 2D la distance entre la projection et le polygone contour
        # On utilise le theroeme de Pytagore pour calculer la distance minimale entre le point et le contour

        projected_pt = point.PlaneProjection3D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1])
        projection_distance = point.point_distance(projected_pt)

        if self.point_on_face(projected_pt):
            if return_other_point:
                return projection_distance, projected_pt
            return projection_distance

        point_2D = point.To2D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1])

        border_distance, other_point = self.polygon2D.PointBorderDistance(point_2D, return_other_point=True)

        other_point = other_point.To3D(self.plane.origin , self.plane.vectors[0], self.plane.vectors[1])

        if return_other_point:
            return (projection_distance**2 + border_distance**2)**0.5, other_point
        return (projection_distance**2 + border_distance**2)**0.5

    def distance_to_face(self, face2, return_points=False):
        ## """
        ## Only works if the surface is planar
        ## TODO : this function does not take into account if Face has holes
        ## TODO : TRAITER LE CAS OU LA DISTANCE LA PLUS COURTE N'EST PAS D'UN SOMMET
        ## """
        # On calcule la distance entre la face 1 et chaque point de la face 2
        # On calcule la distance entre la face 2 et chaque point de la face 1

        if self.face_intersection(face2) is not None:
            return 0, None, None

        polygon1_points_3D = [Point3D(p.vector) for p in self.contours3d[0].tessel_points]
        polygon2_points_3D = [Point3D(p.vector) for p in face2.contours3d[0].tessel_points]

        distances = []
        if not return_points:
            d_min = face2.distance_to_point(polygon1_points_3D[0])
            for point1 in polygon1_points_3D[1:]:
                d = face2.distance_to_point(point1)
                if d < d_min:
                    d_min = d
            for point2 in polygon2_points_3D:
                d = self.distance_to_point(point2)
                if d < d_min:
                    d_min = d
            return d_min

        else:
            for point1 in polygon1_points_3D:
                d, other_point = face2.distance_to_point(point1, return_other_point=True)
                distances.append((d, point1, other_point))
            for point2 in polygon2_points_3D:
                d, other_point = self.distance_to_point(point2, return_other_point=True)
                distances.append((d, point2, other_point))
        
        d_min, point_min, other_point_min = distances[0]
        for distance in distances[1:]:
            if distance[0] < d_min:
                d_min = distance[0]
                point_min = distance[1]
                other_point_min = distance[2]
        
        return d_min, point_min, other_point_min

    def point_on_face(self, point):
        """
        Tells you if a point is on the 3D face and inside its contour
        """
        ## Only works if the surface is planar
        ## TODO : this function does not take into account if Face has holes
        
        point_on_plane = self.plane.point_on_plane(point)
        # The point is not in the same plane
        if not point_on_plane:
            print('point not on plane so not on face')
            return False

        point_2D = point.To2D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1])

        ### PLOT ###
#        ax = self.polygon2D.MPLPlot()
#        point_2D.MPLPlot(ax)
        ############

        if not self.polygon2D.PointBelongs(point_2D):
            ### PLOT ###
#            ax = self.polygon2D.MPLPlot()
#            point_2D.MPLPlot(ax)
#            ax.set_title('DEHORS')
#            ax.set_aspect('equal')
            ############
            return False
        ### PLOT ###
#        ax = self.polygon2D.MPLPlot()
#        point_2D.MPLPlot(ax)
#        ax.set_title('DEDANS')
#        ax.set_aspect('equal')
        ############
        return True

    def edge_intersection(self, edge):
        linesegment = LineSegment3D(*edge.points)
        intersection_point = self.plane.linesegment_intersection(linesegment)

        if intersection_point is None:
            return None

        point_on_face_boo = self.point_on_face(intersection_point)
        if not point_on_face_boo:
            return None

        return intersection_point

    def linesegment_intersection(self, linesegment, abscissea=False):
        if abscissea:
            intersection_point, intersection_abscissea = self.plane.linesegment_intersection(linesegment, True)
        else:
            intersection_point = self.plane.linesegment_intersection(linesegment)

        if intersection_point is None:
            if abscissea:
                return None, None
            return None
        
        point_on_face_boo = self.point_on_face(intersection_point)
        # ## PLOT ###
        # if point_on_face_boo:
        #     ax = self.plot()
        #     linesegment.MPLPlot(ax)
        #     Point3D(intersection_point.vector).MPLPlot(ax)
        #     self.plane.MPLPlot(ax)
        #     ax.set_aspect('equal')

        #     print('=>', self.plane.normal.Dot(intersection_point-self.plane.origin))
        #     print('point_on_face_boo', point_on_face_boo)
        # ###########
        if not point_on_face_boo:
            if abscissea:
                return None, None
            return None

        if abscissea:
            return intersection_point, intersection_abscissea
        return intersection_point

    def face_intersection(self, face2):
        ## """
        ## Only works if the surface is planar
        ## TODO : this function does not take into account if Face has holes
        ## """
        bbox1 = self.bounding_box
        bbox2 = face2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
            return None

        intersection_points = []

        for edge2 in face2.contours3d[0].edges:
            intersection_point = self.edge_intersection(edge2)
            if intersection_point is not None:
                intersection_points.append(intersection_point)

        for edge1 in self.contours3d[0].edges:
            intersection_point = face2.edge_intersection(edge1)
            if intersection_point is not None:
                intersection_points.append(intersection_point)

        if not intersection_points:
            return None

        return intersection_points

    def plot(self, ax=None):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

#        x = [p[0] for edge in self.contours[0].edges for p in edge.points]
#        y = [p[1] for edge in self.contours[0].edges for p in edge.points]
#        z = [p[2] for edge in self.contours[0].edges for p in edge.points]
#        print(x,y,z)
        x = [p[0] for p in self.contours[0].tessel_points]
        y = [p[1] for p in self.contours[0].tessel_points]
        z = [p[2] for p in self.contours[0].tessel_points]

        ax.scatter(x, y, z)
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        for edge in self.contours[0].edges:
            for point1, point2 in (edge.points, edge.points[1:]+[edge.points[0]]):
                xs = [point1[0], point2[0]]
                ys = [point1[1], point2[1]]
                zs = [point1[2], point2[2]]
                line = mpl_toolkits.mplot3d.art3d.Line3D(xs, ys, zs)
                ax.add_line(line)

        plt.show()
        return ax
    
    def minimum_distance(self, other_face, return_points=False) :
        if other_face.__class__ is CylindricalFace3D :
            p1, p2 = other_face.minimum_distance_points_cyl(self)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
        
        if other_face.__class__ is PlaneFace3D : 
            dmin, p1, p2 = self.distance_to_face(other_face, return_points=True)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
            
        if other_face.__class__ is ToroidalFace3D : 
            p1, p2 = other_face.minimum_distance_points_plane(self)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
            
        else :
            return NotImplementedError 


class CylindricalFace3D(Face3D): 
    """
    :param contours2d: The cylinder's contour2D 
    :type contours2d: Contour2D
    :param cylindricalsurface3d: Information about the Cylinder
    :type cylindricalsurface3d: CylindricalSurface3D
    :param points: contours2d's point
    :type points: List of Point2D
    
    :Example: 
        >>> contours2d is rectangular and will create a classic cylinder with x= 2*pi*radius, y=h
    """      
    def __init__(self, contours2d, cylindricalsurface3d, points=None, name=''):
        self.radius = cylindricalsurface3d.radius
        self.center = cylindricalsurface3d.frame.origin
        self.normal = cylindricalsurface3d.frame.w
        edge1 = Circle3D(cylindricalsurface3d.frame.origin, self.radius, cylindricalsurface3d.frame.w )
        edge2 = Circle3D(cylindricalsurface3d.frame.origin+contours2d[0].primitives[0].points[1][1]*cylindricalsurface3d.frame.w, self.radius, cylindricalsurface3d.frame.w)
        ctr = [Contour3D([edge1, edge2], name='')]
        Face3D.__init__(self, ctr)
        
        self.contours2d = contours2d 
        self.cylindricalsurface3d = cylindricalsurface3d 
        if points is None:
            self.points = self.contours2d[0].tessel_points 
        else:
            self.points = points 
        self.name = name 
        
        # # CHECK
        # for pt in self.points:
        #     if not self.frame.point_on_plane(pt):
        #         print('WARNING', pt, 'not on', self.frame.__dict__)
        #         print('dot =', self.frame.normal.Dot(pt-self.frame.origin))
        #         raise ValueError
    
    @classmethod
    def from_contour3d(cls, contours3d, cylindricalsurface3d, name=''):
        """
        :param contours3d: The cylinder's contour3D
        :type contours3d: Contour3D
        :param cylindricalsurface3d: Information about the Cylinder
        :type cylindricalsurface3d: CylindricalSurface3D
        
        :Example:
            >>> contours3d is [Arc3D, LineSegment3D, Arc3D], the cylinder's bones
        """
        
        frame = cylindricalsurface3d.frame
        radius = cylindricalsurface3d.radius
        size = len(contours3d[0].edges)
        
        if contours3d[0].edges[0].__class__ is LineSegment3D and contours3d[0].edges[1].__class__ is Arc3D :
            return CylindricalFace3D.from_arc3d(contours3d[0].edges[0], contours3d[0].edges[1], cylindricalsurface3d)
        
        if contours3d[0].edges[0].__class__ is Arc3D and contours3d[0].edges[2].__class__ is Arc3D and size <= 4:
            
            arc1, arc2 = contours3d[0].edges[0], contours3d[0].edges[2]
            c1, c2 = frame.NewCoordinates(arc1.center), frame.NewCoordinates(arc2.center)
            hmin, hmax = min(c1.vector[2], c2.vector[2]), max(c1.vector[2], c2.vector[2])
            n1 = arc1.normal
            if n1 == -frame.w :
                arc1.setup_arc(arc1.start, arc1.interior, arc1.end, -arc1.normal)
            start1, end1 = arc1.start, arc1.end
            theta1_1, theta1_2 = posangle_arc(start1, end1, radius, frame)
            if not(math.isclose(arc1.angle,abs(theta1_1-theta1_2), abs_tol=1e-4)) :
                if math.isclose(theta1_1, 0, abs_tol=1e-4) :
                    theta1_1 = 2*math.pi
                elif math.isclose(theta1_2, 0, abs_tol=1e-4) :
                    theta1_2 = 2*math.pi
                else :
                    # if theta1_2 > theta1_1 :
                    #     theta1_2 -= math.pi
                    # else :
                    #     theta1_1 -= math.pi
                    print('arc1.angle', arc1.angle)
                    print('theta1_1, theta1_2', theta1_1, theta1_2)
                    raise NotImplementedError
                
            offset1, angle1 = offset_angle(arc1.is_trigo, theta1_1, theta1_2)
            pt1, pt2, pt3, pt4 = Point2D((offset1, hmin)), Point2D((offset1, hmax)), Point2D((offset1+angle1, hmax)), Point2D((offset1+angle1, hmin))
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]    
            points = contours2d[0].tessel_points
            
        else : 
            contours2d = CylindricalFace3D.contours3d_to2d(contours3d, cylindricalsurface3d)
            points = contours2d[0].tessel_points
            
        return cls(contours2d, cylindricalsurface3d, points, name=name)
    
    @classmethod 
    def from_arc3d(cls, lineseg, arc, cylindricalsurface3d): #Work with 2D too
        """
        :param lineseg: The segment which represent the extrusion of the arc
        :type lineseg: LineSegment3D/2D
        :param arc: The Arc circle to extrude
        :type arc: Arc3D/2D, Circle3D/2D
        :param cylindricalsurface3d: Information about the Cylinder
        :type cylindricalsurface3d: CylindricalSurface3D
        
        Particularity : the frame is the base of the cylinder, it begins there and go in the normal direction
        """
        radius = cylindricalsurface3d.radius
        frame = cylindricalsurface3d.frame
        normal, center = frame.w, frame.origin
        offset = 0
        if arc.__class__.__name__ == 'Circle3D' or arc.__class__.__name__ == 'Circle2D' :
            
            frame_adapt = cylindricalsurface3d.frame
            theta = arc.angle
                
        else :
            point12d = arc.start
            if point12d.__class__ is Point3D :
                point12d = point12d.To2D(center, frame.u, frame.v) #Using it to put arc.start at the same height 
            point13d = point12d.To3D(center, frame.u, frame.v)
            if arc.start.__class__ is Point2D :
                u_g2d = Vector2D((arc.start - arc.center).vector)
                u = u_g2d.To3D(center, frame.u, frame.v)
                u.Normalize()
            else :
                u = Vector3D((point13d - center).vector)
                u.Normalize()
            v = normal.Cross(u)
            v.Normalize()
            
            point_last = arc.end
            if point_last.__class__ is Point3D :
                point_last = point_last.To2D(center, u, v)
                
            x, y = point_last.vector[0], point_last.vector[1]
        
            theta = math.atan2(y, x)
            if theta < 0 or math.isclose(theta, 0, abs_tol=1e-9):
                if arc.angle>math.pi :
                    theta += 2*math.pi
                else :
                    offset = theta
                    theta = -theta
                    
            frame_adapt = Frame3D(center, u, v, normal)

        cylindersurface3d = CylindricalSurface3D(frame_adapt, radius)
        segbh = LineSegment2D(Point2D((offset,0)), Point2D((offset,lineseg.Length())))
        circlestart = LineSegment2D(segbh.points[1], segbh.points[1]+Point2D((theta,0)))
        seghb = LineSegment2D(circlestart.points[1],circlestart.points[1]-segbh.points[1]+segbh.points[0])
        circlend = LineSegment2D(seghb.points[1],segbh.points[0])
        
        edges = [segbh, circlestart, seghb, circlend]
        return cls([Contour2D(edges)], cylindersurface3d, points=None, name='')
    
    def points2d_to3d(self, all_contours_points, radius, frame3d) :
        Points3D = []
        for listpt in  all_contours_points: 
            for enum, pt in enumerate(listpt) :
                Points3D.append(Point3D(Vector3D([radius*math.cos(pt[0]),radius*math.sin(pt[0]),pt[1]])))
        Points_3D = [frame3d.OldCoordinates(point) for point in Points3D]
        return Points_3D
    
    def points3d_to2d(points3d, radius) :
        points_2D = []
        for pt in points3d :
            x, y, h = pt[0], pt[1], pt[2]
            
            u1, u2 = x/radius, y/radius
            theta = sin_cos_angle(u1, u2)
            
            points_2D.append(Point2D([theta, h])) 
        i0, i2pi, iangle, ih = 0, 0, 0, 0
        for enum,point in enumerate(points_2D) :
            if enum == 0 :
                h = point[1]
                angle = point[0]
            if math.isclose(point[0], 0, abs_tol=5e-2) :
                i0 += 1
            elif math.isclose(point[0], 2*math.pi, abs_tol=5e-2) :
                i2pi += 1
            elif math.isclose(point[0], angle, abs_tol=5e-2) :
                iangle += 1
            elif math.isclose(point[1], h, abs_tol=1e-6) :
                ih += 1
        
        if i2pi/len(points_2D)>0.5 or i0/len(points_2D)>0.5:
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((2*math.pi, point.vector[1])))
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif iangle/len(points_2D)>0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((angle, point.vector[1]))) 
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif ih/len(points_2D)>0.5 :
            new_points2d = [pt.copy() for pt in points_2D]
            new_points2d.sort(key=lambda pt: pt[0])
            if i2pi == 2 :
                pt1, pt2 = new_points2d[0], new_points2d[-1]
                pt1.vector[0] = 0-0.0000001
                pt2.vector[0] = 2*math.pi+0.0000001
                points_2D = [pt1, pt2]
            else :
                points_2D = [new_points2d[0], new_points2d[-1]]
        
        plus_7pi4, moins_7pi4 = [], []
        for enum, pt in enumerate(points_2D) :
            if pt.vector[0]>7*math.pi/4 :
                plus_7pi4.append(enum)
            elif pt.vector[0]<math.pi/4 :
                moins_7pi4.append(enum)
            
        if len(moins_7pi4) > 3*len(plus_7pi4) :
            for pos in plus_7pi4 :
                new_pt = points_2D[pos].copy() - Point2D((2*math.pi, 0))
                points_2D[pos] = new_pt
        
        return points_2D

    def contours3d_to2d(contours3d, cylindricalsurface3d) :
        frame = cylindricalsurface3d.frame
        n = frame.w
        radius = cylindricalsurface3d.radius
        
        primitives, start_end, all_points = [], [], []
        for edge in contours3d[0].edges :
            new_points = [frame.NewCoordinates(pt) for pt in edge.points]
            if edge.__class__ is Arc3D :
                if edge.normal == n or edge.normal == -n:
                    start2d, end2d = CylindricalFace3D.points3d_to2d(new_points, radius)
                    angle2d = abs(end2d[0]-start2d[0])
                    if math.isclose(edge.angle, 2*math.pi, abs_tol=1e-6):
                        if start2d == end2d :
                            if math.isclose(start2d.vector[0], 2*math.pi, abs_tol = 1e-6) :
                                end2d = end2d - Point2D((2*math.pi, 0))
                            else :
                                end2d = end2d + Point2D((2*math.pi, 0))
                    elif not(math.isclose(edge.angle, angle2d, abs_tol=1e-2)) :
                        # if math.isclose(angle2d, 2*math.pi, abs_tol=1e-2) :
                        if start2d[0] < end2d[0] :
                            end2d = start2d + Point2D((edge.angle,0))
                        else :
                            end2d = start2d - Point2D((edge.angle,0))
                    ls_toadd = LineSegment2D(start2d, end2d)
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        primitives.append([ls_toadd])
                        all_points.extend(ls_toadd.points)
                        start_end.append(ls_toadd.points)
                    
                else :
                    points2d = CylindricalFace3D.points3d_to2d(new_points, radius)
                    lines = []
                    for k in range(0, len(points2d)-1) :
                        lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                    points, prim_list = [], []
                    for ls_toadd in lines :        
                        same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                        if same is False :
                            prim_list.append(ls_toadd)
                            points.extend(ls_toadd.points)
                    if len(points) > 0 : 
                        all_points.extend(points)
                        primitives.append(prim_list)
                        start_end.append([points[0], points[-1]])
                    
            elif edge.__class__ is LineSegment3D :
                start2d, end2d = CylindricalFace3D.points3d_to2d(new_points, radius)
                ls_toadd = LineSegment2D(start2d, end2d)
                same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                if same is False :
                    primitives.append([ls_toadd])
                    all_points.extend(ls_toadd.points)
                    start_end.append(ls_toadd.points)
            
            # elif edge.__class__ is BSplineCurve3D :
                
            #     points2d = CylindricalFace3D.points3d_to2d(new_points, radius)
            #     print('edge.knot_multiplicities', edge.knot_multiplicities)
            #     prim_list = [BSplineCurve2D(edge.degree, points2d, edge.knot_multiplicities, edge.knots)]
            #     all_points.extend(points2d)
            #     primitives.append(prim_list)
            #     start_end.append([points2d[0], points2d[-1]])
            
            else :
                points2d = CylindricalFace3D.points3d_to2d(new_points, radius)
                lines = []
                for k in range(0, len(points2d)-1) :
                    lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                points, prim_list = [], []
                for ls_toadd in lines :        
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        prim_list.append(ls_toadd)
                        points.extend(ls_toadd.points)
                if len(points) > 0 : 
                    all_points.extend(points)
                    primitives.append(prim_list)
                    start_end.append([points[0], points[-1]])
                    
        # fig, ax = plt.subplots() 
        # [pt.MPLPlot(ax=ax) for pt in all_points]
        
        # fig, ax = plt.subplots() 
        # for list_prim in primitives :
        #     [prim.MPLPlot(ax=ax) for prim in list_prim]  
        # maxi, mini = max(pt.vector[0] for pt in all_points), min(pt.vector[0] for pt in all_points)
        # if (mini > 0 and maxi < math.pi) or (mini > math.pi and maxi < 2*math.pi) :
        #     fig, ax = plt.subplots() 
        #     for enum, list_prim in enumerate(primitives) :
        #         [prim.MPLPlot(ax=ax) for prim in list_prim]
        #         start_end[enum][0].MPLPlot(ax=ax, color='g')
        #         start_end[enum][1].MPLPlot(ax=ax, color='r')
               
        #     for n, s_e in enumerate(start_end) :
        #         if n == 0 :
        #             start1, end1 = s_e
        #             start2, end2 = start_end[n+1]
        #         elif n == len(start_end)-1 :
        #             continue
        #         else :
        #             start2, end2 = start_end[n+1]
                
        #         d1, d2 = start1.point_distance(start2), start1.point_distance(end2)
        #         d3, d4 = end1.point_distance(start2), end1.point_distance(end2)
        #         posmin, dmin = min_pos([d1, d2, d3, d4])
        #         if posmin == 2 or posmin == 1:
        #             continue
        #         else:
        #             new_prim = []
        #             for prim in primitives[n+1] :
        #                 new_prim.append(prim.reverse())
        #             primitives[n+1] = new_prim
        #             new_start = end2
        #             end2 = start2
        #             start2 = new_start
        #         start1, end1 = start2, end2
              
        #     points_primitives = []
        #     for list_prim in primitives :
        #         for prim in list_prim :
        #             points_primitives.extend(prim.points)
        #     # points = delete_double_point(points_primitives)
        #     primitives = Face3D.create_primitives(points_primitives)
        #     fig, ax = plt.subplots() 
        #     [p.MPLPlot(ax=ax) for p in primitives]
        #     contour2d = [Contour2D(primitives)]
            
        #     # prim2 = Face3D.create_primitives(contour2d[0].tessel_points)
        #     # fig, ax = plt.subplots() 
        #     # [p.MPLPlot(ax=ax) for p in prim2]
            
        #     # contour2d = [Contour2D(prim2)]
        #     # raise NotImplementedError
            
        # else :
        points_se, primitives_se = [], []
        for double in start_end : 
            primitives_se.append(LineSegment2D(double[0], double[1]))
            points_se.extend(double)
        poly_se = Polygon2D(points_se)
        
        xmax, xmin = max(pt[0] for pt in points_se), min(pt[0] for pt in points_se)
        ymax, ymin = max(pt[1] for pt in points_se), min(pt[1] for pt in points_se) 
        pt1, pt2, pt3, pt4 = Point2D((xmin, ymin)), Point2D((xmin, ymax)), Point2D((xmax, ymin)), Point2D((xmax, ymax))
        diag1, diag2 = LineSegment2D(pt1, pt4), LineSegment2D(pt2, pt3)
        diag1_cut, diag2_cut = [], []
        diag1_pointcut, diag2_pointcut = [], []
        for enum,l in enumerate(primitives_se) :
            cut1 = diag1.line_intersection(l)
            cut2 = diag2.line_intersection(l)
            if cut1 is not None : 
                diag1_cut.append(enum)
                diag1_pointcut.append(cut1)
            if cut2 is not None : 
                diag2_cut.append(enum)
                diag2_pointcut.append(cut2)
        
        points_common = []
        for enum1,pos1 in enumerate(diag1_cut) :
            for enum2,pos2 in enumerate(diag2_cut) :
                if pos1 == pos2 :
                    points_common.append(primitives_se[pos1].points)
        
        # fig, ax = plt.subplots() 
        # [l.MPLPlot(ax=ax) for l in primitives_se]
        # [pt.MPLPlot(ax=ax, color='g') for pt in points_se]
        # pt1.MPLPlot(ax=ax, color='r')
        # pt2.MPLPlot(ax=ax, color='r')            
        # pt3.MPLPlot(ax=ax, color='r')
        # pt4.MPLPlot(ax=ax, color='r')   
        # diag1.MPLPlot(ax=ax)
        # diag2.MPLPlot(ax=ax)
        # for couple in points_common :
        #     [pt.MPLPlot(ax=ax, color='b') for pt in couple] 
        
        if len(points_common) >= 1 :
            solve = False
            for couple in points_common :
                check1, check2 = poly_se.PointBelongs(couple[0]), poly_se.PointBelongs(couple[1])
                start, end = couple[0].vector[0], couple[1].vector[0]
                if math.isclose(start, end, abs_tol = 5e-2) :
                    intersect = min(start, end)
                    if math.isclose(intersect, math.pi, abs_tol = 5e-2) :
                        # all_points = check_singularity(all_points)
                        # intersect = 0
                        
                        ##################### NEW
                        points_sing = check_singularity(all_points)
                        pt0, pt2pi = 0, 0
                        for pt in points_sing :
                            if math.isclose(pt.vector[0], 0, abs_tol=1e-2):
                                pt0 += 1
                            elif math.isclose(pt.vector[0], 2*math.pi, abs_tol=1e-2):
                                pt2pi += 1
                        points_sing.sort(key=lambda pt: pt[1])
                        points_sing.sort(key=lambda pt: pt[0])
                        if pt2pi !=0 and pt0 == 0 :
                            points = [pt.copy() for pt in points_sing[::-1]]
                            points_sing = points
                        points_range = CylindricalFace3D.range_closest(points_sing, radius, frame)
                        all_points = delete_double_point(points_range)
                        break
                        #######################
                        
                    # if math.isclose(intersect, 0, abs_tol = 1e-6) or math.isclose(intersect, 2*math.pi, abs_tol = 1e-6) or (not check1 or not check2):
                    elif math.isclose(intersect, 0, abs_tol = 1e-6) or math.isclose(intersect, 2*math.pi, abs_tol = 1e-6) or (not check1 or not check2):
                        all_points = check_singularity(all_points)
                        
                        points_cleaned = delete_double_point(all_points)
                        all_points = [pt.copy() for pt in points_cleaned]
                        all_points.sort(key=lambda pt: pt[0])
                        d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                        if d2 < d1 :
                            last = all_points[-1].copy()
                            all_points[-1] = all_points[-2].copy()
                            all_points[-2] = last
                        break
                    else :
                        points = []
                        for list_prim in primitives :
                            for k,prim in enumerate(list_prim) :
                                new_list_points = []
                                change = 0
                                for pt in prim.points :
                                    if pt[0] < intersect :
                                        change += 1 
                                        if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                            new_list_points.append(Point2D((intersect + 2*math.pi, pt[1])))
                                        else :
                                            new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    elif math.isclose(pt[0], intersect, abs_tol=1e-1) :
                                        change += 1
                                        new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    else :
                                        new_list_points.append(pt)
                                if change > 0 :
                                    points.extend(new_list_points)
                                    # list_prim[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                                else :
                                    points.extend(prim.points)
                                    continue
                        points_cleaned = delete_double_point(points)
                        all_points = Face3D.range_trigo(points_cleaned)
                    solve = True
                else :
                    points_cleaned = delete_double_point(all_points)
                    all_points = [pt.copy() for pt in points_cleaned]
                    all_points.sort(key=lambda pt: pt[0])
                    d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                    if d2 < d1 :
                        last = all_points[-1].copy()
                        all_points[-1] = all_points[-2].copy()
                        all_points[-2] = last
        else :
            points_cleaned = delete_double_point(all_points)
            all_points = [pt.copy() for pt in points_cleaned]
            all_points.sort(key=lambda pt: pt[0])
            d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
            if d2 < d1 :
                last = all_points[-1].copy()
                all_points[-1] = all_points[-2].copy()
                all_points[-2] = last
            
        primitives = Face3D.create_primitives(all_points)
        
        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax, color='g') for pt in all_points]
        # all_points[0].MPLPlot(ax=ax, color='m')
        # all_points[1].MPLPlot(ax=ax, color='r')
        # all_points[-1].MPLPlot(ax=ax)
        # all_points[-2].MPLPlot(ax=ax, color='b')
        # [p.MPLPlot(ax=ax) for p in primitives]
        
        l_vert = LineSegment2D((pt2+pt4)/2, (pt1+pt3)/2)
        solve = False
        for prim in primitives :
            if solve :
                break
            intersect = prim.line_intersection(l_vert)
            if intersect is not None : 
                x_intersect = intersect.vector[0]
                y_intersect = intersect.vector[1]
                value1, value2 = ymax-0.2*(ymax-ymin), ymin+0.2*(ymax-ymin)
                if y_intersect < max(value1, value2) and y_intersect > min(value1, value2) :
                    points = []
                    for k, prim in enumerate(primitives) :
                        new_list_points, change = [], 0
                        for pt in prim.points :
                            if pt[0] < x_intersect :
                                change += 1
                                if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                    new_list_points.append(Point2D((x_intersect + 2*math.pi, pt[1])))
                                else :
                                    new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                            else :
                                new_list_points.append(pt)
                        if change > 0:
                            points.extend(new_list_points)
                            primitives[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                        else :
                            points.extend(prim.points)
                            continue
                    solve = True
                    points_cleaned = delete_double_point(points)
                    all_points = Face3D.range_trigo(points_cleaned)
                    primitives = Face3D.create_primitives(all_points)
                        
            # fig, ax = plt.subplots()
            # [pt.MPLPlot(ax=ax) for pt in all_points]
            # [p.MPLPlot(ax=ax) for p in primitives]
            # # # intersect.MPLPlot(ax=ax, color='r')
            # # l_vert.MPLPlot(ax=ax)
            contour2d = [Contour2D(primitives)]
        return contour2d
    
    def range_closest(list_point, radius, frame) :
        points_set = delete_double_point(list_point)
        points_set3D = CylindricalFace3D.points2d_to3d(None, [points_set], radius, frame)
        
        points_3dint = [points_set3D[0]]
        points_2dint = [points_set[0]]
        s = 1
        for k in range(1, len(points_set)) :
            closest = points_set3D[s]
            while closest is None : 
                s += 1
                closest = points_set3D[s]
            dist_min = (points_3dint[-1] - closest).Norm()
            pos = s
            for i in range (s+1, len(points_set3D)) :
                close_test = points_set3D[i]
                if close_test is None : 
                    continue
                else :
                    dist_test = (points_3dint[-1] - close_test).Norm()
                    if dist_test <= dist_min:
                        dist_min = dist_test
                        closest = close_test
                        pos = i
            points_2dint.append(points_set[pos])
            points_set3D[pos] = None
        
        return points_2dint        
    
    def triangulation(self, resolution=31):
        radius = self.radius
        frame3d = self.cylindricalsurface3d.frame
        all_contours_points = self.cut_contours(self.contours2d, resolution)
        
        ##########
        # fig, ax = plt.subplots()
        # [prim.MPLPlot(ax=ax) for prim in self.contours2d[0].primitives]
        # for listpt in all_contours_points :
        #     for pt in listpt :
        # # #         if pt is None : 
        # # #             continue
        # # #         else :
        #         pt.MPLPlot(ax=ax) 
        # [pt.MPLPlot(ax=ax, color='r') for pt in self.contours2d[0].tessel_points]
        #############
        
        Triangles, ts = [], []
        for k, listpt in enumerate(all_contours_points) :
            vertices=[]
            segments=[]
            for i, pt in enumerate(listpt):
                vertices.append(pt.vector)
                segments.append([i,i+1])
            
            segments[-1]=(len(listpt)-1,0) 
            tri = {'vertices': vertices, 'segments': segments}
            t = triangle.triangulate(tri, 'p')
            ts.append(t)
        
        seuil, k = 0, 0
        for t in ts :
            if 'triangles' in t:
                triangles = t['triangles'].tolist()
                for n,tri in enumerate(triangles):
                    for i in range (0,3):
                        tri[i]=tri[i]+seuil
                seuil += len(all_contours_points[k])
                k += 1
                Triangles.append(triangles)        
            else:
                Triangles.append(None)
        
        Points_3D = self.points2d_to3d(all_contours_points, radius, frame3d)
        pt3d, tangle = delete_double_pos(Points_3D, Triangles)
        
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # [pt.MPLPlot(ax=ax) for pt in pt3d]
        
        return pt3d, tangle 

    def MPLPlotpoints(self, ax=None, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        x, y, z = [], [], []
        for px, py, pz in self.points :
            x.append(px)
            y.append(py)
            z.append(pz)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color)
        return ax
    
    def MPLPlotcontours(self, ax=None, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        x = []
        y = []
        z = []
        for px, py, pz in self.contours[0].tessel_points :
            x.append(px)
            y.append(py)
            z.append(pz)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color)
        return ax
    
    def frame_mapping(self, frame, side, copy=True) :
        if copy:
            new_cylindricalsurface3d = CylindricalSurface3D.frame_mapping(frame, side, copy)
            return CylindricalFace3D(self.contours2d, new_cylindricalsurface3d, points=self.points, name=self.name)
        else:
            self.cylindricalsurface3d.frame_mapping(frame, side, copy=False)
            
    def minimum_maximum(self, contour2d, radius) :
        points = contour2d.tessel_points
        
        min_h, min_theta = min([pt[1] for pt in points]), min([pt[0] for pt in points])
        max_h, max_theta = max([pt[1] for pt in points]), max([pt[0] for pt in points]) 
        return min_h, min_theta, max_h, max_theta
            
    def minimum_distance_points_cyl(self, other_cyl) :
        r1, r2 = self.radius, other_cyl.radius
        min_h1, min_theta1, max_h1, max_theta1 = self.minimum_maximum(self.contours2d[0], r1)
        
        n1 = self.normal
        u1 = self.cylindricalsurface3d.frame.u
        v1 = self.cylindricalsurface3d.frame.v
        frame1 = Frame3D(self.center, u1, v1, n1)
        # st1 = Point3D((r1*math.cos(min_theta1), r1*math.sin(min_theta1), min_h1))
        # start1 = frame1.OldCoordinates(st1)
        
        min_h2, min_theta2, max_h2, max_theta2 = self.minimum_maximum(other_cyl.contours2d[0], r2)
             
        n2 = other_cyl.normal
        u2 = other_cyl.cylindricalsurface3d.frame.u
        v2 = other_cyl.cylindricalsurface3d.frame.v
        frame2 = Frame3D(other_cyl.center, u2, v2, n2)
        # st2 = Point3D((r2*math.cos(min_theta2), r2*math.sin(min_theta2), min_h2))
        # start2 = frame2.OldCoordinates(st2)
        
        w = other_cyl.center - self.center
        

        n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.Dot(n1), n1.Dot(u1), n1.Dot(v1), n1.Dot(n2), n1.Dot(u2), n1.Dot(v2)
        u1u1, u1v1, u1n2, u1u2, u1v2 = u1.Dot(u1), u1.Dot(v1), u1.Dot(n2), u1.Dot(u2), u1.Dot(v2)
        v1v1, v1n2, v1u2, v1v2 = v1.Dot(v1), v1.Dot(n2), v1.Dot(u2), v1.Dot(v2)
        n2n2, n2u2, n2v2 = n2.Dot(n2), n2.Dot(u2), n2.Dot(v2)
        u2u2, u2v2, v2v2 = u2.Dot(u2), u2.Dot(v2), v2.Dot(v2)
        
        w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.Dot(w), w.Dot(n1), w.Dot(u1), w.Dot(v1), w.Dot(n2), w.Dot(u2), w.Dot(v2)
        
        # x = (theta1, h1, theta2, h2)
        def distance_squared(x):
            return (n1n1*(x[1]**2) + u1u1*((math.cos(x[0]))**2)*(r1**2) + v1v1*((math.sin(x[0]))**2)*(r1**2)
                    + w2 + n2n2*(x[3]**2) + u2u2*((math.cos(x[2]))**2)*(r2**2) + v2v2*((math.sin(x[2]))**2)*(r2**2)
                    + 2*x[1]*r1*math.cos(x[0])*n1u1 + 2*x[1]*r1*math.sin(x[0])*n1v1 -2*x[1]*wn1 
                    - 2*x[1]*x[3]*n1n2 - 2*x[1]*r2*math.cos(x[2])*n1u2 - 2*x[1]*r2*math.sin(x[2])*n1v2
                    + 2*math.cos(x[0])*math.sin(x[0])*u1v1*(r1**2) - 2*r1*math.cos(x[0])*wu1 
                    - 2*r1*x[3]*math.cos(x[0])*u1n2 - 2*r1*r2*math.cos(x[0])*math.cos(x[2])*u1u2
                    - 2*r1*r2*math.cos(x[0])*math.sin(x[2])*u1v2 -2*r1*math.sin(x[0])*wv1
                    - 2*r1*x[3]*math.sin(x[0])*v1n2 - 2*r1*r2*math.sin(x[0])*math.cos(x[2])*v1u2
                    - 2*r1*r2*math.sin(x[0])*math.sin(x[2])*v1v2 + 2*x[3]*wn2 + 2*r2*math.cos(x[2])*wu2
                    + 2*r2*math.sin(x[2])*wv2 + 2*x[3]*r2*math.cos(x[2])*n2u2 + 2*x[3]*r2*math.sin(x[2])*n2v2
                    + 2*math.cos(x[2])*math.sin(x[2])*u2v2*(r2**2))
        
        x01 = npy.array([(min_theta1+max_theta1)/2, (min_h1+max_h1)/2,
                         (min_theta2+max_theta2)/2, (min_h2+max_h2)/2])
        x02 = npy.array([min_theta1, (min_h1+max_h1)/2,
                          min_theta2, (min_h2+max_h2)/2])
        x03 = npy.array([max_theta1, (min_h1+max_h1)/2,
                          max_theta2, (min_h2+max_h2)/2])
        
        minimax = [(min_theta1, min_h1, min_theta2, min_h2), (max_theta1, max_h1, max_theta2, max_h2)]
        
        res1 = scp.optimize.least_squares(distance_squared, x01, bounds=minimax)
        res2 = scp.optimize.least_squares(distance_squared, x02, bounds=minimax)
        res3 = scp.optimize.least_squares(distance_squared, x03, bounds=minimax)
        
        pt1 = Point3D((r1*math.cos(res1.x[0]),r1*math.sin(res1.x[0]),res1.x[1]))
        p1 = frame1.OldCoordinates(pt1)
        pt2 = Point3D((r2*math.cos(res1.x[2]),r2*math.sin(res1.x[2]),res1.x[3]))
        p2 = frame2.OldCoordinates(pt2)
        d = p1.point_distance(p2)
        result = res1
        
        res = [res2, res3]
        for couple in res :
            pttest1 = Point3D((r1*math.cos(couple.x[0]),r1*math.sin(couple.x[0]),couple.x[1]))
            pttest2 = Point3D((r2*math.cos(couple.x[2]),r2*math.sin(couple.x[2]),couple.x[3]))
            ptest1 = frame1.OldCoordinates(pttest1)
            ptest2 = frame2.OldCoordinates(pttest2)
            dtest = ptest1.point_distance(ptest2)
            if dtest < d :
                result = couple
                p1, p2 = ptest1, ptest2
        
        pt1_2d, pt2_2d = Point2D((result.x[0], result.x[1])), Point2D((result.x[2], result.x[3]))
        
        if not(self.contours2d[0].point_belongs(pt1_2d)) :
            #Find the closest one
            points_contours1 = self.contours2d[0].tessel_points
            
            poly1 = Polygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d, return_other_point=True)
            pt1 = Point3D((r1*math.cos(new_pt1_2d.vector[0]),
                           r1*math.sin(new_pt1_2d.vector[0]),
                           new_pt1_2d.vector[1]))
            p1 = frame1.OldCoordinates(pt1)        
        
        if not(other_cyl.contours2d[0].point_belongs(pt2_2d)) :
            #Find the closest one
            points_contours2 = other_cyl.contours2d[0].tessel_points
            
            poly2 = Polygon2D(points_contours2)
            d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d, return_other_point=True)
            pt2 = Point3D((r2*math.cos(new_pt2_2d.vector[0]),
                           r2*math.sin(new_pt2_2d.vector[0]),
                           new_pt2_2d.vector[1]))
            p2 = frame2.OldCoordinates(pt2)        
            
        return p1, p2
    
    def minimum_distance_points_plane(self, planeface): #Planeface with contour2D
        #### ADD THE FACT THAT PLANEFACE.CONTOURS : [0] = contours totale, le reste = trous
        r = self.radius
        min_h1, min_theta1, max_h1, max_theta1 = self.minimum_maximum(self.contours2d[0], r)
        
        n1 = self.normal
        u1 = self.cylindricalsurface3d.frame.u
        v1 = self.cylindricalsurface3d.frame.v
        frame1 = Frame3D(self.center, u1, v1, n1)
        # st1 = Point3D((r*math.cos(min_theta1), r*math.sin(min_theta1), min_h1))
        # start1 = frame1.OldCoordinates(st1)
        
        poly2d = planeface.polygon2D
        pfpoints = poly2d.points
        xmin, ymin = min([pt[0] for pt in pfpoints]), min([pt[1] for pt in pfpoints])
        xmax, ymax = max([pt[0] for pt in pfpoints]), max([pt[1] for pt in pfpoints])
        origin, vx, vy = planeface.plane.origin, planeface.plane.vectors[0], planeface.plane.vectors[1] 
        pf1_2d, pf2_2d = Point2D((xmin, ymin)), Point2D((xmin, ymax))
        pf3_2d, pf4_2d = Point2D((xmax, ymin)), Point2D((xmax, ymax))
        pf1, pf2 = pf1_2d.To3D(origin, vx, vy), pf2_2d.To3D(origin, vx, vy) 
        pf3, _ = pf3_2d.To3D(origin, vx, vy), pf4_2d.To3D(origin, vx, vy)
        
        u, v = (pf3-pf1), (pf2-pf1)
        u.Normalize()
        v.Normalize()
        
        w = pf1 - self.center
        
        n1n1, n1u1, n1v1, n1u, n1v = n1.Dot(n1), n1.Dot(u1), n1.Dot(v1), n1.Dot(u), n1.Dot(v)
        u1u1, u1v1, u1u, u1v = u1.Dot(u1), u1.Dot(v1), u1.Dot(u), u1.Dot(v)
        v1v1, v1u, v1v = v1.Dot(v1), v1.Dot(u), v1.Dot(v)
        uu, uv, vv = u.Dot(u), u.Dot(v), v.Dot(v)
        
        w2, wn1, wu1, wv1, wu, wv = w.Dot(w), w.Dot(n1), w.Dot(u1), w.Dot(v1), w.Dot(u), w.Dot(v)
        
        # x = (h, theta, x, y)
        def distance_squared(x):
            return(n1n1*(x[0]**2) + ((math.cos(x[1]))**2)*u1u1*(r**2) + ((math.sin(x[1]))**2)*v1v1*(r**2)
                   + w2 + uu*(x[2]**2) + vv*(x[3]**2) + 2*x[0]*math.cos(x[1])*r*n1u1
                   + 2*x[0]*math.sin(x[1])*r*n1v1 - 2*x[0]*wn1 - 2*x[0]*x[2]*n1u
                   - 2*x[0]*x[3]*n1v + 2*math.sin(x[1])*math.cos(x[1])*u1v1*(r**2)
                   - 2*r*math.cos(x[1])*wu1 - 2*r*x[2]*math.cos(x[1])*u1u 
                   - 2*r*x[3]*math.sin(x[1])*u1v - 2*r*math.sin(x[1])*wv1
                   - 2*r*x[2]*math.sin(x[1])*v1u - 2*r*x[3]*math.sin(x[1])*v1v
                   + 2*x[2]*wu + 2*x[3]*wv + 2*x[2]*x[3]*uv )
        
        x01 = npy.array([(min_h1+max_h1)/2, (min_theta1+max_theta1)/2,
                         (xmax-xmin)/2, (ymax-ymin)/2])

        minimax = [(min_h1, min_theta1, 0, 0), (max_h1, max_theta1, xmax-xmin, ymax-ymin)]
        
        res1 = scp.optimize.least_squares(distance_squared, x01, bounds=minimax)   

        pt1 = Point3D((r*math.cos(res1.x[1]),r*math.sin(res1.x[1]),res1.x[0]))
        p1 = frame1.OldCoordinates(pt1)
        p2 = pf1 + res1.x[2]*u + res1.x[3]*v
        pt1_2d = Point2D((res1.x[1], res1.x[0]))
        pt2_2d = p2.To2D(pf1,u, v)
        
        if not(self.contours2d[0].point_belongs(pt1_2d)) :
            #Find the closest one
            points_contours1 = self.contours2d[0].tessel_points
            
            poly1 = Polygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d, return_other_point=True)
            pt1 = Point3D((r*math.cos(new_pt1_2d.vector[0]),
                           r*math.sin(new_pt1_2d.vector[0]),
                           new_pt1_2d.vector[1]))
            p1 = frame1.OldCoordinates(pt1)
        
        if not(planeface.contours[0].point_belongs(pt2_2d)) :
            #Find the closest one
            d2, new_pt2_2d = planeface.polygon2D.PointBorderDistance(pt2_2d, return_other_point=True)
            
            p2 = new_pt2_2d.To3D(pf1, u, v)
        
        return p1, p2
            
    def minimum_distance(self, other_face, return_points=False) :
        if other_face.__class__ is CylindricalFace3D :
            p1, p2 = self.minimum_distance_points_cyl(other_face)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
        
        if other_face.__class__ is PlaneFace3D : 
            p1, p2 = self.minimum_distance_points_plane(other_face)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
            
        if other_face.__class__ is ToroidalFace3D : 
            p1, p2 = other_face.minimum_distance_points_cyl(self)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
        
        else :
            return NotImplementedError 
            

class ToroidalFace3D (Face3D) :
    """
    :param contours2d: The Tore's contour2D 
    :type contours2d: Contour2D
    :param toroidalsurface3d: Information about the Tore
    :type toroidalsurface3d: ToroidalSurface3D
    :param points: Angle's Tore
    :type points: List of float
    
    :Example: 
        >>> contours2d is rectangular and will create a classic tore with x:2*pi, y:2*pi
        x is for exterior, and y for the circle to revolute
        >>> points = [pi, 2*pi] for an half tore
    """      
    
    def __init__(self, contours2d, toroidalsurface3d, param, name=''):
        
        self.rcenter = toroidalsurface3d.rcenter
        self.rcircle = toroidalsurface3d.rcircle
        self.toroidalsurface3d = toroidalsurface3d 
        
        self.center = self.toroidalsurface3d.frame.origin
        self.normal = self.toroidalsurface3d.frame.w
        vec1, vec2 = self.toroidalsurface3d.frame.u, self.toroidalsurface3d.frame.v
        ptext = self.center + Point3D((self.rcenter*vec1).vector)
        ccircle = ptext - Point3D((self.rcircle*vec1).vector)
        c1 = Arc3D(ptext, self.center+Point3D((self.rcenter*vec1*math.cos(param[0]/2)+self.rcenter*vec2*math.sin(param[0]/2)).vector), self.center+Point3D((self.rcenter*vec1*math.cos(param[0])+self.rcenter*vec2*math.sin(param[0])).vector), self.normal) 
        c2 = Arc3D(ptext, ptext.Rotation(ccircle, vec2, param[1]/2), ptext.Rotation(ccircle, vec2, param[1]), vec2)
            
        edges = [c1, c2]
        ctr = [Contour3D(edges)]
        
        Face3D.__init__(self, ctr)
        self.contours2d = contours2d 
        self.param = param
        self.name = name 
    
    
    @classmethod
    def from_contour3d(cls, contours3d, toroidalsurface3d, name=''):
        """
        :param contours3d: The Tore's contour3D
        :type contours3d: Contour3D
        :param toroidalsurface3d: Information about the Tore
        :type toroidalsurface3d: ToroidalSurface3D
        
        :Example:
            >>> contours3d is [Arc3D, Arc3D, Arc3D], the tore's bones
        """
        frame = toroidalsurface3d.frame
        rcenter, rcircle = toroidalsurface3d.rcenter, toroidalsurface3d.rcircle
        # center = frame.origin
        
        if contours3d[0].__class__ is Point3D : #If it is a complete tore
            angle = 2*math.pi
            pt1, pt2, pt3, pt4 = Point2D((0, 0)), Point2D((0, angle)), Point2D((angle, angle)), Point2D((angle, 0))
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]
            param = [angle, angle]
            
        elif contours3d[0].edges[0].__class__ is Arc3D and contours3d[0].edges[1].__class__ is Arc3D : #Portion of Tore
            theta = contours3d[0].edges[0].angle
            phi1 = contours3d[0].edges[1].angle #arc start
            phi2 = phi1 #contours3d[0].edges[2].angle # if using two different arc at each side of the tore
            #Creation of the window
            pt1, pt2, pt3, pt4 = Point2D((0, 0)), Point2D((0, phi1)), Point2D((theta, phi2)), Point2D((theta, 0))
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]
            param = [theta, phi1]
            
        elif contours3d[0].edges[0].__class__ is Arc3D and contours3d[0].edges[2].__class__ is Arc3D : #if it is a Tore but not linked by Arc3D
            points2d = []
            for edge in contours3d[0].edges :
                new_points = [frame.NewCoordinates(pt) for pt in edge.points]
                points = ToroidalFace3D.points3d_to2d(new_points, rcenter, rcircle)
                points2d.extend(points)
                
            points_cleaned = delete_double_point(points2d)
            all_points = Face3D.range_trigo(points_cleaned)
            primitives = Face3D.create_primitives(all_points)
            contours2d = [Contour2D(primitives)]
            theta = max(pt[0] for pt in contours2d[0].tessel_points) - min(pt[0] for pt in contours2d[0].tessel_points)
            phi = max(pt[1] for pt in contours2d[0].tessel_points) - min(pt[1] for pt in contours2d[0].tessel_points)
            param = [theta, phi]
            
        else:
            contours2d = ToroidalFace3D.contours3d_to2d(contours3d, toroidalsurface3d)
            theta = max(pt[0] for pt in contours2d[0].tessel_points) - min(pt[0] for pt in contours2d[0].tessel_points)
            phi = max(pt[1] for pt in contours2d[0].tessel_points) - min(pt[1] for pt in contours2d[0].tessel_points)
            param = [theta, phi]
        
        return cls(contours2d, toroidalsurface3d, param, name=name)
    
    @classmethod
    def from_arc3d(cls, arc, arcgen):
        """
        :param arc: The arc which is extruded by the arcgen
        :type arc: Arc3D/2D, Circle3D/2D
        :param arcgen: The Arc generator
        :type arcgen: Arc3D/2D, Circle3D/2D
        """
        rcircle = arc.radius
        rcenter = arcgen.radius
        
        center = arcgen.center
        normal = arcgen.normal
        normal.Normalize()
        
        center1 = arcgen.start
        u = Vector3D((center1 - center).vector)
        u.Normalize()
        v = normal.Cross(u)
        # print('arc angle', arc.angle)
        offset1 = 0
        if arcgen.__class__.__name__ == 'Circle3D' or arcgen.__class__.__name__ == 'Circle2D' : 
            theta = 2*math.pi
       
        else : #Offset for the Arcgen
            point_last = arcgen.end
            if point_last.__class__ is Point3D :
                point_last = arcgen.end.To2D(center, u, v)
            x1, y1 = point_last.vector[0], point_last.vector[1]
            
            theta = math.atan2(y1, x1)
            if theta < 0 or math.isclose(theta, 0, abs_tol=1e-9):
                if arcgen.angle>math.pi :
                    theta += 2*math.pi
                else :
                    offset1 = theta
                    theta = -theta
       
        offset2 = 0
        
        if arc.__class__.__name__ == 'Circle3D' or arc.__class__.__name__ == 'Circle2D' : 
            phi = 2*math.pi
        
        else : # Offset for the Arc 
            
            center_generated = arc.center
            n = v
            if arc.start.__class__ is Point2D :
                u_g2d = Vector2D((arc.start - center_generated).vector)
                u_g = u_g2d.To3D(center, normal, u)
                u_g.Normalize()
            else :
                u_g = Vector3D((arc.start - center_generated).vector)
                u_g.Normalize()
            v_g = n.Cross(u_g)
            v_g.Normalize()
            
            _, last_generated, c2d = arc.end, arc.start, center_generated
            if last_generated.__class__ is Point3D :
                # first_generated = arc.start.To2D(center_generated, u_g, v_g)
                c2d = center.To2D(center_generated, u_g, v_g)
                last_generated = last_generated.To2D(center_generated, u_g, v_g)
            x2, y2 = last_generated.vector[0], last_generated.vector[1]
            phi = math.atan2(y2, x2)
            # Calculate angle between first point of arcgen and arc
          
            angle_offset = math.atan2(c2d.vector[1], c2d.vector[0]) 
            if arc.__class__.__name__ == 'Arc2D' :
                offset2 = 0
                phi = arc.angle
            else :
                if n == -arc.normal :
                    offset2 += -phi+math.pi+angle_offset
                else : 
                    offset2 += -phi-math.pi+angle_offset

        frame3d = Frame3D(center, u, v, normal)
        # print('rcenter rcircle')
        # print(rcenter, rcircle)
        toroidalsurface3d = ToroidalSurface3D(frame3d, rcenter, rcircle)
        # print('tp', theta, phi)
        pt1, pt2, pt3, pt4 = Point2D((offset1, offset2)), Point2D((offset1, phi+offset2)), Point2D((theta+offset1, phi+offset2)), Point2D((theta+offset1, offset2))
        seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
        edges = [seg1, seg2, seg3, seg4]
        contours2d =  [Contour2D(edges)]
        param = [theta, phi]
        print('theta phi', theta, phi)
        
        return cls(contours2d, toroidalsurface3d, param, name='')
    
    def points_resolution(self, line, pos, resolution) : #With a resolution wished
        points = []
        points.append(line.points[0])
        limit = line.points[1].vector[pos]
        start = line.points[0].vector[pos]
        vec = [0,0]
        vec[pos] = start
        echelon = [line.points[0].vector[0] - vec[0], line.points[0].vector[1] - vec[1]]
        flag = start + resolution
        while flag < limit :
            echelon[pos] = flag
            flag += resolution
            points.append(Point2D(echelon))
        points.append(line.points[1])
        return points
    
    def points2d_to3d(self, points2d, rcenter, rcircle, frame3d) :
        # source wikipedia Tore
        points3D = []
        R = rcenter
        for pt in points2d :
            phi, theta = pt[1], pt[0] 
            x = (R+rcircle*math.cos(phi))*math.cos(theta)
            y = (R+rcircle*math.cos(phi))*math.sin(theta)
            z = rcircle*math.sin(phi)
            points3D.append(Point3D([x,y,z]))
        Points_3D = [frame3d.OldCoordinates(point) for point in points3D]
        return Points_3D
    
    def points3d_to2d(points3d, rcenter, rcircle) :
        points_2D = []
        R, r = rcenter, rcircle 
        for pt in points3d :
            x, y, z = pt[0], pt[1], pt[2]
            if z<-r :
                z = -r
            elif z>r :
                z = r
                
            zr = z/r
            phi = math.asin(zr)
            
            u = R + math.sqrt((r**2)-(z**2))
            u1, u2 = round(x/u, 5), round(y/u, 5)
            #cos(theta)=u1, sin(theta)=u2 :
            theta = sin_cos_angle(u1, u2)
                    
            points_2D.append(Point2D([theta, phi]))
            
        i0, i2pi, iangle, ih = 0, 0, 0, 0
        for enum,point in enumerate(points_2D) :
            if enum == 0 :
                h = point[1]
                angle = point[0]                
            if math.isclose(point[0], 0, abs_tol=1e-2) :
                i0 += 1
            elif math.isclose(point[0], 2*math.pi, abs_tol=1e-3) :
                i2pi += 1
            elif math.isclose(point[0], angle, abs_tol=1e-4) :
                iangle += 1
            elif math.isclose(point[1], h, abs_tol=1e-4) :
                ih += 1
                
        if i2pi/len(points_2D)>0.5 or i0/len(points_2D)>0.5:
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((2*math.pi, point.vector[1])))
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        # elif i0/len(points_2D)>0.5 :
        #     new_points2d = []
        #     for point in points_2D :
        #         new_points2d.append(Point2D((2*math.pi, point.vector[1]))) 
        #     new_points2d.sort(key=lambda pt: pt[1])
        #     points_2D = [new_points2d[0], new_points2d[-1]]
        elif iangle/len(points_2D)>0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((angle, point.vector[1]))) 
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif ih/len(points_2D)>0.5 :
            new_points2d = [pt.copy() for pt in points_2D]
            new_points2d.sort(key=lambda pt: pt[0])
            if i2pi == 2 or i0 == 2 or (i0 >= 1 and iangle == 2) or (i2pi >= 1 and iangle == 2): 
                pt1, pt2 = new_points2d[0], new_points2d[-1]
                pt1.vector[0] = 0-0.0000001
                pt2.vector[0] = 2*math.pi+0.0000001
                points_2D = [pt1, pt2]
            else :
                points_2D = [new_points2d[0], new_points2d[-1]]
        
        plus_7pi4, moins_7pi4 = [], []
        for enum, pt in enumerate(points_2D) :
            if pt.vector[0]>7*math.pi/4 :
                plus_7pi4.append(enum)
            elif pt.vector[0]<math.pi/4 :
                moins_7pi4.append(enum)
            
        if len(moins_7pi4) > 3*len(plus_7pi4) :
            for pos in plus_7pi4 :
                new_pt = points_2D[pos].copy() - Point2D((2*math.pi, 0))
                points_2D[pos] = new_pt  
                
        return points_2D
    
    def contours3d_to2d(contours3d, toroidalsurface3d) :
        frame = toroidalsurface3d.frame
        n = frame.w
        # center = frame.origin
        rcenter, rcircle = toroidalsurface3d.rcenter, toroidalsurface3d.rcircle
        
        # print('contours3d[0].edges', contours3d[0].edges)
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # [edge.MPLPlot(ax=ax) for edge in contours3d[0].edges]
        
        primitives, start_end, all_points = [], [], []
        for edge in contours3d[0].edges :
            new_points = [frame.NewCoordinates(pt) for pt in edge.points]
            if edge.__class__ is Arc3D :
                if edge.normal == n or edge.normal == -n:
                    start2d, end2d = ToroidalFace3D.points3d_to2d(new_points, rcenter, rcircle)
                    # if math.isclose(edge.angle, 2*math.pi, abs_tol=1e-6):
                    #     if start2d == end2d :
                    #         if math.isclose(start2d.vector[0], 2*math.pi, abs_tol = 1e-6) :
                    #             end2d = end2d - Point2D((2*math.pi, 0))
                    #         else :
                    #             end2d = end2d + Point2D((2*math.pi, 0))
                    
                    ######################New
                    angle2d = abs(end2d[0]-start2d[0])
                    if math.isclose(edge.angle, 2*math.pi, abs_tol=1e-6):
                        if start2d == end2d :
                            if math.isclose(start2d.vector[0], 2*math.pi, abs_tol = 1e-6) :
                                end2d = end2d - Point2D((2*math.pi, 0))
                            else :
                                end2d = end2d + Point2D((2*math.pi, 0))
                    elif not(math.isclose(edge.angle, angle2d, abs_tol=1e-2)) :
                        # if math.isclose(angle2d, 2*math.pi, abs_tol=1e-2) :
                        if start2d[0] < end2d[0] :
                            end2d = start2d + Point2D((edge.angle,0))
                        else :
                            end2d = start2d - Point2D((edge.angle,0))
                    #####################
                    
                    ls_toadd = LineSegment2D(start2d, end2d)
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        primitives.append([ls_toadd])
                        all_points.extend(ls_toadd.points)
                        start_end.append(ls_toadd.points)
                    
                else :
                    points2d = ToroidalFace3D.points3d_to2d(new_points, rcenter, rcircle)
                    lines = []
                    for k in range(0, len(points2d)-1) :
                        lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                    points, prim_list = [], []
                    for ls_toadd in lines :        
                        same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                        if same is False :
                            prim_list.append(ls_toadd)
                            points.extend(ls_toadd.points)
                    if len(points) > 0 : 
                        all_points.extend(points)
                        primitives.append(prim_list)
                        start_end.append([points[0], points[-1]])
                    
            elif edge.__class__ is LineSegment3D :
                start2d, end2d = ToroidalFace3D.points3d_to2d(new_points, rcenter, rcircle)
                ls_toadd = LineSegment2D(start2d, end2d)
                same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                if same is False :
                    primitives.append([ls_toadd])
                    all_points.extend(ls_toadd.points)
                    start_end.append(ls_toadd.points)
            
            else :
                points2d = ToroidalFace3D.points3d_to2d(new_points, rcenter, rcircle)
                lines = []
                for k in range(0, len(points2d)-1) :
                    lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                points, prim_list = [], []
                for ls_toadd in lines :        
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        prim_list.append(ls_toadd)
                        points.extend(ls_toadd.points)
                if len(points) > 0 : 
                    all_points.extend(points)
                    primitives.append(prim_list)
                    start_end.append([points[0], points[-1]])
        
        # fig, ax = plt.subplots() 
        # [pt.MPLPlot(ax=ax) for pt in all_points]
        # for list_prim in primitives :
        #     fig, ax = plt.subplots() 
        #     [prim.MPLPlot(ax=ax) for prim in list_prim]  
        
        # raise NotImplementedError
        
        points_se, primitives_se = [], []
        for double in start_end : 
            primitives_se.append(LineSegment2D(double[0], double[1]))
            points_se.extend(double)
        poly_se = Polygon2D(points_se)
            
        xmax, xmin = max(pt[0] for pt in points_se), min(pt[0] for pt in points_se)
        ymax, ymin = max(pt[1] for pt in points_se), min(pt[1] for pt in points_se) 
        pt1, pt2, pt3, pt4 = Point2D((xmin, ymin)), Point2D((xmin, ymax)), Point2D((xmax, ymin)), Point2D((xmax, ymax))
        diag1, diag2 = LineSegment2D(pt1, pt4), LineSegment2D(pt2, pt3)
        diag1_cut, diag2_cut = [], []
        diag1_pointcut, diag2_pointcut = [], []
        for enum,l in enumerate(primitives_se) :
            cut1 = diag1.line_intersection(l)
            cut2 = diag2.line_intersection(l)
            if cut1 is not None : 
                diag1_cut.append(enum)
                diag1_pointcut.append(cut1)
            if cut2 is not None : 
                diag2_cut.append(enum)
                diag2_pointcut.append(cut2)
        
        points_common = []
        for enum1,pos1 in enumerate(diag1_cut) :
            for enum2,pos2 in enumerate(diag2_cut) :
                if pos1 == pos2 :
                    points_common.append(primitives_se[pos1].points)
        # fig, ax = plt.subplots() 
        # [l.MPLPlot(ax=ax) for l in lines]
        # [pt.MPLPlot(ax=ax, color='g') for pt in points_se]
        # pt1.MPLPlot(ax=ax, color='r')
        # pt2.MPLPlot(ax=ax, color='r')            
        # pt3.MPLPlot(ax=ax, color='r')
        # pt4.MPLPlot(ax=ax, color='r')   
        # diag1.MPLPlot(ax=ax)
        # diag2.MPLPlot(ax=ax)
        # for couple in points_common :
        #     [pt.MPLPlot(ax=ax, color='b') for pt in couple]        
        
        if len(points_common) >= 1 :
            solve = False
            for couple in points_common :
                if solve :
                    break
                check1, check2 = poly_se.PointBelongs(couple[0]), poly_se.PointBelongs(couple[1])
                start, end = couple[0].vector[0], couple[1].vector[0]
                if math.isclose(start, end, abs_tol = 5e-2) :
                    intersect = min(start, end)
                    if math.isclose(intersect, math.pi, abs_tol = 5e-2) :
                        # all_points = check_singularity(all_points)
                        # intersect = 0
                        
                        ##################### NEW
                        points_sing = check_singularity(all_points)
                        pt0, pt2pi = 0, 0
                        for pt in points_sing :
                            if math.isclose(pt.vector[0], 0, abs_tol=1e-2):
                                pt0 += 1
                            elif math.isclose(pt.vector[0], 2*math.pi, abs_tol=1e-2):
                                pt2pi += 1
                        points_sing.sort(key=lambda pt: pt[1])
                        points_sing.sort(key=lambda pt: pt[0])
                        if pt2pi !=0 and pt0 == 0 :
                            points = [pt.copy() for pt in points_sing[::-1]]
                            points_sing = points
                        points_range = Face3D.range_closest(points_sing)
                        all_points = delete_double_point(points_range)
                        break
                        #######################
                    
                    # if math.isclose(intersect, 0, abs_tol = 1e-6) or math.isclose(intersect, 2*math.pi, abs_tol = 1e-6) or (not check1 or not check2):
                    elif math.isclose(intersect, 0, abs_tol = 1e-6) or math.isclose(intersect, 2*math.pi, abs_tol = 1e-6) or (not check1 or not check2):
                        all_points = check_singularity(all_points)
                        
                        points_cleaned = delete_double_point(all_points)
                        all_points = [pt.copy() for pt in points_cleaned]
                        all_points.sort(key=lambda pt: pt[0])
                        d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                        if d2 < d1 :
                            last = all_points[-1].copy()
                            all_points[-1] = all_points[-2].copy()
                            all_points[-2] = last
                        break
                    else :
                        points = []
                        for list_prim in primitives :
                            for k,prim in enumerate(list_prim) :
                                new_list_points = []
                                change = 0
                                for pt in prim.points :
                                    if pt[0] < intersect :
                                        change += 1 
                                        if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                            new_list_points.append(Point2D((intersect + 2*math.pi, pt[1])))
                                        else :
                                            new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    elif math.isclose(pt[0], intersect, abs_tol=1e-1) :
                                        change += 1
                                        new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    else :
                                        new_list_points.append(pt)
                                if change > 0 :
                                    points.extend(new_list_points)
                                    # list_prim[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                                else :
                                    points.extend(prim.points)
                                    continue
                        points_cleaned = delete_double_point(points)
                        all_points = Face3D.range_trigo(points_cleaned)
                    solve = True
                else :
                    points_cleaned = delete_double_point(all_points)
                    all_points = [pt.copy() for pt in points_cleaned]
                    all_points.sort(key=lambda pt: pt[0])
                    d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                    if d2 < d1 :
                        last = all_points[-1].copy()
                        all_points[-1] = all_points[-2].copy()
                        all_points[-2] = last
        else :
            points_cleaned = delete_double_point(all_points)
            all_points = [pt.copy() for pt in points_cleaned]
            all_points.sort(key=lambda pt: pt[0])
            d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
            if d2 < d1 :
                last = all_points[-1].copy()
                all_points[-1] = all_points[-2].copy()
                all_points[-2] = last 
        
        primitives = Face3D.create_primitives(all_points)
        
        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax, color='g') for pt in all_points]
        # all_points[0].MPLPlot(ax=ax, color='m')
        # all_points[1].MPLPlot(ax=ax, color='r')
        # all_points[-1].MPLPlot(ax=ax)
        # all_points[-2].MPLPlot(ax=ax, color='b')
        # [p.MPLPlot(ax=ax) for p in primitives]
        
        l_vert = LineSegment2D((pt2+pt4)/2, (pt1+pt3)/2)
        solve = False
        for prim in primitives :
            if solve :
                break
            intersect = prim.line_intersection(l_vert)
            if intersect is not None : 
                x_intersect = intersect.vector[0]
                y_intersect = intersect.vector[1]
                value1, value2 = ymax-0.2*(ymax-ymin), ymin+0.2*(ymax-ymin)
                if y_intersect < max(value1, value2) and y_intersect > min(value1, value2) :
                    points = []
                    for k, prim in enumerate(primitives) :
                        new_list_points, change = [], 0
                        for pt in prim.points :
                            if pt[0] < x_intersect :
                                change += 1
                                if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                    new_list_points.append(Point2D((x_intersect + 2*math.pi, pt[1])))
                                else :
                                    new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                            else :
                                new_list_points.append(pt)
                        if change > 0:
                            points.extend(new_list_points)
                            primitives[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                        else :
                            points.extend(prim.points)
                            continue
                    solve = True
                    points_cleaned = delete_double_point(points)
                    all_points = Face3D.range_trigo(points_cleaned)
                    primitives = Face3D.create_primitives(all_points)
                    
        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax) for pt in all_points]
        # [p.MPLPlot(ax=ax) for p in primitives]
        # intersect.MPLPlot(ax=ax, color='r')
        # l_vert.MPLPlot(ax=ax)
        contour2d = [Contour2D(primitives)]
        return contour2d
    
    def triangulation(self, resolution=30):
        rcenter = self.rcenter
        rcircle = self.rcircle
        
        frame3d = self.toroidalsurface3d.frame
        
        centerota = Point2D((0,0))
        
        angle_theta = self.param[0]
        angle_phi = self.param[1]
        pas_theta = 2*math.pi/30 #Step of 12 degrees
        pas_phi = pas_theta
        
        resolution_theta = abs(int(angle_theta/pas_theta))
        resolution_phi = abs(int(angle_phi/pas_phi))
        
        if resolution_phi < 5 :
            resolution_phi = 5
        if resolution_theta < 5 :
            resolution_theta = 5
        
        ctr_pt1 = self.cut_contours(self.contours2d, resolution_theta)
        
        all_contours_points = []
        for listpt in ctr_pt1 :
            bandept = []
            for pt in listpt :
                bandept.append(pt.Rotation(centerota, -math.pi/2))
            edges = []
            for k in range (0,len(bandept)) :
                if k == len(bandept)-1 :
                    edges.append(LineSegment2D(bandept[k], bandept[0]))
                else :
                    edges.append(LineSegment2D(bandept[k], bandept[k+1]))
            ctr2d = [Contour2D(edges)]
            all_contours_points.extend(self.cut_contours(ctr2d, resolution_phi))
        
        pts_frame1 = []
        for listpt in all_contours_points :
            for pt in listpt : 
                pts_frame1.append(pt.Rotation(centerota, math.pi/2))
        
        all_points = self.delete_double(pts_frame1) # All points necessary to triangulate
        all_points.sort(key=lambda pt: pt[0]) 
        ptvert, pts = [], []
        for k in range(0, resolution_theta +1) :
            ptvert.append([point for point in all_points[k*(resolution_phi+1):(k+1)*(resolution_phi+1)]])
            ptvert[k].sort(key=lambda pt: pt[1])
            pts.extend(ptvert[k])
        
        # A = dict(vertices=[pt.vector for pt in pts]) #in all_points
        # B = triangle.triangulate(A, 'cp')
        # Triangles = list(B['triangles'])
        # triangle.compare(plt, A, B)
        # print(B['triangles'].shape)
        # print(B)
        # pts3d = self.points2d_to3d(pts, rcenter, rcircle, frame3d)
        # pt3d, tangle = delete_double_pos(pts3d, [Triangles])
        
        Triangles, ts = [], []
        if math.isclose(angle_phi, 2*math.pi, abs_tol=1e-6) or math.isclose(angle_theta, 2*math.pi, abs_tol=1e-6):
            step = resolution_phi
            for k in range (0,len(pts)-resolution_phi-2) :
                vertices, segments = [], []
                
                if k%step == 0 and k!=0:
                    step += resolution_phi+1
                    continue
                
                listpt = [pts[k], pts[k+1], pts[k+1+resolution_phi+1], pts[k+resolution_phi+1]]
                # print('>>>>>>listpt', listpt)
                listindice = [k, k+1, k+1+resolution_phi+1, k+resolution_phi+1]
                for i, pt in enumerate(listpt):
                    vertices.append(pt.vector)
                    segments.append([i,i+1])
                
                segments[-1]=(len(listpt)-1,0) 
                tri = {'vertices': vertices, 'segments': segments}
                t = triangle.triangulate(tri, 'p')
                if 'triangles' in t:
                    triangles = t['triangles'].tolist()
                    for n,tri in enumerate(triangles):
                        # print('tri', tri)
                        for i in range (0,3):
                            tri[i]=listindice[tri[i]]
                    Triangles.append(triangles)        
                else:
                    Triangles.append(None)
                ts.append(t)
                
                pts3d = self.points2d_to3d(pts, rcenter, rcircle, frame3d)
                pt3d, tangle = delete_double_pos(pts3d, Triangles)
        else :
            A = dict(vertices=[pt.vector for pt in pts]) #in all_points
            B = triangle.triangulate(A, 'cp')
            Triangles = list(B['triangles'])
            # triangle.compare(plt, A, B)
            # print(B['triangles'].shape)
            # print(B)
            pts3d = self.points2d_to3d(pts, rcenter, rcircle, frame3d)
            pt3d, tangle = delete_double_pos(pts3d, [Triangles]) 
        
        
        
        return pt3d, tangle
    
    def frame_mapping(self, frame, side, copy=True) :
        if copy:
            new_toroidalsurface3d = ToroidalSurface3D.frame_mapping(frame, side, copy)
            return ToroidalFace3D(self.contours2d, new_toroidalsurface3d, points=self.points, name=self.name)
        else:
            self.toroidalsurface3d.frame_mapping(frame, side, copy=False)
            
    def minimum_maximum_tore(self, contour2d) :
        points = contour2d.tessel_points
        
        min_phi, min_theta = min([pt[1] for pt in points]), min([pt[0] for pt in points])
        max_phi, max_theta = max([pt[1] for pt in points]), max([pt[0] for pt in points])
        return min_phi, min_theta, max_phi, max_theta
    
    def minimum_distance_points_tore(self, other_tore) :
        R1, r1, R2, r2 = self.rcenter, self.rcircle, other_tore.rcenter, other_tore.rcircle
        
        min_phi1, min_theta1, max_phi1, max_theta1 = self.minimum_maximum_tore(self.contours2d[0])
        
        # start1 = self.start
        n1 = self.normal
        u1 = self.toroidalsurface3d.frame.u
        v1 = self.toroidalsurface3d.frame.v
        frame1 = Frame3D(self.center, u1, v1, n1)
        # start1 = self.points2d_to3d([[min_theta1, min_phi1]], R1, r1, frame1)
        
        min_phi2, min_theta2, max_phi2, max_theta2 = self.minimum_maximum_tore(other_tore.contours2d[0])
             
        # start2 = other_tore.start
        n2 = other_tore.normal
        u2 = other_tore.toroidalsurface3d.frame.u
        v2 = other_tore.toroidalsurface3d.frame.v
        frame2 = Frame3D(other_tore.center, u2, v2, n2)
        # start2 = other_tore.points2d_to3d([[min_theta2, min_phi2]], R2, r2, frame2)
        
        w = other_tore.center - self.center
        

        n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.Dot(n1), n1.Dot(u1), n1.Dot(v1), n1.Dot(n2), n1.Dot(u2), n1.Dot(v2)
        u1u1, u1v1, u1n2, u1u2, u1v2 = u1.Dot(u1), u1.Dot(v1), u1.Dot(n2), u1.Dot(u2), u1.Dot(v2)
        v1v1, v1n2, v1u2, v1v2 = v1.Dot(v1), v1.Dot(n2), v1.Dot(u2), v1.Dot(v2)
        n2n2, n2u2, n2v2 = n2.Dot(n2), n2.Dot(u2), n2.Dot(v2)
        u2u2, u2v2, v2v2 = u2.Dot(u2), u2.Dot(v2), v2.Dot(v2)
        
        w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.Dot(w), w.Dot(n1), w.Dot(u1), w.Dot(v1), w.Dot(n2), w.Dot(u2), w.Dot(v2)
        
        # x = (phi1, theta1, phi2, theta2)
        def distance_squared(x):
            return (u1u1*(((R1+r1*math.cos(x[0]))*math.cos(x[1]))**2)
                    + v1v1*(((R1+r1*math.cos(x[0]))*math.sin(x[1]))**2)
                    + n1n1*((math.sin(x[0]))**2)*(r1**2) + w2
                    + u2u2*(((R2+r2*math.cos(x[2]))*math.cos(x[3]))**2)
                    + v2v2*(((R2+r2*math.cos(x[2]))*math.sin(x[3]))**2)
                    + n2n2*((math.sin(x[2]))**2)*(r2**2)
                    + 2*u1v1*math.cos(x[1])*math.sin(x[1])*((R1+r1*math.cos(x[0]))**2)
                    + 2*(R1+r1*math.cos(x[0]))*math.cos(x[1])*r1*math.sin(x[0])*n1u1
                    - 2*(R1+r1*math.cos(x[0]))*math.cos(x[1])*wu1
                    - 2*(R1+r1*math.cos(x[0]))*(R2+r2*math.cos(x[2]))*math.cos(x[1])*math.cos(x[3])*u1u2
                    - 2*(R1+r1*math.cos(x[0]))*(R2+r2*math.cos(x[2]))*math.cos(x[1])*math.sin(x[3])*u1v2
                    - 2*(R1+r1*math.cos(x[0]))*math.cos(x[1])*r2*math.sin(x[2])*u1n2
                    + 2*(R1+r1*math.cos(x[0]))*math.sin(x[1])*r1*math.sin(x[0])*n1v1
                    - 2*(R1+r1*math.cos(x[0]))*math.sin(x[1])*wv1
                    - 2*(R1+r1*math.cos(x[0]))*(R2+r2*math.cos(x[2]))*math.sin(x[1])*math.cos(x[3])*v1u2
                    - 2*(R1+r1*math.cos(x[0]))*(R2+r2*math.cos(x[2]))*math.sin(x[1])*math.sin(x[3])*v1v2
                    - 2*(R1+r1*math.cos(x[0]))*math.sin(x[1])*r2*math.sin(x[2])*v1n2
                    - 2*r1*math.sin(x[0])*wn1 
                    - 2*r1*math.sin(x[0])*(R2+r2*math.cos(x[2]))*math.cos(x[3])*n1u2
                    - 2*r1*math.sin(x[0])*(R2+r2*math.cos(x[2]))*math.sin(x[3])*n1v2
                    - 2*r1*r2*math.sin(x[0])*math.sin(x[2])*n1n2
                    + 2*(R2+r2*math.cos(x[2]))*math.cos(x[3])*wu2
                    + 2*(R2+r2*math.cos(x[2]))*math.sin(x[3])*wv2
                    + 2*r2*math.sin(x[2])*wn2 
                    + 2*u2v2*math.cos(x[3])*math.sin(x[3])*((R2+r2*math.cos(x[2]))**2)
                    + 2*math.cos(x[3])*(R2+r2*math.cos(x[2]))*r2*math.sin(x[2])*n2u2
                    + 2*math.sin(x[3])*(R2+r2*math.cos(x[2]))*r2*math.sin(x[2])*n2v2)
                    
        
        x01 = npy.array([(min_phi1+max_phi1)/2, (min_theta1+max_theta1)/2,
                         (min_phi2+max_phi2)/2, (min_theta2+max_theta2)/2 ])
        x02 = npy.array([min_phi1, min_theta1,
                         min_phi2, min_theta2])
        x03 = npy.array([max_phi1, max_theta1,
                         max_phi2, max_theta2])
        
        minimax = [(min_phi1, min_theta1, min_phi2, min_theta2), (max_phi1, max_theta1, max_phi2, max_theta2)]
        
        res1 = scp.optimize.least_squares(distance_squared, x01, bounds=minimax)
        res2 = scp.optimize.least_squares(distance_squared, x02, bounds=minimax)
        res3 = scp.optimize.least_squares(distance_squared, x03, bounds=minimax)
        
        # frame1, frame2 = Frame3D(self.center, u1, v1, n1), Frame3D(other_tore.center, u2, v2, n2)
        pt1 = self.points2d_to3d([[res1.x[1], res1.x[0]]], R1, r1, frame1)
        pt2 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R2, r2, frame2)
        p1, p2 = pt1[0], pt2[0]
        d = p1.point_distance(p2)
        result = res1
        
        res = [res2, res3]
        for couple in res :
            ptest1 = self.points2d_to3d([[couple.x[1], couple.x[0]]], R1, r1, frame1)
            ptest2 = self.points2d_to3d([[couple.x[3], couple.x[2]]], R2, r2, frame2)
            dtest = ptest1[0].point_distance(ptest2[0])
            if dtest < d :
                result = couple
                p1, p2 = ptest1[0], ptest2[0]
        
        pt1_2d, pt2_2d = Point2D((result.x[1], result.x[0])), Point2D((result.x[3], result.x[2]))
               
        if not(self.contours2d[0].point_belongs(pt1_2d)) :
            #Find the closest one
            points_contours1 = self.contours2d[0].tessel_points
            
            poly1 = Polygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d, return_other_point=True)
                    
            pt1 = self.points2d_to3d([new_pt1_2d], R1, r1, frame1)
            p1 = pt1[0]
        
        
        if not(other_tore.contours2d[0].point_belongs(pt2_2d)) :
            #Find the closest one
            points_contours2 = other_tore.contours2d[0].tessel_points
            
            poly2 = Polygon2D(points_contours2)
            d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d, return_other_point=True)
                    
            pt2 = self.points2d_to3d([new_pt2_2d], R2, r2, frame2)
            p2 = pt2[0]
            
        return p1, p2
    
    def minimum_distance_points_cyl(self, cyl) :
        R2, r2, r = self.rcenter, self.rcircle, cyl.radius
        
        min_h, min_theta, max_h, max_theta = cyl.minimum_maximum(cyl.contours2d[0], r)
             
        n1 = cyl.normal
        u1 = cyl.cylindricalsurface3d.frame.u
        v1 = cyl.cylindricalsurface3d.frame.v
        frame1 = Frame3D(cyl.center, u1, v1, n1)
        # st1 = Point3D((r*math.cos(min_theta), r*math.sin(min_theta), min_h))
        # start1 = frame1.OldCoordinates(st1)
        
        min_phi2, min_theta2, max_phi2, max_theta2 = self.minimum_maximum_tore(self.contours2d[0])
        
        n2 = self.normal
        u2 = self.toroidalsurface3d.frame.u
        v2 = self.toroidalsurface3d.frame.v
        frame2 = Frame3D(self.center, u2, v2, n2)
        # start2 = self.points2d_to3d([[min_theta2, min_phi2]], R2, r2, frame2)
        
        w = self.center - cyl.center
        

        n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.Dot(n1), n1.Dot(u1), n1.Dot(v1), n1.Dot(n2), n1.Dot(u2), n1.Dot(v2)
        u1u1, u1v1, u1n2, u1u2, u1v2 = u1.Dot(u1), u1.Dot(v1), u1.Dot(n2), u1.Dot(u2), u1.Dot(v2)
        v1v1, v1n2, v1u2, v1v2 = v1.Dot(v1), v1.Dot(n2), v1.Dot(u2), v1.Dot(v2)
        n2n2, n2u2, n2v2 = n2.Dot(n2), n2.Dot(u2), n2.Dot(v2)
        u2u2, u2v2, v2v2 = u2.Dot(u2), u2.Dot(v2), v2.Dot(v2)
        
        w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.Dot(w), w.Dot(n1), w.Dot(u1), w.Dot(v1), w.Dot(n2), w.Dot(u2), w.Dot(v2)
        
        # x = (theta, h, phi2, theta2)
        def distance_squared(x):
            return (u1u1*((math.cos(x[0])*r)**2) + v1v1*((math.sin(x[0])*r)**2)
                    + n1n1*(x[1]**2) + w2 
                    + u2u2*(((R2+r2*math.cos(x[2]))*math.cos(x[3]))**2)
                    + v2v2*(((R2+r2*math.cos(x[2]))*math.sin(x[3]))**2)
                    + n2n2*((math.sin(x[2]))**2)*(r2**2)
                    + 2*u1v1*math.cos(x[0])*math.sin(x[0])*(r**2)
                    + 2*r*math.cos(x[0])*x[1]*n1u1 - 2*r*math.cos(x[0])*wu1
                    - 2*r*math.cos(x[0])*(R2+r2*math.cos(x[2]))*math.cos(x[3])*u1u2
                    - 2*r*math.cos(x[0])*(R2+r2*math.cos(x[2]))*math.sin(x[3])*u1v2
                    - 2*r*math.cos(x[0])*r2*math.sin(x[2])*u1n2 
                    + 2*r*math.sin(x[0])*x[1]*n1v1- 2*r*math.sin(x[0])*wv1
                    - 2*r*math.sin(x[0])*(R2+r2*math.cos(x[2]))*math.cos(x[3])*v1u2
                    - 2*r*math.sin(x[0])*(R2+r2*math.cos(x[2]))*math.sin(x[3])*v1v2
                    - 2*r*math.sin(x[0])*r2*math.sin(x[2])*v1n2 - 2*x[1]*wn1
                    - 2*x[1]*(R2+r2*math.cos(x[2]))*math.cos(x[3])*n1u2
                    - 2*x[1]*(R2+r2*math.cos(x[2]))*math.sin(x[3])*n1v2
                    - 2*x[1]*r2*math.sin(x[2])*n1n2
                    + 2*(R2+r2*math.cos(x[2]))*math.cos(x[3])*wu2
                    + 2*(R2+r2*math.cos(x[2]))*math.sin(x[3])*wv2
                    + 2*r2*math.sin(x[2])*wn2 
                    + 2*u2v2*math.cos(x[3])*math.sin(x[3])*((R2+r2*math.cos(x[2]))**2)
                    + 2*math.cos(x[3])*(R2+r2*math.cos(x[2]))*r2*math.sin(x[2])*n2u2
                    + 2*math.sin(x[3])*(R2+r2*math.cos(x[2]))*r2*math.sin(x[2])*n2v2)
                    
        
        x01 = npy.array([(min_theta+max_theta)/2, (min_h+max_h)/2,
                         (min_phi2+max_phi2)/2, (min_theta2+max_theta2)/2 ])
        x02 = npy.array([min_theta, min_h,
                         min_phi2, min_theta2])
        x03 = npy.array([max_theta, max_h,
                         max_phi2, max_theta2])
        
        minimax = [(min_theta, min_h, min_phi2, min_theta2), (max_theta, max_h, max_phi2, max_theta2)]
        
        res1 = scp.optimize.least_squares(distance_squared, x01, bounds=minimax)
        res2 = scp.optimize.least_squares(distance_squared, x02, bounds=minimax)
        res3 = scp.optimize.least_squares(distance_squared, x03, bounds=minimax)
        
        pt1 = Point3D((r*math.cos(res1.x[0]),r*math.sin(res1.x[0]),res1.x[1]))
        p1 = frame1.OldCoordinates(pt1)
        pt2 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R2, r2, frame2)
        p2 = pt2[0]
        d = p1.point_distance(p2)
        result = res1
        
        res = [res2, res3]
        for couple in res :
            pttest1 = Point3D((r*math.cos(couple.x[0]),r*math.sin(couple.x[0]),couple.x[1]))
            ptest1 = frame1.OldCoordinates(pttest1)
            ptest2 = self.points2d_to3d([[couple.x[3], couple.x[2]]], R2, r2, frame2)
            dtest = ptest1.point_distance(ptest2[0])
            if dtest < d :
                result = couple
                p1, p2 = ptest1, ptest2[0]
        
        # pt1_2d, pt2_2d = Point2D((result.x[0]*r, result.x[1])), Point2D((result.x[3], result.x[2]))
        pt1_2d, pt2_2d = Point2D((result.x[0], result.x[1])), Point2D((result.x[3], result.x[2]))
        
               
        # fig, ax = plt.subplots()
        # [prim.MPLPlot(ax=ax) for prim in cyl.contours2d[0].primitives]
        # Point2D((min_theta*r, min_h)).MPLPlot(ax=ax)
        # Point2D((min_theta*r, max_h)).MPLPlot(ax=ax)
        # Point2D((max_theta*r, min_h)).MPLPlot(ax=ax)
        # Point2D((max_theta*r, max_h)).MPLPlot(ax=ax)
        # pt1_2d.MPLPlot(ax=ax, color='r')
        
        # fig, ax = plt.subplots()
        # [prim.MPLPlot(ax=ax) for prim in self.contours2d[0].primitives]
        # Point2D((min_theta2, min_phi2)).MPLPlot(ax=ax)
        # Point2D((min_theta2, max_phi2)).MPLPlot(ax=ax)
        # Point2D((max_theta2, min_phi2)).MPLPlot(ax=ax)
        # Point2D((max_theta2, max_phi2)).MPLPlot(ax=ax)
        # pt2_2d.MPLPlot(ax=ax, color='g')
        
        if not(self.contours2d[0].point_belongs(pt2_2d)) :
            #Find the closest one
            points_contours2 = self.contours2d[0].tessel_points
            
            poly2 = Polygon2D(points_contours2)
            d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d, return_other_point=True)
                    
            pt2 = self.points2d_to3d([new_pt2_2d], R2, r2, frame2)
            p2 = pt2[0]
        
        
        if not(cyl.contours2d[0].point_belongs(pt1_2d)) :
            #Find the closest one
            points_contours1 = cyl.contours2d[0].tessel_points
            
            poly1 = Polygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d, return_other_point=True)
                    
            pt1 = Point3D((r*math.cos(new_pt1_2d.vector[0]),
                           r*math.sin(new_pt1_2d.vector[0]),
                           new_pt1_2d.vector[1]))
            p1 = frame1.OldCoordinates(pt1)     
            
        return p1, p2
    
    
    def minimum_distance_points_plane(self, planeface): #Planeface with contour2D
        #### ADD THE FACT THAT PLANEFACE.CONTOURS : [0] = contours totale, le reste = trous
        
        poly2d = planeface.polygon2D
        pfpoints = poly2d.points
        xmin, ymin = min([pt[0] for pt in pfpoints]), min([pt[1] for pt in pfpoints])
        xmax, ymax = max([pt[0] for pt in pfpoints]), max([pt[1] for pt in pfpoints])
        origin, vx, vy = planeface.plane.origin, planeface.plane.vectors[0], planeface.plane.vectors[1] 
        pf1_2d, pf2_2d = Point2D((xmin, ymin)), Point2D((xmin, ymax))
        pf3_2d, pf4_2d = Point2D((xmax, ymin)), Point2D((xmax, ymax))
        pf1, pf2 = pf1_2d.To3D(origin, vx, vy), pf2_2d.To3D(origin, vx, vy) 
        pf3, _ = pf3_2d.To3D(origin, vx, vy), pf4_2d.To3D(origin, vx, vy)
        
        u, v = (pf3-pf1), (pf2-pf1)
        u.Normalize()
        v.Normalize()
        
        R1, r1 = self.rcenter, self.rcircle
        min_phi1, min_theta1, max_phi1, max_theta1 = self.minimum_maximum_tore(self.contours2d[0])
        
        n1 = self.normal
        u1 = self.toroidalsurface3d.frame.u
        v1 = self.toroidalsurface3d.frame.v
        frame1 = Frame3D(self.center, u1, v1, n1)
        # start1 = self.points2d_to3d([[min_theta1, min_phi1]], R1, r1, frame1)
        
        w = self.center - pf1
        
        n1n1, n1u1, n1v1, n1u, n1v = n1.Dot(n1), n1.Dot(u1), n1.Dot(v1), n1.Dot(u), n1.Dot(v)
        u1u1, u1v1, u1u, u1v = u1.Dot(u1), u1.Dot(v1), u1.Dot(u), u1.Dot(v)
        v1v1, v1u, v1v = v1.Dot(v1), v1.Dot(u), v1.Dot(v)
        uu, uv, vv = u.Dot(u), u.Dot(v), v.Dot(v)
        
        w2, wn1, wu1, wv1, wu, wv = w.Dot(w), w.Dot(n1), w.Dot(u1), w.Dot(v1), w.Dot(u), w.Dot(v)
        
        # x = (x, y, phi1, theta1)
        def distance_squared(x):
            return( uu*(x[0]**2) + vv*(x[1]**2) + w2 
                   + u1u1*(((R1+r1*math.cos(x[2]))*math.cos(x[3]))**2)
                   + v1v1*(((R1+r1*math.cos(x[2]))*math.sin(x[3]))**2)
                   + n1n1*((math.sin(x[2]))**2)*(r1**2)
                   + 2*x[0]*x[1]*uv - 2*x[0]*wu
                   - 2*x[0]*(R1+r1*math.cos(x[2]))*math.cos(x[3])*u1u
                   - 2*x[0]*(R1+r1*math.cos(x[2]))*math.sin(x[3])*v1u
                   - 2*x[0]*math.sin(x[2])*r1*n1u - 2*x[1]*wv
                   - 2*x[1]*(R1+r1*math.cos(x[2]))*math.cos(x[3])*u1v
                   - 2*x[1]*(R1+r1*math.cos(x[2]))*math.sin(x[3])*v1v
                   - 2*x[1]*math.sin(x[2])*r1*n1v
                   + 2*(R1+r1*math.cos(x[2]))*math.cos(x[3])*wu1
                   + 2*(R1+r1*math.cos(x[2]))*math.sin(x[3])*wv1
                   + 2*math.sin(x[2])*r1*wn1
                   + 2*u1v1*math.cos(x[3])*math.sin(x[3])*((R1+r1*math.cos(x[2]))**2)
                   + 2*(R1+r1*math.cos(x[2]))*math.cos(x[3])*r1*math.sin(x[2])*n1u1
                   + 2*(R1+r1*math.cos(x[2]))*math.sin(x[3])*r1*math.sin(x[2])*n1v1 )
                
                
        
        x01 = npy.array([(xmax-xmin)/2, (ymax-ymin)/2,
                         (min_phi1+max_phi1)/2, (min_theta1+max_theta1)/2 ])

        minimax = [(0, 0, min_phi1, min_theta1), (xmax-xmin, ymax-ymin, max_phi1, max_theta1)]
        
        res1 = scp.optimize.least_squares(distance_squared, x01, bounds=minimax)   

        # frame1 = Frame3D(self.center, u1, v1, n1)
        pt1 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R1, r1, frame1)
        p1 = pt1[0]
        p2 = pf1 + res1.x[2]*u + res1.x[3]*v
        
        pt1_2d = Point2D((res1.x[3], res1.x[2]))
        pt2_2d = p2.To2D(pf1,u, v)
        
        if not(self.contours2d[0].point_belongs(pt1_2d)) :
            #Find the closest one
            points_contours1 = self.contours2d[0].tessel_points
            
            poly1 = Polygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d, return_other_point=True)
                    
            pt1 = self.points2d_to3d([new_pt1_2d], R1, r1, frame1)
            p1 = pt1[0]
        
        if not(planeface.contours[0].point_belongs(pt2_2d)) :
            #Find the closest one
            d2, new_pt2_2d = planeface.polygon2D.PointBorderDistance(pt2_2d, return_other_point=True)
            
            p2 = new_pt2_2d.To3D(pf1, u, v)
        
        return p1, p2
    
    def minimum_distance(self, other_face, return_points=False) :
        if other_face.__class__ is ToroidalFace3D :
            p1, p2 = self.minimum_distance_points_tore(other_face)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
        
        if other_face.__class__ is CylindricalFace3D :
            p1, p2 = self.minimum_distance_points_cyl(other_face)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
        
        if other_face.__class__ is PlaneFace3D : 
            p1, p2 = self.minimum_distance_points_plane(other_face)
            if return_points : 
                return p1.point_distance(p2), p1, p2
            else :
                return p1.point_distance(p2)
        else :
            return NotImplementedError 
 
class ConicalFace3D (Face3D) :
    """
    :param contours2d: The Cone's contour2D 
    :type contours2d: Contour2D
    :param conicalsurface3d: Information about the Cone
    :type conicalsurface3d: ConicalSurface3D
    :param points: Contour2d's parameter Cone
    :type points: List of float
    
    """      
    
    def __init__(self, contours2d, conicalsurface3d, param, name=''):
        # self.rb = conicalsurface3d.radius
        # param : rb theta1 rt theta2 hmin hmax
        
        self.rb = param[0]
        self.rt = param[2]
        self.conicalsurface3d = conicalsurface3d 
        
        self.center = self.conicalsurface3d.frame.origin
        self.normal = self.conicalsurface3d.frame.w
        vec1 = self.conicalsurface3d.frame.u
        ptext1 = self.center + Point3D((self.rb*vec1).vector) + self.normal*param[4]
        center2 = self.center + param[5]*self.normal
        ptext2 =  center2 + Point3D((self.rt*vec1).vector) + param[5]*self.normal 
        if math.isclose(self.rb, 0, abs_tol = 1e-6) :
            c1 = LineSegment3D(self.center, center2)
        else :
            c1 = Arc3D(ptext1, ptext1.Rotation(self.center, self.normal, param[3]/2), ptext1.Rotation(self.center, self.normal, param[3]), self.normal) 
        if math.isclose(self.rt, 0, abs_tol = 1e-6) :
            c2 = LineSegment3D(center2, self.center)
        else :
            c2 = Arc3D(ptext2, ptext2.Rotation(center2, self.normal, param[1]/2), ptext2.Rotation(center2, self.normal, param[1]), self.normal)
            
        edges = [c1, c2]
        ctr = [Contour3D(edges)]
        
        Face3D.__init__(self, ctr)
        self.contours2d = contours2d 
        self.name = name 
        self.param = param
    
    @classmethod
    def from_contour3d(cls, contours3d, conicalsurface3d, name=''):
        """
        :param contours3d: The Cone's contour3D
        :type contours3d: Contour3D
        :param conicalsurface3d: Information about the Cone
        :type conicalsurface3d: ConicalSurface3D
        
        :Example:
            >>> contours3d is [Arc3D, LineSegment3D, Arc3D or LineSegment3D], the Cone's bones
        """
        frame = conicalsurface3d.frame
        new_pts = []
        for pt in contours3d[0].tessel_points:
            new_pts.append(frame.NewCoordinates(pt))
        size = len(contours3d[0].edges)
        
        if contours3d[0].edges[0].__class__ is Arc3D and contours3d[0].edges[1].__class__ is LineSegment3D and size<=4: #Portion of Cone
            
            hmax, hmin, Rmax, Rmin = ConicalFace3D.from_points3d_param(contours3d[0].tessel_points, frame)
            if contours3d[0].edges[2].__class__ is LineSegment3D :
                radius_max = max([Rmax, Rmin])
                arc = contours3d[0].edges[0]
                if arc.normal == -frame.w :
                    arc.setup_arc(arc.start, arc.interior, arc.end, -arc.normal)
                theta1, theta2 = posangle_arc(arc.start, arc.end, radius_max, frame)
                offset, angle = offset_angle(arc.is_trigo, theta1, theta2)
                pt1, pt2, pt3, pt4 = Point2D((offset, hmin)), Point2D((offset, hmax)), Point2D((offset+angle, hmax)), Point2D((offset+angle, hmin))
                    
            else :
                arc1, arc2 = contours3d[0].edges[0], contours3d[0].edges[2]
                r1, r2, n1, n2 = arc1.radius, arc2.radius, arc1.normal, arc2.normal
                if n1 == -frame.w :
                    arc1.setup_arc(arc1.start, arc1.interior, arc1.end, -arc1.normal)
                if n2 == -frame.w :
                    arc2.setup_arc(arc2.start, arc2.interior, arc2.end, -arc2.normal)
                start1, end1 = arc1.start, arc1.end
                start2, end2 = arc2.start, arc2.end
                if math.isclose(r1, Rmin, abs_tol=1e-6) :
                    theta1_1, theta1_2 = posangle_arc(start1, end1, Rmin, frame)
                    offset1, angle1 = offset_angle(arc1.is_trigo, theta1_1, theta1_2)
                    
                    theta2_1, theta2_2 = posangle_arc(start2, end2, Rmax, frame)
                    offset2, angle2 = offset_angle(arc2.is_trigo, theta2_1, theta2_2)
                    
                    pt1, pt2, pt3, pt4 = Point2D((offset1, hmin)), Point2D((offset2, hmax)), Point2D((offset2+angle2, hmax)), Point2D((offset1+angle1, hmin))
                elif math.isclose(r1, Rmax, abs_tol=1e-6) :
                    theta1_1, theta1_2 = posangle_arc(start2, end2, Rmin, frame)
                    offset1, angle1 = offset_angle(arc2.is_trigo, theta1_1, theta1_2)
                    
                    theta2_1, theta2_2 = posangle_arc(start1, end1, Rmax, frame)
                    offset2, angle2 = offset_angle(arc1.is_trigo, theta2_1, theta2_2)
                    
                    pt1, pt2, pt3, pt4 = Point2D((offset1, hmin)), Point2D((offset2, hmax)), Point2D((offset2+angle2, hmax)), Point2D((offset1+angle1, hmin))
                else :
                    print('r1, r2, Rmin, Rmax', r1, r2, Rmin, Rmax)
                    print('You should minimize the abs_tol')
                    raise NotImplementedError
                angle=max(angle1, angle2)
            #Creation of the window
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]
            param = [Rmin,angle,Rmax,angle,hmin,hmax]
            
            # fig, ax = plt.subplots()
            # [seg.MPLPlot(ax=ax) for seg in primitives]
            # [pt.MPLPlot(ax=ax, color='r') for pt in contours2d[0].tessel_points]
            
        elif contours3d[0].edges[1].__class__ is Arc3D and contours3d[0].edges[0].__class__ is LineSegment3D and size<=4: 
            hmax, hmin, Rmax, Rmin = ConicalFace3D.from_points3d_param(contours3d[0].tessel_points, frame)
            
            if size > 3 :
                if contours3d[0].edges[3].__class__ is Arc3D :
                    arc1, arc2 = contours3d[0].edges[1], contours3d[0].edges[3]
                    r1, r2, n1, n2 = arc1.radius, arc2.radius, arc1.normal, arc2.normal
                    if n1 == -frame.w :
                        arc1.setup_arc(arc1.start, arc1.interior, arc1.end, -arc1.normal)
                    if n2 == -frame.w :
                        arc2.setup_arc(arc2.start, arc2.interior, arc2.end, -arc2.normal)
                    start1, end1 = arc1.start, arc1.end
                    start2, end2 = arc2.start, arc2.end
                    if math.isclose(r1, Rmin, abs_tol=1e-6) :
                        theta1_1, theta1_2 = posangle_arc(start1, end1, Rmin, frame)
                        offset1, angle1 = offset_angle(arc1.is_trigo, theta1_1, theta1_2)
                        
                        theta2_1, theta2_2 = posangle_arc(start2, end2, Rmax, frame)
                        offset2, angle2 = offset_angle(arc2.is_trigo, theta2_1, theta2_2)
                        
                        pt1, pt2, pt3, pt4 = Point2D((offset1, hmin)), Point2D((offset2, hmax)), Point2D((offset2+angle2, hmax)), Point2D((offset1+angle1, hmin))
                    elif math.isclose(r1, Rmax, abs_tol=1e-6) :
                        theta1_1, theta1_2 = posangle_arc(start2, end2, Rmin, frame)
                        offset1, angle1 = offset_angle(arc2.is_trigo, theta1_1, theta1_2)
                            
                        theta2_1, theta2_2 = posangle_arc(start1, end1, Rmax, frame)
                        offset2, angle2 = offset_angle(arc1.is_trigo, theta2_1, theta2_2)
                        
                        pt1, pt2, pt3, pt4 = Point2D((offset1, hmin)), Point2D((offset2, hmax)), Point2D((offset2+angle2, hmax)), Point2D((offset1+angle1, hmin))
                    else :
                        print('r1, r2, Rmin, Rmax', r1, r2, Rmin, Rmax)
                        print('You should minimize the abs_tol')
                        raise NotImplementedError
                    angle=max(angle1, angle2)
                else :
                    radius_max = max([Rmax, Rmin])
                    arc = contours3d[0].edges[1]
                    if arc.normal == -frame.w :
                        arc.setup_arc(arc.start, arc.interior, arc.end, -arc.normal)
                    theta1, theta2 = posangle_arc(arc.start, arc.end, radius_max, frame)
                    offset, angle = offset_angle(arc.is_trigo, theta1, theta2)
                    pt1, pt2, pt3, pt4 = Point2D((offset, hmin)), Point2D((offset, hmax)), Point2D((offset+angle, hmax)), Point2D((offset+angle, hmin))
            else :
                radius_max = max([Rmax, Rmin])
                arc = contours3d[0].edges[1]
                if arc.normal == -frame.w :
                    arc.setup_arc(arc.start, arc.interior, arc.end, -arc.normal)
                theta1, theta2 = posangle_arc(arc.start, arc.end, radius_max, frame)
                offset, angle = offset_angle(arc.is_trigo, theta1, theta2)
                pt1, pt2, pt3, pt4 = Point2D((offset, hmin)), Point2D((offset, hmax)), Point2D((offset+angle, hmax)), Point2D((offset+angle, hmin))
            
            #Creation of the window
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]
            param = [Rmin,angle,Rmax,angle,hmin,hmax]
            
            # fig, ax = plt.subplots()
            # [seg.MPLPlot(ax=ax) for seg in primitives]
            # [pt.MPLPlot(ax=ax, color='g') for pt in contours2d[0].tessel_points]
            
        elif contours3d[0].edges[0].__class__ is LineSegment3D and contours3d[0].edges[1].__class__ is LineSegment3D and size==3 :
            hmax, hmin, Rmax, Rmin = ConicalFace3D.from_points3d_param(contours3d[0].tessel_points, frame)
            
            if contours3d[0].edges[2].__class__ is Arc3D:
                radius_max = max([Rmax, Rmin])
                arc = contours3d[0].edges[2]
                if arc.normal == -frame.w :
                    arc.setup_arc(arc.start, arc.interior, arc.end, -arc.normal)
                theta1, theta2 = posangle_arc(arc.start, arc.end, radius_max, frame)
                offset, angle = offset_angle(arc.is_trigo, theta1, theta2)
                #Creation of the window
                pt1, pt2, pt3, pt4 = Point2D((offset, hmin)), Point2D((offset, hmax)), Point2D((offset+angle, hmax)), Point2D((offset+angle, hmin))
                seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
                primitives = [seg1, seg2, seg3, seg4]
                contours2d =  [Contour2D(primitives)]
                param = [Rmin,angle,Rmax,angle,hmin,hmax]
                
            else :
                contours2d = ConicalFace3D.contours3d_to2d(contours3d, conicalsurface3d)
                theta = max(pt[0] for pt in contours2d[0].tessel_points) - min(pt[0] for pt in contours2d[0].tessel_points)
                param = [Rmin,angle,Rmax,angle,hmin,hmax]

        elif contours3d[0].edges[0].__class__ is Arc3D and contours3d[0].edges[2].__class__ is Arc3D and size<=4: 
            r1 = contours3d[0].edges[0].radius
            t1 = contours3d[0].edges[0].angle
            r2 = contours3d[0].edges[2].radius
            t2 = contours3d[0].edges[2].angle
            
            hmax, hmin, Rmax, Rmin = ConicalFace3D.from_points3d_param(contours3d[0].tessel_points, frame)
            if math.isclose(Rmin, r1, abs_tol=1e-4) :
                theta1, theta2 = t1, t2
            else : 
                theta2, theta1 = t1, t2
            
            offset = (theta1 - theta2)/2
            
            #Creation of the window
            pt1, pt2, pt3, pt4 = Point2D((offset, hmin)), Point2D((offset, hmax)), Point2D((offset+theta2, hmax)), Point2D((theta1, hmin))
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]
            param = [Rmin,theta1,Rmax,theta2,hmin,hmax]
            
        else:
            # fig = plt.figure()
            # ax = fig.add_subplot(111, projection='3d')
            # [edge.MPLPlot(ax=ax) for edge in contours3d[0].edges]
            contours2d = ConicalFace3D.contours3d_to2d(contours3d, conicalsurface3d)
            hmax, hmin, Rmax, Rmin = ConicalFace3D.from_points3d_param(contours3d[0].tessel_points, frame)
            theta = max(pt[0] for pt in contours2d[0].tessel_points) - min(pt[0] for pt in contours2d[0].tessel_points)
            param = [Rmin,theta,Rmax,theta,hmin,hmax]
            # print('contours3d edges', contours3d[0].edges)
            # raise NotImplementedError
         
        # fig, ax = plt.subplots()
        # [seg.MPLPlot(ax=ax) for seg in contours2d[0].primitives]
        # [pt.MPLPlot(ax=ax, color='g') for pt in contours2d[0].tessel_points]   
         
        return cls(contours2d, conicalsurface3d, param, name=name)
    
    def points2d_to3d(self, all_contours_points, rb, rt, hmin, hmax, frame3d) :
        Points3D = []
        for listpt in  all_contours_points: 
            for pt in listpt :
                ratio = (pt[1]-hmin)/(hmax-hmin)
                radius = (1-ratio)*rb + ratio*rt
                Points3D.append(Point3D((radius*math.cos(pt[0]),radius*math.sin(pt[0]),pt[1])))            
        Points_3D = [frame3d.OldCoordinates(point) for point in Points3D]
        return Points_3D
    
    def points3d_to2d(points3d, rb, rt, hmin, hmax) :
        points_2D = []
        for pt in points3d :
            x, y, h = pt[0], pt[1], pt[2]
            ratio = (h-hmin)/(hmax-hmin)
            radius = (1-ratio)*rb + ratio*rt
            u1, u2 = x/radius, y/radius
            theta = sin_cos_angle(u1, u2)
            
            points_2D.append(Point2D([theta, h]))  
            
        i0, i2pi, ih, iangle = 0, 0, 0, 0
        for enum,point in enumerate(points_2D) :
            if enum == 0 :
                h = point[1]
                angle = point[0]
            if math.isclose(point[0], 0, abs_tol=1e-6) :
                i0 += 1
            elif math.isclose(point[0], 2*math.pi, abs_tol=1e-6) :
                i2pi += 1
            elif math.isclose(point[0], angle, abs_tol=5e-2) :
                iangle += 1
            elif math.isclose(point[1], h, abs_tol=1e-6) :
                ih += 1
        
        if i2pi/len(points_2D)>0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((2*math.pi, point.vector[1])))
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif i0/len(points_2D)>0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((2*math.pi, point.vector[1]))) 
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif iangle/len(points_2D)>0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((angle, point.vector[1]))) 
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif ih/len(points_2D)>0.5 :
            new_points2d = [pt.copy() for pt in points_2D]
            new_points2d.sort(key=lambda pt: pt[0])
            if i2pi == 2 :
                pt1, pt2 = new_points2d[0], new_points2d[-1]
                pt1.vector[0] = 0-0.0000001
                pt2.vector[0] = 2*math.pi+0.0000001
                points_2D = [pt1, pt2]
            else :
                points_2D = [new_points2d[0], new_points2d[-1]]
        
        plus_7pi4, moins_7pi4 = [], []
        for enum, pt in enumerate(points_2D) :
            if pt.vector[0]>7*math.pi/4 :
                plus_7pi4.append(enum)
            elif pt.vector[0]<math.pi/4 :
                moins_7pi4.append(enum)
            
        if len(moins_7pi4) > 3*len(plus_7pi4) :
            for pos in plus_7pi4 :
                new_pt = points_2D[pos].copy() - Point2D((2*math.pi, 0))
                points_2D[pos] = new_pt
        # points_2D = check_singularity(points_2D)     
        
        return points_2D
    
    def from_points3d_param(globalpoints, conical_frame):
        #globalpoints should be points in global frame
        new_pts = []
        for pt in globalpoints:
            new_pts.append(conical_frame.NewCoordinates(pt))
        hmax, hmin = max(pt[2] for pt in new_pts), min(pt[2] for pt in new_pts)  
        pts_max, pts_min = [], []
        for pt in new_pts :
            if math.isclose(pt[2], hmin, abs_tol = 1e-6) : 
                pts_min.append(pt)
            elif math.isclose(pt[2], hmax, abs_tol = 1e-6) :
                pts_max.append(pt)
            else :
                continue
        list_rmax, list_rmin = [], []
        center = conical_frame.origin
        new_center = conical_frame.NewCoordinates(center)
        cmax = Point3D((new_center.vector[0], new_center.vector[1], hmax))
        cmin = Point3D((new_center.vector[0], new_center.vector[1], hmin))
        for pt in pts_max : 
            list_rmax.append((pt-cmax).Norm())
        Rmax = sum(list_rmax)/len(list_rmax)
        for pt in pts_min :
            list_rmin.append((pt-cmin).Norm())
        Rmin = sum(list_rmin)/len(list_rmin)
        
        return hmax, hmin, Rmax, Rmin
     
    def contours3d_to2d(contours3d, conicalsurface3d) :
        frame = conicalsurface3d.frame
        n = frame.w
        # center = frame.origin
        
        # print('\n contours3d[0].edges', contours3d[0].edges)
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # [edge.MPLPlot(ax=ax) for edge in contours3d[0].edges]
            
        hmax, hmin, Rmax, Rmin = ConicalFace3D.from_points3d_param(contours3d[0].tessel_points, frame)
        
        if math.isclose(Rmin, 0, abs_tol=1e-6) :
            Rmin = 1e-6
        
        primitives, start_end, all_points = [], [], []
        for edge in contours3d[0].edges :
            new_points = [frame.NewCoordinates(pt) for pt in edge.points]
            if edge.__class__ is Arc3D :
                if edge.normal == n or edge.normal == -n:
                    start2d, end2d = ConicalFace3D.points3d_to2d(new_points, edge.radius, edge.radius, hmin, hmax)
                    angle2d = abs(end2d[0]-start2d[0])
                    if math.isclose(edge.angle, 2*math.pi, abs_tol=1e-6):
                        if start2d == end2d :
                            if math.isclose(start2d.vector[0], 2*math.pi, abs_tol = 1e-6) :
                                end2d = end2d - Point2D((2*math.pi, 0))
                            else :
                                end2d = end2d + Point2D((2*math.pi, 0))
                    elif not(math.isclose(edge.angle, angle2d, abs_tol=1e-2)) :
                        # if math.isclose(angle2d, 2*math.pi, abs_tol=1e-2) :
                        if start2d[0] < end2d[0] :
                            end2d = start2d + Point2D((edge.angle,0))
                        else :
                            end2d = start2d - Point2D((edge.angle,0))
                    ls_toadd = LineSegment2D(start2d, end2d)
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        primitives.append([ls_toadd])
                        all_points.extend(ls_toadd.points)
                        start_end.append(ls_toadd.points)
                    
                else :
                    points2d = ConicalFace3D.points3d_to2d(new_points, edge.radius, edge.radius, hmin, hmax)
                    lines = []
                    for k in range(0, len(points2d)-1) :
                        lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                    points, prim_list = [], []
                    for ls_toadd in lines :        
                        same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                        if same is False :
                            prim_list.append(ls_toadd)
                            points.extend(ls_toadd.points)
                    if len(points) > 0 : 
                        all_points.extend(points)
                        primitives.append(prim_list)
                        start_end.append([points[0], points[-1]])
                    
            elif edge.__class__ is LineSegment3D :
                start2d, end2d = ConicalFace3D.points3d_to2d(new_points, Rmin, Rmax, hmin, hmax)
                ls_toadd = LineSegment2D(start2d, end2d)
                same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                if same is False :
                    primitives.append([ls_toadd])
                    all_points.extend(ls_toadd.points)
                    start_end.append(ls_toadd.points)
            
            else :
                points2d = ConicalFace3D.points3d_to2d(new_points, Rmin, Rmax, hmin, hmax)
                lines = []
                for k in range(0, len(points2d)-1) :
                    lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                points, prim_list = [], []
                for ls_toadd in lines :        
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        prim_list.append(ls_toadd)
                        points.extend(ls_toadd.points)
                if len(points) > 0 : 
                    all_points.extend(points)
                    primitives.append(prim_list)
                    start_end.append([points[0], points[-1]])
        
        # fig, ax = plt.subplots() 
        # [pt.MPLPlot(ax=ax) for pt in all_points]
        
        points_se, primitives_se = [], []
        for double in start_end : 
            primitives_se.append(LineSegment2D(double[0], double[1]))
            points_se.extend(double)
        poly_se = Polygon2D(points_se)
            
        xmax, xmin = max(pt[0] for pt in points_se), min(pt[0] for pt in points_se)
        ymax, ymin = max(pt[1] for pt in points_se), min(pt[1] for pt in points_se) 
        pt1, pt2, pt3, pt4 = Point2D((xmin, ymin)), Point2D((xmin, ymax)), Point2D((xmax, ymin)), Point2D((xmax, ymax))
        diag1, diag2 = LineSegment2D(pt1, pt4), LineSegment2D(pt2, pt3)
        diag1_cut, diag2_cut = [], []
        diag1_pointcut, diag2_pointcut = [], []
        for enum,l in enumerate(primitives_se) :
            cut1 = diag1.line_intersection(l)
            cut2 = diag2.line_intersection(l)
            if cut1 is not None : 
                diag1_cut.append(enum)
                diag1_pointcut.append(cut1)
            if cut2 is not None : 
                diag2_cut.append(enum)
                diag2_pointcut.append(cut2)
        
        points_common = []
        for enum1,pos1 in enumerate(diag1_cut) :
            for enum2,pos2 in enumerate(diag2_cut) :
                if pos1 == pos2 :
                    points_common.append(primitives_se[pos1].points)
        # fig, ax = plt.subplots() 
        # [l.MPLPlot(ax=ax) for l in primitives_se]
        # [pt.MPLPlot(ax=ax, color='g') for pt in points_se]
        # pt1.MPLPlot(ax=ax, color='r')
        # pt2.MPLPlot(ax=ax, color='r')            
        # pt3.MPLPlot(ax=ax, color='r')
        # pt4.MPLPlot(ax=ax, color='r')   
        # diag1.MPLPlot(ax=ax)
        # diag2.MPLPlot(ax=ax)
        # for couple in points_common :
        #     [pt.MPLPlot(ax=ax, color='b') for pt in couple]        
        
        if len(points_common) >= 1 :
            solve = False
            for couple in points_common :
                if solve :
                    break
                check1, check2 = poly_se.PointBelongs(couple[0]), poly_se.PointBelongs(couple[1])
                start, end = couple[0].vector[0], couple[1].vector[0]
                if math.isclose(start, end, abs_tol = 5e-2) :
                    intersect = min(start, end)
                    if math.isclose(intersect, math.pi, abs_tol = 5e-2) :
                        points_sing = check_singularity(all_points)
                        pt0, pt2pi = 0, 0
                        for pt in points_sing :
                            if math.isclose(pt.vector[0], 0, abs_tol=1e-2):
                                pt0 += 1
                            elif math.isclose(pt.vector[0], 2*math.pi, abs_tol=1e-2):
                                pt2pi += 1
                        points_sing.sort(key=lambda pt: pt[1])
                        points_sing.sort(key=lambda pt: pt[0])
                        if pt2pi !=0 and pt0 == 0 :
                            points = [pt.copy() for pt in points_sing[::-1]]
                            points_sing = points
                        points_range = Face3D.range_closest(points_sing)
                        all_points = delete_double_point(points_range)
                        # intersect = 0
                        break
                    elif math.isclose(intersect, 0, abs_tol = 1e-6) or math.isclose(intersect, 2*math.pi, abs_tol = 1e-6) or (not check1 or not check2):
                        all_points = check_singularity(all_points)
                        points_cleaned = delete_double_point(all_points)
                        all_points = [pt.copy() for pt in points_cleaned]
                        all_points.sort(key=lambda pt: pt[0])
                        d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                        if d2 < d1 :
                            last = all_points[-1].copy()
                            all_points[-1] = all_points[-2].copy()
                            all_points[-2] = last
                        break
                    else :
                        points = []
                        for list_prim in primitives :
                            for k,prim in enumerate(list_prim) :
                                new_list_points = []
                                change = 0
                                for pt in prim.points :
                                    if pt[0] < intersect :
                                        change += 1 
                                        if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                            new_list_points.append(Point2D((intersect + 2*math.pi, pt[1])))
                                        else :
                                            new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    elif math.isclose(pt[0], intersect, abs_tol=1e-4) :
                                        change += 1
                                        new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    else :
                                        new_list_points.append(pt)
                                if change > 0 :
                                    points.extend(new_list_points)
                                    list_prim[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                                else :
                                    points.extend(prim.points)
                                    continue
                        points_cleaned = delete_double_point(points)
                        all_points = Face3D.range_trigo(points_cleaned)
                    solve = True
                else :
                    points_cleaned = delete_double_point(all_points)
                    all_points = [pt.copy() for pt in points_cleaned]
                    all_points.sort(key=lambda pt: pt[0])
                    d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                    if d2 < d1 :
                        last = all_points[-1].copy()
                        all_points[-1] = all_points[-2].copy()
                        all_points[-2] = last
        else :
            points_cleaned = delete_double_point(all_points)
            all_points = [pt.copy() for pt in points_cleaned]
            all_points.sort(key=lambda pt: pt[0])
            d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
            if d2 < d1 :
                last = all_points[-1].copy()
                all_points[-1] = all_points[-2].copy()
                all_points[-2] = last
                    
        primitives = Face3D.create_primitives(all_points)
        
        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax, color='g') for pt in all_points]
        # all_points[0].MPLPlot(ax=ax, color='m')
        # all_points[1].MPLPlot(ax=ax, color='r')
        # all_points[-1].MPLPlot(ax=ax)
        # all_points[-2].MPLPlot(ax=ax, color='b')
        # [p.MPLPlot(ax=ax) for p in primitives]
        
        l_vert = LineSegment2D((pt2+pt4)/2, (pt1+pt3)/2)
        solve = False
        for prim in primitives :
            if solve :
                break
            intersect = prim.line_intersection(l_vert)
            if intersect is not None : 
                x_intersect = intersect.vector[0]
                y_intersect = intersect.vector[1]
                value1, value2 = ymax-0.2*(ymax-ymin), ymin+0.2*(ymax-ymin)
                if y_intersect < max(value1, value2) and y_intersect > min(value1, value2) :
                    points = []
                    for k, prim in enumerate(primitives) :
                        new_list_points, change = [], 0
                        for pt in prim.points :
                            if pt[0] < x_intersect :
                                change += 1
                                if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                    new_list_points.append(Point2D((x_intersect + 2*math.pi, pt[1])))
                                else :
                                    new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                            else :
                                new_list_points.append(pt)
                        if change > 0:
                            points.extend(new_list_points)
                            # primitives[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                        else :
                            points.extend(prim.points)
                            continue
                    solve = True
                    points_cleaned = delete_double_point(points)
                    all_points = Face3D.range_trigo(points_cleaned)
                    primitives = Face3D.create_primitives(all_points)
                    
        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax) for pt in all_points]
        # [p.MPLPlot(ax=ax) for p in primitives]
        # intersect.MPLPlot(ax=ax, color='r')
        # l_vert.MPLPlot(ax=ax)
        contour2d = [Contour2D(primitives)]
        return contour2d
    
    def create_triangle(self, all_contours_points, part) :
        Triangles, ts = [], []
        pts, h_list = [], []
        for listpt in all_contours_points :
            for pt in listpt :
                pts.append(pt)
                h_list.append(pt[1])
        if part == 'bot' :        
            h_concerned = min(h_list)
        else :
            h_concerned = max(h_list)
        peak_list, other = [], []
        for pt in pts : 
            if pt[1]==h_concerned :
                peak_list.append(pt)
            else : 
                other.append(pt)
        points = [peak_list[0]] + other
        
        for i in range(1, len(points)):
            if i == len(points)-1 :
                vertices = [points[i].vector, points[0].vector, points[1].vector]
                segments = [[0,1],[1,2],[2,0]]
                listindice = [i,0,1]
            else :
                vertices = [points[i].vector, points[0].vector, points[i+1].vector]
                segments = [[0,1],[1,2],[2,0]]
                listindice = [i,0,i+1]
            tri = {'vertices': vertices, 'segments': segments}
            t = triangle.triangulate(tri, 'p')
            if 'triangles' in t:
                triangles = t['triangles'].tolist()
                triangles[0] = listindice
                Triangles.append(triangles)        
            else:
                Triangles.append(None)
            ts.append(t)
            
        return points, Triangles
    
    def triangulation(self, resolution=30):
        hmin = self.param[4]
        hmax = self.param[5]
        
        frame3d = self.conicalsurface3d.frame
        
        angle_theta1 = self.param[1]
        angle_theta2 = self.param[3]
        pas_theta = 2*math.pi/30 #Step of 12 degrees
         
        resolution_theta = abs(int(max(angle_theta1,angle_theta2)/pas_theta))
        if resolution_theta < 5 :
            resolution_theta = 5

        all_contours_points = self.cut_contours(self.contours2d, resolution_theta)
        # print('all_contours_points', all_contours_points)

        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax, color='r') for pt in self.contours2d[0].tessel_points]
        # [prim.MPLPlot(ax=ax) for prim in self.contours2d[0].primitives]

        # fig, ax = plt.subplots()
        # for list_pt in all_contours_points :
        #     [pt.MPLPlot(ax=ax) for pt in list_pt]

        if self.rb == 0 :
            # print('>>>>>>>1')
            points, Triangles = self.create_triangle(all_contours_points, part='bot')
            Points_3D = self.points2d_to3d([points], self.rb, self.rt, hmin, hmax, frame3d)
        elif self.rt == 0 :
            # print('>>>>>>>2')
            points, Triangles = self.create_triangle(all_contours_points, part='top')
            Points_3D = self.points2d_to3d([points], self.rb, self.rt, hmin, hmax, frame3d)
        else :    
            # print('>>>>>>>3')
            Triangles, ts = [], []
            for k, listpt in enumerate(all_contours_points) :
                vertices=[]
                segments=[]
                for i, pt in enumerate(listpt):
                    vertices.append(pt.vector)
                    if i == len(listpt)-1 :
                        segments.append([i, 0])
                    else :
                        segments.append([i,i+1])
                segments[-1]=(len(listpt)-1,0) 
                tri = {'vertices': vertices, 'segments': segments}
                t = triangle.triangulate(tri, 'p')
                ts.append(t)
            
            seuil, k = 0, 0
            for t in ts :
                if 'triangles' in t:
                    triangles = t['triangles'].tolist()
                    for n,tri in enumerate(triangles):
                        for i in range (0,3):
                            tri[i]=tri[i]+seuil
                    seuil += len(all_contours_points[k])
                    k += 1
                    Triangles.append(triangles)        
                else:
                    Triangles.append(None)
        
            Points_3D = self.points2d_to3d(all_contours_points, self.rb, self.rt, hmin, hmax, frame3d)
        
        pt3d, tangle = delete_double_pos(Points_3D, Triangles)
        
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # [pt.MPLPlot(ax=ax, color='r') for pt in pt3d]  
        
        # print('>>>>len(Points_3D)', len(Points_3D))
        # print('>>>>>len(pt3d)', len(pt3d)) 
        # print('>>>>>>len(tangle)', len(tangle))
        
        return pt3d, tangle 
    
    def frame_mapping(self, frame, side, copy=True) :
        if copy:
            new_conicalsurface3d = ConicalSurface3D.frame_mapping(frame, side, copy)
            return ConicalFace3D(self.contours2d, new_conicalsurface3d, points=self.points, name=self.name)
        else:
            self.conicalsurface3d.frame_mapping(frame, side, copy=False)
  
class SphericalFace3D (Face3D) :
    """
    :param contours2d: The Sphere's contour2D 
    :type contours2d: Contour2D
    :param sphericalsurface3d: Information about the Sphere
    :type sphericalsurface3d: SphericalSurface3D
    :param points: Angle's Sphere
    :type points: List of float
    
    :Example: 
        >>> contours2d is rectangular and will create a classic tore with x:2*pi, y:2*pi
        
        >>> points = [pi, 2*pi] for an half sphere
    """      
    
    def __init__(self, contours2d, sphericalsurface3d, points=None, name=''):
        self.radius = sphericalsurface3d.radius
        self.sphericalsurface3d = sphericalsurface3d 
        
        self.center = self.sphericalsurface3d.frame.origin
        self.normal = self.sphericalsurface3d.frame.w
        vec1 = self.sphericalsurface3d.frame.u
        ptext1 = self.center + Point3D((self.radius*vec1).vector) 
        ptext2 = self.center + Point3D((self.radius*vec1).vector) 
        c1 = Arc3D(ptext1, ptext1.Rotation(self.center, self.normal, points[0]/2), ptext1.Rotation(self.center, self.normal, points[0]), self.normal) 
        c2 = Arc3D(ptext2, ptext2.Rotation(self.center, self.normal, points[1]/2), ptext2.Rotation(self.center, self.normal, points[1]), self.normal)
            
        edges = [c1, c2]
        ctr = [Contour3D(edges)]
        
        Face3D.__init__(self, ctr)
        self.contours2d = contours2d 
        
        self.points = points 
        self.name = name 
    
        # pts3d, t = self.triangulation()
        # self.start = pts3d[0]
    
    @classmethod
    def from_contour3d(cls, contours3d, sphericalsurface3d, name=''):
        """
        :param contours3d: The Sphere's contour3D
        :type contours3d: Contour3D
        :param sphericalsurface3d: Information about the Sphere
        :type sphericalsurface3d: SphericalSurface3D
        
        """
        frame = sphericalsurface3d.frame
        center = frame.origin
        normal = frame.w
        # r = sphericalsurface3d.radius
        
        if contours3d[0].__class__ is Point3D : #If it is a complete sphere
            angle = 2*math.pi
            pt1, pt2, pt3, pt4 = Point2D((0, 0)), Point2D((0, angle)), Point2D((angle, angle)), Point2D((angle, 0))
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]
            points = [angle, angle]
         
        elif contours3d[0].edges[0].__class__ is Arc3D and contours3d[0].edges[1].__class__ is Arc3D : #Portion of Sphere
            ## we supposed a contours with 4 edges maximum here
            # fig = plt.figure()
            # ax = fig.add_subplot(111, projection='3d')
            # [edge.MPLPlot(ax=ax) for edge in contours3d[0].edges]
            # print('contours3d[0].edges', contours3d[0].edges)
            # print('normal', normal)
            arc_base, arc_link = [], []
            range_list = []
            for edge in contours3d[0].edges :
                # print('edge.normal', edge.normal)
                if edge.normal == normal :
                    arc_base.append(edge)
                    range_list.append(1)
                else :
                    arc_link.append(edge)
                    range_list.append(2)
            if len(arc_base) == 0 :
                raise NotImplementedError
                 
            # radius = sphericalsurface3d.radius 
            theta = arc_base[0].angle
            phi = arc_link[0].angle
            if range_list[-1] == range_list[-2] :
                offset_phi = -math.pi/2
            else :
                pos = len(arc_base)-1
                c1 = arc_base[pos].center
                vec1, vec2 = c1 - center, arc_base[pos].start - center
                offset_phi = -math.pi/2 + vectors3d_angle(vec1, vec2)
            offset_theta = 0
            
            # Creation of the window
            pt1, pt2, pt3, pt4 = Point2D((offset_theta, offset_phi)), Point2D((offset_theta, offset_phi + phi)), Point2D((offset_theta+theta, offset_phi + phi)), Point2D((offset_theta+theta, offset_phi))
            seg1, seg2, seg3, seg4 = LineSegment2D(pt1, pt2), LineSegment2D(pt2, pt3), LineSegment2D(pt3, pt4), LineSegment2D(pt4, pt1) 
            primitives = [seg1, seg2, seg3, seg4]
            contours2d =  [Contour2D(primitives)]
            points = [theta, phi]
            
            # fig, ax = plt.subplots()
            # [prim.MPLPlot(ax=ax) for prim in primitives]
                
        else: 
            contours2d = SphericalFace3D.contours3d_to2d(contours3d, sphericalsurface3d)
            theta = max(pt[0] for pt in contours2d[0].tessel_points) - min(pt[0] for pt in contours2d[0].tessel_points)
            phi = max(pt[1] for pt in contours2d[0].tessel_points) - min(pt[1] for pt in contours2d[0].tessel_points)
            points = [theta, phi]
            # print('contours3d edges', contours3d[0].edges)
            # raise NotImplementedError
            
        return cls(contours2d, sphericalsurface3d, points, name=name)
 
    def points2d_to3d(self, points2d, r, frame3d) :
        # source mathcurve.com/surfaces/sphere
        # -pi<theta<pi, -pi/2<phi<pi/2
        points3D = []
        for pt in points2d :
            theta, phi = pt[0], pt[1]
            x = r*math.cos(phi)*math.cos(theta)
            y = r*math.cos(phi)*math.sin(theta)
            z = r*math.sin(phi)
            points3D.append(Point3D([x,y,z]))
        Points_3D = [frame3d.OldCoordinates(point) for point in points3D]
        return Points_3D
    
    def points3d_to2d(points3d, r) :
        points_2D = []
        for pt in points3d :
            x, y, z = pt[0], pt[1], pt[2]
            if z<-r :
                z = -r
            elif z>r :
                z = r
                
            zr = z/r
            phi = math.asin(zr)
            
            u = math.sqrt((r**2)-(z**2))
            if u == 0 :
                u1, u2 = x, y
            else :
                u1, u2 = round(x/u, 5), round(y/u, 5)
            #cos(theta)=u1, sin(theta)=u2 :
            theta = sin_cos_angle(u1, u2)
                    
            points_2D.append(Point2D([theta, phi]))
            
        i0, i2pi, iangle, ih = 0, 0, 0, 0
        for enum,point in enumerate(points_2D) :
            if enum == 0 :
                h = point[1]
                angle = point[0]                
            if math.isclose(point[0], 0, abs_tol=1e-6) :
                i0 += 1
            elif math.isclose(point[0], 2*math.pi, abs_tol=1e-6) :
                i2pi += 1
            elif math.isclose(point[0], angle, abs_tol=5e-2) :
                iangle += 1
            elif math.isclose(point[1], h, abs_tol=5e-2) :
                ih += 1
    
        if i2pi/len(points_2D)>=0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((2*math.pi, point.vector[1])))
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif i0/len(points_2D)>=0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((2*math.pi, point.vector[1]))) 
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif iangle/len(points_2D)>=0.5 :
            new_points2d = []
            for point in points_2D :
                new_points2d.append(Point2D((angle, point.vector[1]))) 
            new_points2d.sort(key=lambda pt: pt[1])
            points_2D = [new_points2d[0], new_points2d[-1]]
        elif ih/len(points_2D)>=0.5 :
            new_points2d = [pt.copy() for pt in points_2D]
            new_points2d.sort(key=lambda pt: pt[0])
            if i2pi == 2 or i0 == 2:
                pt1, pt2 = new_points2d[0], new_points2d[-1]
                pt1.vector[0] = 0-0.0000001
                pt2.vector[0] = 2*math.pi+0.0000001
                points_2D = [pt1, pt2]
            else :
                points_2D = [new_points2d[0], new_points2d[-1]]
                
        plus_7pi4, moins_7pi4 = [], []
        for enum, pt in enumerate(points_2D) :
            if pt.vector[0]>7*math.pi/4 :
                plus_7pi4.append(enum)
            elif pt.vector[0]<math.pi/4 :
                moins_7pi4.append(enum)
        
        if len(moins_7pi4) > 3*len(plus_7pi4) :
            for pos in plus_7pi4 :
                new_pt = points_2D[pos].copy() - Point2D((2*math.pi, 0))
                points_2D[pos] = new_pt
                
        return points_2D
    
    def contours3d_to2d(contours3d, sphericalsurface3d) :
        frame = sphericalsurface3d.frame
        n = frame.w
        # center = frame.origin
        radius = sphericalsurface3d.radius
        
        # print('contours3d[0].edges', contours3d[0].edges)
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # [edge.MPLPlot(ax=ax) for edge in contours3d[0].edges]
        
        primitives, start_end, all_points = [], [], []
        for edge in contours3d[0].edges :
            new_points = [frame.NewCoordinates(pt) for pt in edge.points]
            if edge.__class__ is Arc3D :
                if edge.normal == n or edge.normal == -n:
                    start2d, end2d = SphericalFace3D.points3d_to2d(new_points, radius)
                    angle2d = abs(end2d[0]-start2d[0])
                    if math.isclose(edge.angle, 2*math.pi, abs_tol=1e-6):
                        if start2d == end2d :
                            if math.isclose(start2d.vector[0], 2*math.pi, abs_tol = 1e-6) :
                                end2d = end2d - Point2D((2*math.pi, 0))
                            else :
                                end2d = end2d + Point2D((2*math.pi, 0))
                    elif not(math.isclose(edge.angle, angle2d, abs_tol=1e-2)) :
                        # if math.isclose(angle2d, 2*math.pi, abs_tol=1e-2) :
                        if start2d[0] < end2d[0] :
                            end2d = start2d + Point2D((edge.angle,0))
                        else :
                            end2d = start2d - Point2D((edge.angle,0))
                    
                    ls_toadd = LineSegment2D(start2d, end2d)
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        primitives.append([ls_toadd])
                        all_points.extend(ls_toadd.points)
                        start_end.append(ls_toadd.points)
                    
                else :
                    points2d = SphericalFace3D.points3d_to2d(new_points, radius)
                    lines = []
                    for k in range(0, len(points2d)-1) :
                        lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                    points, prim_list = [], []
                    for ls_toadd in lines :        
                        same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                        if same is False :
                            prim_list.append(ls_toadd)
                            points.extend(ls_toadd.points)
                    if len(points) > 0 : 
                        all_points.extend(points)
                        primitives.append(prim_list)
                        start_end.append([points[0], points[-1]])
                    
            elif edge.__class__ is LineSegment3D :
                start2d, end2d = SphericalFace3D.points3d_to2d(new_points, radius)
                ls_toadd = LineSegment2D(start2d, end2d)
                same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                if same is False :
                    primitives.append([ls_toadd])
                    all_points.extend(ls_toadd.points)
                    start_end.append(ls_toadd.points)
            
            else :
                points2d = SphericalFace3D.points3d_to2d(new_points, radius)
                lines = []
                for k in range(0, len(points2d)-1) :
                    lines.append(LineSegment2D(points2d[k], points2d[k+1]))
                points, prim_list = [], []
                for ls_toadd in lines :        
                    same = Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False :
                        prim_list.append(ls_toadd)
                        points.extend(ls_toadd.points)
                if len(points) > 0 : 
                    all_points.extend(points)
                    primitives.append(prim_list)
                    start_end.append([points[0], points[-1]])
        
        # fig, ax = plt.subplots() 
        # [pt.MPLPlot(ax=ax) for pt in all_points]
        # for list_prim in primitives :
        #     fig, ax = plt.subplots() 
        #     [prim.MPLPlot(ax=ax) for prim in list_prim]  
        
        # raise NotImplementedError
        
        points_se, primitives_se = [], []
        for double in start_end : 
            primitives_se.append(LineSegment2D(double[0], double[1]))
            points_se.extend(double)
        poly_se = Polygon2D(points_se)
        
        xmax, xmin = max(pt[0] for pt in points_se), min(pt[0] for pt in points_se)
        ymax, ymin = max(pt[1] for pt in points_se), min(pt[1] for pt in points_se) 
        pt1, pt2, pt3, pt4 = Point2D((xmin, ymin)), Point2D((xmin, ymax)), Point2D((xmax, ymin)), Point2D((xmax, ymax))
        diag1, diag2 = LineSegment2D(pt1, pt4), LineSegment2D(pt2, pt3)
        diag1_cut, diag2_cut = [], []
        diag1_pointcut, diag2_pointcut = [], []
        for enum,l in enumerate(primitives_se) :
            cut1 = diag1.line_intersection(l)
            cut2 = diag2.line_intersection(l)
            if cut1 is not None : 
                diag1_cut.append(enum)
                diag1_pointcut.append(cut1)
            if cut2 is not None : 
                diag2_cut.append(enum)
                diag2_pointcut.append(cut2)
        
        points_common = []
        for enum1,pos1 in enumerate(diag1_cut) :
            for enum2,pos2 in enumerate(diag2_cut) :
                if pos1 == pos2 :
                    points_common.append(primitives_se[pos1].points)
        # fig, ax = plt.subplots() 
        # [l.MPLPlot(ax=ax) for l in primitives_se]
        # [pt.MPLPlot(ax=ax, color='g') for pt in points_se]
        # pt1.MPLPlot(ax=ax, color='r')
        # pt2.MPLPlot(ax=ax, color='r')            
        # pt3.MPLPlot(ax=ax, color='r')
        # pt4.MPLPlot(ax=ax, color='r')   
        # diag1.MPLPlot(ax=ax)
        # diag2.MPLPlot(ax=ax)
        # for couple in points_common :
        #     [pt.MPLPlot(ax=ax, color='b') for pt in couple]        
        
        if len(points_common) >= 1 :
            solve = False
            for couple in points_common :
                if solve :
                    break
                check1, check2 = poly_se.PointBelongs(couple[0]), poly_se.PointBelongs(couple[1])
                start, end = couple[0].vector[0], couple[1].vector[0]
                if math.isclose(start, end, abs_tol = 5e-2) :
                    intersect = min(start, end)
                    if math.isclose(intersect, math.pi, abs_tol = 5e-2) :
                        points_sing = check_singularity(all_points)
                        pt0, pt2pi = 0, 0
                        for pt in points_sing :
                            if math.isclose(pt.vector[0], 0, abs_tol=1e-2):
                                pt0 += 1
                            elif math.isclose(pt.vector[0], 2*math.pi, abs_tol=1e-2):
                                pt2pi += 1
                        points_sing.sort(key=lambda pt: pt[1])
                        points_sing.sort(key=lambda pt: pt[0])
                        if pt2pi !=0 and pt0 == 0 :
                            points = [pt.copy() for pt in points_sing[::-1]]
                            points_sing = points
                        points_range = Face3D.range_closest(points_sing)
                        all_points = delete_double_point(points_range)
                        break
                        
                    elif math.isclose(intersect, 0, abs_tol = 1e-6) or math.isclose(intersect, 2*math.pi, abs_tol = 1e-6) or (not check1 or not check2):
                        all_points = check_singularity(all_points)
                        
                        points_cleaned = delete_double_point(all_points)
                        all_points = [pt.copy() for pt in points_cleaned]
                        all_points.sort(key=lambda pt: pt[0])
                        d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                        if d2 < d1 :
                            last = all_points[-1].copy()
                            all_points[-1] = all_points[-2].copy()
                            all_points[-2] = last
                        break
                    else :
                        points = []
                        for list_prim in primitives :
                            for k,prim in enumerate(list_prim) :
                                new_list_points = []
                                change = 0
                                for enum, pt in enumerate(prim.points) :
                                    if pt[0] < intersect :
                                        change += 1 
                                        if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                            new_list_points.append(Point2D((intersect + 2*math.pi, pt[1])))
                                        else :
                                            new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    elif math.isclose(pt[0], intersect, abs_tol=1e-1) :
                                        change += 1
                                        new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                                    else :
                                        new_list_points.append(pt)
                                if change > 0 :
                                    points.extend(new_list_points)
                                    # list_prim[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                                else :
                                    points.extend(prim.points)
                                    continue
                        points_cleaned = delete_double_point(points)
                        all_points = Face3D.range_trigo(points_cleaned)
                    solve = True
                else :
                    points_cleaned = delete_double_point(all_points)
                    all_points = [pt.copy() for pt in points_cleaned]
                    all_points.sort(key=lambda pt: pt[0])
                    d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
                    if d2 < d1 :
                        last = all_points[-1].copy()
                        all_points[-1] = all_points[-2].copy()
                        all_points[-2] = last
        else :
            points_cleaned = delete_double_point(all_points)
            all_points = [pt.copy() for pt in points_cleaned]
            all_points.sort(key=lambda pt: pt[0])
            d1, d2 = (all_points[0]-all_points[-1]).Norm(), (all_points[0]-all_points[-2]).Norm()
            if d2 < d1 :
                last = all_points[-1].copy()
                all_points[-1] = all_points[-2].copy()
                all_points[-2] = last
        
        primitives = Face3D.create_primitives(all_points)
        
        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax, color='g') for pt in all_points]
        # all_points[0].MPLPlot(ax=ax, color='m')
        # all_points[1].MPLPlot(ax=ax, color='r')
        # all_points[-1].MPLPlot(ax=ax)
        # all_points[-2].MPLPlot(ax=ax, color='b')
        # [p.MPLPlot(ax=ax) for p in primitives]
        
        l_vert = LineSegment2D((pt2+pt4)/2, (pt1+pt3)/2)
        solve = False
        for prim in primitives :
            if solve :
                break
            intersect = prim.line_intersection(l_vert)
            if intersect is not None : 
                x_intersect = intersect.vector[0]
                y_intersect = intersect.vector[1]
                value1, value2 = ymax-0.2*(ymax-ymin), ymin+0.2*(ymax-ymin)
                if y_intersect < max(value1, value2) and y_intersect > min(value1, value2) :
                    points = []
                    for k, prim in enumerate(primitives) :
                        new_list_points, change = [], 0
                        for pt in prim.points :
                            if pt[0] < x_intersect :
                                change += 1
                                if math.isclose(pt[0], 0, abs_tol = 1e-1) :
                                    new_list_points.append(Point2D((x_intersect + 2*math.pi, pt[1])))
                                else :
                                    new_list_points.append(Point2D((2*math.pi + pt[0], pt[1])))
                            else :
                                new_list_points.append(pt)
                        if change > 0:
                            points.extend(new_list_points)
                            primitives[k] = LineSegment2D(new_list_points[0], new_list_points[1])
                        else :
                            points.extend(prim.points)
                            continue
                    solve = True
                    points_cleaned = delete_double_point(points)
                    all_points = Face3D.range_trigo(points_cleaned)
                    primitives = Face3D.create_primitives(all_points)
                else :
                    points_cleaned = delete_double_point(all_points)
                    all_points = Face3D.range_trigo(points_cleaned)
                    primitives = Face3D.create_primitives(all_points)
                    
        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax) for pt in all_points]
        # [p.MPLPlot(ax=ax) for p in primitives]
        # intersect.MPLPlot(ax=ax, color='r')
        # l_vert.MPLPlot(ax=ax)
        contour2d = [Contour2D(primitives)]
        return contour2d
    
    def triangulation(self, resolution=30):
        r = self.radius
        
        frame3d = self.sphericalsurface3d.frame
        
        centerota = Point2D((0,0))
        
        angle_theta = self.points[0]
        angle_phi = self.points[1]
        pas_theta = 2*math.pi/30 #Step of 12 degrees
        pas_phi = pas_theta
        
        resolution_theta = abs(int(angle_theta/pas_theta))
        resolution_phi = abs(int(angle_phi/pas_phi))
        
        if resolution_phi < 5 :
            resolution_phi = 5
        if resolution_theta < 5 :
            resolution_theta = 5
        
        ctr_pt1 = self.cut_contours(self.contours2d, resolution_theta)
        
        
        all_contours_points = []
        for listpt in ctr_pt1 :
            bandept = []
            for pt in listpt :
                bandept.append(pt.Rotation(centerota, -math.pi/2))
            edges = []
            for k in range (0,len(bandept)) :
                if k == len(bandept)-1 :
                    edges.append(LineSegment2D(bandept[k], bandept[0]))
                else :
                    edges.append(LineSegment2D(bandept[k], bandept[k+1]))
            ctr2d = [Contour2D(edges)]
            all_contours_points.extend(self.cut_contours(ctr2d, resolution_phi))
        
        pts_frame1 = []
        for listpt in all_contours_points :
            for pt in listpt : 
                pts_frame1.append(pt.Rotation(centerota, math.pi/2))
        
        all_points = self.delete_double(pts_frame1) # All points necessary to triangulate
        all_points.sort(key=lambda pt: pt[0]) 
        ptvert, pts = [], []
        for k in range(0, resolution_theta +1) :
            ptvert.append([point for point in all_points[k*(resolution_phi+1):(k+1)*(resolution_phi+1)]])
            ptvert[k].sort(key=lambda pt: pt[1])
            pts.extend(ptvert[k])
        
        Triangles, ts = [], []
        if math.isclose(angle_phi, 2*math.pi, abs_tol=1e-6) or math.isclose(angle_theta, 2*math.pi, abs_tol=1e-6):
            step = resolution_phi
            for k in range (0,len(pts)-resolution_phi-2) :
                vertices=[]
                segments=[]
                
                if k%step == 0 and k!=0:
                    step += resolution_phi+1
                    continue
                
                listpt = [pts[k], pts[k+1], pts[k+1+resolution_phi+1], pts[k+resolution_phi+1]]
                listindice = [k, k+1, k+1+resolution_phi+1, k+resolution_phi+1]
                for i, pt in enumerate(listpt):
                    vertices.append(pt.vector)
                    segments.append([i,i+1])
                
                segments[-1]=(len(listpt)-1,0) 
                tri = {'vertices': vertices, 'segments': segments}
                t = triangle.triangulate(tri, 'p')
                if 'triangles' in t:
                    triangles = t['triangles'].tolist()
                    for n,tri in enumerate(triangles):
                        for i in range (0,3):
                            tri[i]=listindice[tri[i]]
                    Triangles.append(triangles)        
                else:
                    Triangles.append(None)
                ts.append(t)
            
            pts3d = self.points2d_to3d(pts, r, frame3d)
            pt3d, tangle = delete_double_pos(pts3d, Triangles)
        
        else :
            A = dict(vertices=[pt.vector for pt in pts]) #in all_points
            B = triangle.triangulate(A, 'cp')
            Triangles = list(B['triangles'])
            # triangle.compare(plt, A, B)
            # print(B['triangles'].shape)
            # print(B)
            pts3d = self.points2d_to3d(pts, r, frame3d)
            pt3d, tangle = delete_double_pos(pts3d, [Triangles])
        
        return pt3d, tangle
         
class BSplineFace3D(Face3D):
    def __init__(self, contours, bspline_shape, points=None, name=''):
        Face3D.__init__(self, contours)
        
        self.contours = contours
        self.bspline_shape = bspline_shape
        self.name = name
        
        self.points2d, self.points3d, self.nbpts_contour2d = self._contour_bspline()
        
    def _contour_bspline(self):
        if self.bspline_shape.__class__ is BSplineExtrusion:
            normal = self.bspline_shape.vectorextru
            
            list_extrudeh = [normal.Dot(Vector3D(pt.vector)) for pt in self.bspline_shape.points]
            #find the bd box's height ,changing the frame of points
            list_h = [normal.Dot(Vector3D(pt.vector)) for pt in self.contours[0].tessel_points[:]]
            mini, maxi = min(list_h), max(list_h)
            
            n = 10
            control_points = []
            for i in range(0,n+1):
                height = mini + (maxi-mini)*i/n
                for enum, pt in enumerate(self.bspline_shape.points) :
                    point_toadd = pt + (height-list_extrudeh[enum])*normal
                    if enum == 0 :
                        pt0 = point_toadd.copy()
                    control_points.append(point_toadd.vector)
                control_points.append(pt0.vector)
            
            surface = BSpline.Surface()
            surface.degree_u = 3 
            surface.degree_v = 1 
            nb_v = len(self.bspline_shape.points)+1
            nb_u = n+1
            surface.set_ctrlpts(control_points, nb_u, nb_v)
            surface.knotvector_u = utilities.generate_knot_vector(surface.degree_u, surface.ctrlpts_size_u)
            surface.knotvector_v = utilities.generate_knot_vector(surface.degree_v, surface.ctrlpts_size_v)
            surface.delta = 0.05
            surface_points = surface.evalpts
            points = [Point3D((vector)) for vector in surface_points]
            
            delta = surface.delta
            nbu, nbv = int(1/delta[0])+1, int(1/delta[1])+1
            if len(self.contours[0].edges) == 4 :
                list_pointsv, list_pointsu = [], []
                for etage in range(0, nbu) :
                    list_pointsv.append(points[etage*nbv:(etage+1)*nbv])
                
                for etage in range(0, nbv) :
                    pointsu = []
                    for k in range(0, nbu) :
                        pointsu.append(points[k*nbv+etage])
                    list_pointsu.append(pointsu)
                    
                points2d, contours2d = [], []
                points3d, contours3d = [], []
                for u, list_pointv in enumerate(list_pointsv) :
                    for v, pt in enumerate(list_pointv) :
                        pt_2d = Point2D((u/(nbu-1), (v/(nbv-1))))
                        if pt_2d.vector[0] == 0 or pt_2d.vector[0] == 1 or pt_2d.vector[1] == 0 or pt_2d.vector[1] == 1 :
                            contours2d.append(pt_2d)
                        else :
                            points2d.append(pt_2d)
                            points3d.append(pt)
                
                range_contours2d = Face3D.range_trigo(contours2d)
                
                all_points2d = range_contours2d + points2d
                
                for pt2d in range_contours2d :
                    u, v = int(pt2d.vector[0]*(nbu-1)), int(pt2d.vector[1]*(nbv-1))
                    contours3d.append(list_pointsv[u][v])
                
                all_points3d = contours3d + points3d
                nbpts_contour2d = len(range_contours2d)
                
            else :
                points_surface, points_contour = points, self.contours[0].tessel_points
                list_pointsv, tessel_points2d = BSplineFace3D.localisation(nbu, nbv, points_surface, points_contour)
                
                poly2d = Polygon2D(tessel_points2d)
                minu, maxu = int(min(pt[0] for pt in tessel_points2d)*(nbu-1)), int(max(pt[0] for pt in tessel_points2d)*(nbu-1))+1
                minv, maxv = int(min(pt[1] for pt in tessel_points2d)*(nbv-1)), int(max(pt[1] for pt in tessel_points2d)*(nbv-1))+1
               
                all_points2d = delete_double_point(tessel_points2d)
                all_points3d = delete_double_point(points_contour)
                nbpts_contour2d = len(points_contour)
                for u, list_pointv in enumerate(list_pointsv[minu:maxu+1]) :
                    for v, pt in enumerate(list_pointv[minv:maxv+1]) :
                        etageu = u + minu
                        etagev = v + minv
                        pt_2d = Point2D((etageu/(nbu-1), (etagev/(nbv-1))))
                        if poly2d.PointBelongs(pt_2d) :
                            all_points2d.append(pt_2d)
                            all_points3d.append(pt)
            # PLOT #
            # from matplotlib import cm
            # from geomdl.visualization import VisMPL
            # surface.vis = VisMPL.VisSurface(ctrlpts=False, legend=False)
            # surface.render(colormap=cm.terrain)
            
            return all_points2d, all_points3d, nbpts_contour2d
            
        elif self.bspline_shape.__class__ is BSplineSurface3D:
            
            delta = self.bspline_shape.surface.delta
            nbu, nbv = int(1/delta[0])+1, int(1/delta[1])+1
            points_surface, points_contour = self.bspline_shape.points, self.contours[0].tessel_points
            list_pointsv, tessel_points2d = BSplineFace3D.localisation(nbu, nbv, points_surface, points_contour)
            
            poly2d = Polygon2D(tessel_points2d)
            minu, maxu = int(min(pt[0] for pt in tessel_points2d)*(nbu-1)), int(max(pt[0] for pt in tessel_points2d)*(nbu-1))+1
            minv, maxv = int(min(pt[1] for pt in tessel_points2d)*(nbv-1)), int(max(pt[1] for pt in tessel_points2d)*(nbv-1))+1
            
            all_points2d = delete_double_point(tessel_points2d)
            all_points3d = delete_double_point(points_contour)
            nbpts_contour2d = len(points_contour)
            for u, list_pointv in enumerate(list_pointsv[minu:maxu+1]) :
                for v, pt in enumerate(list_pointv[minv:maxv+1]) :
                    etageu = u + minu
                    etagev = v + minv
                    pt_2d = Point2D((etageu/(nbu-1), (etagev/(nbv-1))))
                    if poly2d.PointBelongs(pt_2d) :
                        all_points2d.append(pt_2d)
                        all_points3d.append(pt)
            
            return all_points2d, all_points3d, nbpts_contour2d          
    
    def localisation(nbu, nbv, points_surface, points_contour) :
        
        list_pointsv, list_pointsu = [], []
        for etage in range(0, nbu) :
            list_pointsv.append(points_surface[etage*nbv:(etage+1)*nbv])
        
        for etage in range(0, nbv) :
            pointsu = []
            for k in range(0, nbu) :
                pointsu.append(points_surface[k*nbv+etage])
            list_pointsu.append(pointsu)
        
        tessel_points2d = []
        four_points, four_points2d = [], []
        four_distances = []
        for point in points_contour :
            all_posv, distance, global_position = [], [], []
            offset = 0
            for list_pointv in list_pointsv :
                point1, point2 = list_pointv[0], list_pointv[1]
                d = [point.point_distance(point1), point.point_distance(point2)]
                posv = [0, 1]
                glo_posv = [offset+0, offset+1]
                for enum, pt in enumerate(list_pointv[2:]) :
                    dtest = point.point_distance(pt)
                    if dtest < max(d) :
                        maxi, maxpos = max_pos(d)
                        d[maxpos] = dtest
                        posv[maxpos] = enum+2
                        glo_posv[maxpos] = offset+enum+2
                all_posv.append(posv)
                distance.append(d)
                global_position.append(glo_posv)
                offset += len(list_pointv)
            #take 4 closest points
            d = [distance[0], distance[1]]
            avg1, avg2 = sum(d[0])/2, sum(d[1])/2
            if avg1 > avg2 :
                avg = avg1
                pos_to_change = 0
            else :
                avg = avg2
                pos_to_change = 1
            posu = [0, 1]
            for k, dist in enumerate(distance[2:]) :
                avg_test = sum(dist)/2
                if avg_test<avg:
                    d[pos_to_change] = dist
                    posu[pos_to_change] = k+2
                    avg1, avg2 = sum(d[0])/2, sum(d[1])/2
                    if avg1 > avg2 :
                        avg = avg1
                        pos_to_change = 0
                    else :
                        avg = avg2
                        pos_to_change = 1
            point1 = points_surface[global_position[posu[0]][0]]
            point1_2d = Point2D((posu[0]/(nbu-1), all_posv[posu[0]][0]/(nbv-1)))
            d1 = distance[posu[0]][0]
            point2 = points_surface[global_position[posu[0]][1]]
            point2_2d = Point2D((posu[0]/(nbu-1), all_posv[posu[0]][1]/(nbv-1)))
            d2 = distance[posu[0]][1]
            point3 = points_surface[global_position[posu[1]][0]]
            point3_2d = Point2D((posu[1]/(nbu-1), all_posv[posu[1]][0]/(nbv-1)))
            d3 = distance[posu[1]][0]
            point4 = points_surface[global_position[posu[1]][1]]
            point4_2d = Point2D((posu[1]/(nbu-1), all_posv[posu[1]][1]/(nbv-1)))
            d4 = distance[posu[1]][1]
            four_points.append([point1, point2, point3, point4])
            four_points2d.append([point1_2d, point2_2d, point3_2d, point4_2d])
            four_distances.append([d1, d2, d3, d4])
            
            D = sum(four_distances[-1])
            point_2d = (point1_2d*((D-d1)/D) + point2_2d*((D-d2)/D) + point3_2d*((D-d3)/D) + point4_2d*((D-d4)/D))/3
            tessel_points2d.append(point_2d)
        return list_pointsv, tessel_points2d
    
    def triangulation(self):
        seg = []
        for k in range(0, self.nbpts_contour2d) :
            if k == self.nbpts_contour2d-1 :
                seg.append([k,0])
            else :
                seg.append([k,k+1])
        
        A = dict(vertices=[pt.vector for pt in self.points2d],
                 segments=seg) 
        B = triangle.triangulate(A, 'p')
        Triangles = list(B['triangles'])
        # triangle.compare(plt, A, B)
        # print(B['triangles'].shape)
        # print(B)
        
        return self.points3d, [Triangles]
                    
class Shell3D(CompositePrimitive3D):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes  = ['bounding_box']
    _non_eq_attributes = ['name', 'color', 'alpha' 'bounding_box', 'contours']
    _non_hash_attributes = []

    def __init__(self, faces, color=None, alpha=1., name=''):
        self.faces = faces
        self.name = name
        self.color = color
        self.alpha = alpha
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return sum([hash(f) for f in self.faces])

    def __eq__(self, other_):
        if self.__class__ != other_.__class__:
            return False
        equal = True
        for face, other_face in zip(self.faces, other_.faces):
            equal = (equal and face == other_face)
        return equal

    @classmethod
    def from_step(cls, arguments, object_dict):
        faces = []
        for face in arguments[1]:
            faces.append(object_dict[int(face[1:])])
        return cls(faces, name = arguments[0][1:-1])

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_faces = [face.Rotation(center, axis, angle, copy=True) for face in self.faces]
            return Shell3D(new_faces, name = self.name)
        else:
            for face in self.faces:
                face.Rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_faces = [face.Translation(offset, copy=True) for face in self.faces]
            return Shell3D(new_faces, name = self.name)
        else:
            for face in self.faces:
                face.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_faces = [face.frame_mapping(frame, side, copy=True) for face in self.faces]
            return Shell3D(new_faces, name = self.name)
        else:
            for face in self.faces:
                face.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        new_faces = [face.copy() for face in self.faces]
        return Shell3D(new_faces, name = self.name)

    def union(self, shell2):
        new_faces = [face for face in self.faces+shell2.faces]
        new_name = self.name+' union '+shell2.name
        new_color = self.color
        return Shell3D(new_faces, name = new_name, color = new_color)
    
    def volume(self):
        """
        Do not consider holes
        """
        volume = 0
        for i, face in enumerate(self.faces):
            points_3D, triangles_indexes = face.triangulation()
            for triangle_indexes in triangles_indexes[0]:
                
                point1 = points_3D[triangle_indexes[0]]
                point2 = points_3D[triangle_indexes[1]]
                point3 = points_3D[triangle_indexes[2]]
                
                v321 = point3[0] * point2[1] * point1[2]
                v231 = point2[0] * point3[1] * point1[2]
                v312 = point3[0] * point1[1] * point2[2]
                v132 = point1[0] * point3[1] * point2[2]
                v213 = point2[0] * point1[1] * point3[2]
                v123 = point1[0] * point2[1] * point3[2]
                volume_tetraedre = 1/6 * (-v321 + v231 + v312 - v132 - v213 + v123)
                # print(volume_tetraedre)
                
                volume += volume_tetraedre
        
        return abs(volume)
                
    def _bounding_box(self):
        """
        Returns the boundary box
        """
        points = []
        for face in self.faces:
            points.extend(face.bounding_box.points)

        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])

        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax, self.name)

    def point_belongs(self, point, nb_rays=1):
        """
        Ray Casting algorithm
        Returns True if the point is inside the Shell, False otherwise
        """
        epsilon = 10

        bbox = self.bounding_box
        # print('bounding_box',bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax, bbox.zmin, bbox.zmax)
        # print(point)
        if point[0] < bbox.xmin or point[0] > bbox.xmax:
            return False
        if point[1] < bbox.ymin or point[1] > bbox.ymax:
            return False
        if point[2] < bbox.zmin or point[2] > bbox.zmax:
            return False

        rays = []
        for k in range (0, nb_rays):
            rays.append(LineSegment3D(point,Point3D((random.uniform(0, 1)*epsilon, random.uniform(0, 1)*epsilon, random.uniform(0, 1)*epsilon))))

        rays = sorted(rays, key=lambda ray: ray.Length())

        rays_intersections = []
        tests = []
        
        # for ray in rays[:3]:
        for ray in rays[:nb_rays]:
            count = 0
            ray_intersection = []
            is_inside = True
            for face in self.faces:
                intersection_point = face.linesegment_intersection(ray)
                if intersection_point is not None:
                    ray_intersection.append(intersection_point)
                    count += 1
            if count%2 == 0:
                is_inside = False
            tests.append(is_inside)
            rays_intersections.append(ray_intersection)
        
        for test1, test2 in zip(tests[:-1], tests[1:]):
            if test1 != test2:
                raise ValueError

        return tests[0]

    def is_inside_shell(self, shell2):
        """
        Returns True if all the points of self are inside shell2 and no face \
        are intersecting
        """
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.is_inside_bbox(bbox2):
            return False

        points = []
        for face in self.faces:
            points.extend(face.contours3d[0].tessel_points)

        for point in points:
            if not shell2.point_belongs(point):
                return False

        # Check if any faces are intersecting
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
#                    print('Two faces are intersecting :', face1, face2)
                    return False

        return True

    def shell_intersection(self, shell2):
        """
        Return None if disjointed
        Return (1, 0) or (0, 1) if one is inside the other
        Return (n1, n2) if intersection

        4 cases :
            (n1, n2) with face intersection             => (n1, n2)
            (0, 0) with face intersection               => (0, 0)
            (0, 0) with no face intersection            => None
            (1, 0) or (0, 1) with no face intersection  => 1
        """
        # Check if boundary boxes don't intersect
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
#            print("No intersection of shells' BBox")
            return None

        # Check if any point of the first shell is in the second shell
        points1 = []
        for face in self.faces:
            points1.extend(face.contours3d[0].tessel_points)
        points2 = []
        for face in shell2.faces:
            points2.extend(face.contours3d[0].tessel_points)

        nb_pts1 = len(points1)
        nb_pts2 = len(points2)
        compteur1 = 0
        compteur2 = 0
        for point1 in points1:
            if shell2.point_belongs(point1):
                compteur1 += 1
        for point2 in points2:
            if self.point_belongs(point2):
                compteur2 += 1

        inter1 = compteur1/nb_pts1
        inter2 = compteur2/nb_pts2
#        print('shell intersection')
#        print('shell1 intersecte shell2 à', inter1*100, '%')
#        print('shell2 intersecte shell1 à', inter2*100, '%')

        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
#                    print('Two faces are intersecting :', face1, face2)
#                    ax = face1.plot()
#                    face2.plot(ax)
                    return inter1, inter2
        if (inter1, inter2) == (0, 0):
            return None
        return 1

    def distance_to_shell(self, shell2, add_to_volumemodel=None):
        """
        Returns a Mesure object if the distance is not zero, otherwise returns None
        """

        if self.shell_intersection(shell2) is not None and self.shell_intersection(shell2) != 1:
            return None

        # distance_min, point1_min, point2_min = self.faces[0].distance_to_face(shell2.faces[0], return_points=True)
        distance_min, point1_min, point2_min = self.faces[0].minimum_distance(shell2.faces[0], return_points=True)
        for face1 in self.faces:
            bbox1 = face1.bounding_box
            for face2 in shell2.faces:
                bbox2 = face2.bounding_box
                bbox_distance = bbox1.distance_to_bbox(bbox2)
                if bbox_distance < distance_min:
                    # distance, point1, point2 = face1.distance_to_face(face2, return_points=True)
                    distance, point1, point2 = face1.minimum_distance(face2, return_points=True)
                    if distance == 0:
                        return None
                    elif distance < distance_min:
                        distance_min, point1_min, point2_min = distance, point1, point2

        mesure = Measure3D(point1_min, point2_min)

        if add_to_volumemodel is not None:
            add_to_volumemodel.primitives.append(mesure)

        return mesure

    def distance_to_point(self, point, add_to_volumemodel=None):
        """
        Computes the distance of a point to a Shell3D, whether it is inside or outside the Shell3D
        """
        distance_min, point1_min = self.faces[0].distance_to_point(point, return_other_point=True)
        for face in self.faces[1:]:
            bbox_distance = self.bounding_box.distance_to_point(point)
            if bbox_distance < distance_min:
                distance, point1 = face.distance_to_point(point, return_other_point=True)
                if distance < distance_min:
                    distance_min, point1_min = distance, point1

        mesure = Measure3D(point, point1_min)

        if add_to_volumemodel is not None:
            add_to_volumemodel.primitives.append(mesure)

        return mesure

    def intersection_internal_aabb_volume(self, shell2):
        """
        aabb made of the intersection points and the points of self internal to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_inside_shell2 = []
        for face in self.faces:
            for point in face.contours3d[0].tessel_points:
                if shell2.point_belongs(point):
                    shell1_points_inside_shell2.append(point)

        if len(intersections_points+shell1_points_inside_shell2) == 0:
            return 0
        bbox = BoundingBox.from_points(intersections_points+shell1_points_inside_shell2)
        return bbox.volume()

    def intersection_external_aabb_volume(self, shell2):
        """
        aabb made of the intersection points and the points of self external to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_outside_shell2 = []
        for face in self.faces:
            for point in face.contours3d[0].tessel_points:
                if not shell2.point_belongs(point):
                    shell1_points_outside_shell2.append(point)

        if len(intersections_points+shell1_points_outside_shell2) == 0:
            return 0
        bbox = BoundingBox.from_points(intersections_points+shell1_points_outside_shell2)
        return bbox.volume()

    # def babylon_meshes(self):
    def triangulation(self):
        positions = []
        indices = []

        nb_points = 0
        for i, face in enumerate(self.faces):
            # print('\n face', face)
            
            # print('face.triangulation', face.triangulation())
            # print('i', i)
            # print()
            points_3D, triangles_indexes = face.triangulation()
            # points_3D_triangles_indexes = face.triangulation()
            # print('=>', points_3D, triangles_indexes)
            # print()
            # print('len pt3d',len(points_3D))
            
            if points_3D is None : 
                continue
            
            
            for point in points_3D:
                # print('===========>',point)
                positions.extend([k for k in round(point, 6)])
            # print('len positions',len(positions))
            for j, indexes in enumerate(triangles_indexes):
                # print('j',j)
                # print('indexes',indexes)
                # print()
                for index in indexes:
                    indices.extend([k+nb_points for k in index])
            nb_points += len(points_3D)
        return positions, indices
    
    def babylon_meshes(self):
        positions, indices = self.triangulation()
        babylon_mesh = {'positions': positions,
                        'indices': indices,
                        'alpha': self.alpha,
                        'name': self.name
                        }
        
        if self.color is None:
            babylon_mesh['color'] = [0.8, 0.8, 0.8]
        else:
            babylon_mesh['color'] = list(self.color)
        
        return [babylon_mesh]

    def babylon_script(self, name='primitive_mesh'):
        s = 'var {} = new BABYLON.Mesh("{}", scene);\n'.format(name, name)

        mesh = self.babylon_meshes()[0]
        
        s += 'var positions = {};\n'.format(mesh['positions'])
        s += 'var indices = {};\n'.format(mesh['indices'])
        s += 'var normals = [];\n'
        s += 'var vertexData = new BABYLON.VertexData();\n'
        s += 'BABYLON.VertexData.ComputeNormals(positions, indices, normals);\n'
        s += 'vertexData.positions = positions;\n'
        s += 'vertexData.indices = indices;\n'
        s += 'vertexData.normals = normals;\n'
        s += 'vertexData.applyToMesh({});\n'.format(name)
        s += '{}.enableEdgesRendering(0.9);\n'.format(name)
        s += '{}.edgesWidth = 0.1;\n'.format(name)
        s += '{}.edgesColor = new BABYLON.Color4(0, 0, 0, 0.6);\n'.format(name)
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
#        s += 'mat.diffuseColor = BABYLON.Color3.Green();\n'
#        s += 'mat.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);\n'
#        s += 'mat.emissiveColor = new BABYLON.Color3(1, 1, 1);\n'
#        s += 'mat.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);\n'
        s += 'mat.backFaceCulling = false;\n'
        s += 'mat.alpha = {};\n'.format(self.alpha)
        s += '{}.material = mat;\n'.format(name)
        if self.color is not None:
            s += 'mat.diffuseColor = new BABYLON.Color3({}, {}, {});\n'.format(*self.color)
        return s

class BoundingBox(dc.DessiaObject):
    """
    An axis aligned boundary box
    """
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax, name=''):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.points = [Point3D((self.xmin, self.ymin, self.zmin)), \
                       Point3D((self.xmax, self.ymin, self.zmin)), \
                       Point3D((self.xmax, self.ymax, self.zmin)), \
                       Point3D((self.xmin, self.ymax, self.zmin)), \
                       Point3D((self.xmin, self.ymin, self.zmax)), \
                       Point3D((self.xmax, self.ymin, self.zmax)), \
                       Point3D((self.xmax, self.ymax, self.zmax)), \
                       Point3D((self.xmin, self.ymax, self.zmax))]
        self.center = (self.points[0]+self.points[-2])/2
        self.name = name

    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def plot(self, ax=None, color=''):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        bbox_edges =  [[self.points[0], self.points[1]], [self.points[0], self.points[3]], \
                       [self.points[0], self.points[4]], [self.points[1], self.points[2]], \
                       [self.points[1], self.points[5]], [self.points[2], self.points[3]], \
                       [self.points[2], self.points[6]], [self.points[3], self.points[7]], \
                       [self.points[4], self.points[5]], [self.points[5], self.points[6]], \
                       [self.points[6], self.points[7]], [self.points[7], self.points[4]]]

        x = [p[0] for p in self.points]
        y = [p[1] for p in self.points]
        z = [p[2] for p in self.points]
        ax.scatter(x, y, z)
        for edge in bbox_edges:
            ax.plot3D([edge[0][0], edge[1][0]],
                      [edge[0][1], edge[1][1]],
                      [edge[0][2], edge[1][2]],
                      'gray')
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.show()
        return ax

    @classmethod
    def from_points(cls, points):
        # if len(points) == 0:
        #     return (0, 0, 0, 0, 0, 0)
        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)

    def volume(self):
        return (self.xmax-self.xmin)*(self.ymax-self.ymin)*(self.zmax-self.zmin)

    def bbox_intersection(self, bbox2):
        return (self.xmin < bbox2.xmax and self.xmax > bbox2.xmin \
                and self.ymin < bbox2.ymax and self.ymax > bbox2.ymin \
                and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin)

    def is_inside_bbox(self, bbox2):
        return (self.xmin > bbox2.xmin and self.xmax < bbox2.xmax \
                and self.ymin > bbox2.ymin and self.ymax < bbox2.ymax \
                and self.zmin > bbox2.zmin and self.zmax < bbox2.zmax)

    def intersection_volume(self, bbox2):
        if not self.bbox_intersection(bbox2):
            return 0
        if self.is_inside_bbox(bbox2) or bbox2.is_inside_bbox(self):
            return min(self.volume(), bbox2.volume())

        lx = min(self.xmax, bbox2.xmax) - max(self.xmin, bbox2.xmin)
        ly = min(self.ymax, bbox2.ymax) - max(self.ymin, bbox2.ymin)
        lz = min(self.zmax, bbox2.zmax) - max(self.zmin, bbox2.zmin)

        return lx * ly * lz

    # def intersection_volume(self, bbox2):
    #     if not self.bbox_intersection(bbox2):
    #         return 0
    #
    #     permute_bbox1 = self
    #     permute_bbox2 = bbox2
    #
    #     if permute_bbox2.xmin < permute_bbox1.xmin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     lx = permute_bbox1.xmax - permute_bbox2.xmin
    #
    #     if permute_bbox2.ymin < permute_bbox1.ymin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     ly = permute_bbox1.ymax - permute_bbox2.ymin
    #
    #     if permute_bbox2.zmin < permute_bbox1.zmin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     lz = permute_bbox1.zmax - permute_bbox2.zmin
    #
    #     return lx*ly*lz

    def distance_to_bbox(self, bbox2):
        if self.bbox_intersection(bbox2):
            return 0

        permute_bbox1 = self
        permute_bbox2 = bbox2

        if permute_bbox2.xmin < permute_bbox1.xmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dx = permute_bbox2.xmin - permute_bbox1.xmax
        if dx < 0:
            dx = 0

        if permute_bbox2.ymin < permute_bbox1.ymin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dy = permute_bbox2.ymin - permute_bbox1.ymax
        if dy < 0:
            dy = 0

        if permute_bbox2.zmin < permute_bbox1.zmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dz = permute_bbox2.zmin - permute_bbox1.zmax
        if dz < 0:
            dz = 0

        return (dx**2 + dy**2 + dz**2)**0.5

    def point_belongs(self, point):
        return self.xmin < point[0] and point[0] < self.xmax \
        and self.ymin < point[1] and point[1] < self.ymax \
        and self.zmin < point[2] and point[2] < self.zmax

    def distance_to_point(self, point):
        if self.point_belongs(point):
            return min([self.xmax-point[0], point[0]-self.xmin,
                        self.ymax-point[1], point[1]-self.ymin,
                        self.zmax-point[2], point[2]-self.zmin])
        else:
            if point[0] < self.xmin:
                dx = self.xmin - point[0]
            elif self.xmax < point[0]:
                dx = point[0] - self.xmax
            else:
                dx = 0

            if point[1] < self.ymin:
                dy = self.ymin - point[1]
            elif self.ymax < point[1]:
                dy = point[1] - self.ymax
            else:
                dy = 0

            if point[2] < self.zmin:
                dz = self.zmin - point[2]
            elif self.zmax < point[2]:
                dz = point[2] - self.zmax
            else:
                dz = 0
        return (dx**2 + dy**2 + dz**2)**0.5

    def distance_between_two_points_on_bbox(self, point1, point2):

        if   math.isclose(point1[0], self.xmin, abs_tol=1e-8):
            face_point1 = 5
        elif math.isclose(point1[0], self.xmax, abs_tol=1e-8):
            face_point1 = 3
        elif math.isclose(point1[1], self.ymin, abs_tol=1e-8):
            face_point1 = 4
        elif math.isclose(point1[1], self.ymax, abs_tol=1e-8):
            face_point1 = 2
        elif math.isclose(point1[2], self.zmin, abs_tol=1e-8):
            face_point1 = 6
        elif math.isclose(point1[2], self.zmax, abs_tol=1e-8):
            face_point1 = 1
        else:
            raise NotImplementedError

        if   math.isclose(point2[0], self.xmin, abs_tol=1e-8):
            face_point2 = 5
        elif math.isclose(point2[0], self.xmax, abs_tol=1e-8):
            face_point2 = 3
        elif math.isclose(point2[1], self.ymin, abs_tol=1e-8):
            face_point2 = 4
        elif math.isclose(point2[1], self.ymax, abs_tol=1e-8):
            face_point2 = 2
        elif math.isclose(point2[2], self.zmin, abs_tol=1e-8):
            face_point2 = 6
        elif math.isclose(point2[2], self.zmax, abs_tol=1e-8):
            face_point2 = 1
        else:
            raise NotImplementedError

        point1_copy = point1.copy()
        point2_copy = point2.copy()
#        print(face_point1, face_point2)
        if face_point1 > face_point2:
#            print('inversion')
            point1, point2 = point2, point1
            face_point1, face_point2 = face_point2, face_point1

        # The points are on the same face
        if face_point1 == face_point2:
            return point1.point_distance(point2)

        deltax = self.xmax - self.xmin
        deltay = self.ymax - self.ymin
        deltaz = self.zmax - self.zmin

        point1_2d_coordinate_dict = {1: Point2D((point1[0]-self.xmin-deltax/2, point1[1]-self.ymin-deltay/2)),
                                     2: Point2D((point1[2]-self.zmin-deltaz/2, point1[0]-self.xmin-deltax/2)),
                                     3: Point2D((point1[1]-self.ymin-deltay/2, point1[2]-self.zmin-deltaz/2)),
                                     4: Point2D((point1[0]-self.xmin-deltax/2, point1[2]-self.zmin-deltaz/2)),
                                     5: Point2D((point1[2]-self.zmin-deltaz/2, point1[1]-self.ymin-deltay/2)),
                                     6: Point2D((point1[1]-self.ymin-deltay/2, point1[0]-self.xmin-deltax/2))}

        point2_2d_coordinate_dict = {1: Point2D((point2[0]-self.xmin-deltax/2, point2[1]-self.ymin-deltay/2)),
                                     2: Point2D((point2[2]-self.zmin-deltaz/2, point2[0]-self.xmin-deltax/2)),
                                     3: Point2D((point2[1]-self.ymin-deltay/2, point2[2]-self.zmin-deltaz/2)),
                                     4: Point2D((point2[0]-self.xmin-deltax/2, point2[2]-self.zmin-deltaz/2)),
                                     5: Point2D((point2[2]-self.zmin-deltaz/2, point2[1]-self.ymin-deltay/2)),
                                     6: Point2D((point2[1]-self.ymin-deltay/2, point2[0]-self.xmin-deltax/2))}

        vertex_2d_coordinate_dict = {1: [Point2D((self.xmin-self.xmin-deltax/2, self.ymin-self.ymin-deltay/2)), Point2D((self.xmin-self.xmin-deltax/2, self.ymax-self.ymin-deltay/2)), Point2D((self.xmax-self.xmin-deltax/2, self.ymax-self.ymin-deltay/2)), Point2D((self.xmax-self.xmin-deltax/2, self.ymin-self.ymin-deltay/2))],
                                     2: [Point2D((self.zmin-self.zmin-deltaz/2, self.xmin-self.xmin-deltax/2)), Point2D((self.zmin-self.zmin-deltaz/2, self.xmax-self.xmin-deltax/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.xmax-self.xmin-deltax/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.xmin-self.xmin-deltax/2))],
                                     3: [Point2D((self.ymin-self.ymin-deltay/2, self.zmin-self.zmin-deltaz/2)), Point2D((self.ymin-self.ymin-deltay/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.ymax-self.ymin-deltay/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.ymax-self.ymin-deltay/2, self.zmin-self.zmin-deltaz/2))],
                                     4: [Point2D((self.xmin-self.xmin-deltax/2, self.zmin-self.zmin-deltaz/2)), Point2D((self.xmin-self.xmin-deltax/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.xmax-self.xmin-deltax/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.xmax-self.xmin-deltax/2, self.zmin-self.zmin-deltaz/2))],
                                     5: [Point2D((self.zmin-self.zmin-deltaz/2, self.ymin-self.ymin-deltay/2)), Point2D((self.zmin-self.zmin-deltaz/2, self.ymax-self.ymin-deltay/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.ymax-self.ymin-deltay/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.ymin-self.ymin-deltay/2))],
                                     6: [Point2D((self.ymin-self.ymin-deltay/2, self.xmin-self.xmin-deltax/2)), Point2D((self.ymin-self.ymin-deltay/2, self.xmax-self.xmin-deltax/2)), Point2D((self.ymax-self.ymin-deltay/2, self.xmax-self.xmin-deltax/2)), Point2D((self.ymax-self.ymin-deltay/2, self.xmin-self.xmin-deltax/2))],}

#        vertex_2d_coordinate_dict2 = {1: [Point2D((self.xmin, self.ymin)), Point2D((self.xmin, self.ymax)), Point2D((self.xmax, self.ymax)), Point2D((self.xmax, self.ymin))],
#                                    2: [Point2D((self.zmin, self.xmin)), Point2D((self.zmin, self.xmax)), Point2D((self.zmax, self.xmax)), Point2D((self.zmax, self.xmin))],
#                                     3: [Point2D((self.ymin, self.zmin)), Point2D((self.ymin, self.zmax)), Point2D((self.ymax, self.zmax)), Point2D((self.ymax, self.zmin))],
#                                     4: [Point2D((self.xmin, self.zmin)), Point2D((self.xmin, self.zmax)), Point2D((self.xmax, self.zmax)), Point2D((self.xmax, self.zmin))],
#                                     5: [Point2D((self.zmin, self.ymin)), Point2D((self.zmin, self.ymax)), Point2D((self.zmax, self.ymax)), Point2D((self.zmax, self.ymin))],
#                                     6: [Point2D((self.ymin, self.xmin)), Point2D((self.ymin, self.xmax)), Point2D((self.ymax, self.xmax)), Point2D((self.ymax, self.xmin))],}


        vertex_to_3d_dict = {1: (2, self.zmax, 0, 1),
                             2: (1, self.ymax, 2, 0),
                             3: (0, self.xmax, 1, 2),
                             4: (1, self.ymin, 0, 2),
                             5: (0, self.xmin, 2, 1),
                             6: (2, self.zmin, 1, 0)}

        offset_dict = {0: self.xmin+deltax/2,
                       1: self.ymin+deltay/2,
                       2: self.zmin+deltaz/2}

        opposite_face_dict = {1: 6, 2: 4, 3: 5, 4: 2, 5: 3, 6: 1}

        combination_dict = {(1, 2): Frame2D(Point2D((0,  deltay/2+deltaz/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (2, 1): Frame2D(Point2D(( deltay/2+deltaz/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (1, 3): Frame2D(Point2D(( deltax/2+deltaz/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (3, 1): Frame2D(Point2D((0,  deltax/2+deltaz/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (1, 4): Frame2D(Point2D((0, -deltay/2-deltaz/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (4, 1): Frame2D(Point2D((-deltay/2-deltaz/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (1, 5): Frame2D(Point2D((-deltax/2-deltaz/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (5, 1): Frame2D(Point2D((0, -deltax/2-deltaz/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (2, 3): Frame2D(Point2D((0,  deltax/2+deltay/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (3, 2): Frame2D(Point2D(( deltax/2+deltay/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (2, 5): Frame2D(Point2D((0, -deltax/2-deltay/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (5, 2): Frame2D(Point2D((-deltax/2-deltay/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (2, 6): Frame2D(Point2D((-deltaz/2-deltay/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (6, 2): Frame2D(Point2D((0, -deltaz/2-deltay/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (3, 4): Frame2D(Point2D((-deltay/2-deltax/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (4, 3): Frame2D(Point2D((0, -deltay/2-deltax/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (3, 6): Frame2D(Point2D((0, -deltaz/2-deltax/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (6, 3): Frame2D(Point2D((-deltaz/2-deltax/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (4, 5): Frame2D(Point2D((-deltax/2-deltay/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (5, 4): Frame2D(Point2D((0, -deltax/2-deltay/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (4, 6): Frame2D(Point2D((0, -deltaz/2-deltay/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (6, 4): Frame2D(Point2D((-deltaz/2-deltay/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (5, 6): Frame2D(Point2D((-deltaz/2-deltax/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (6, 5): Frame2D(Point2D((0, -deltaz/2-deltax/2)), Vector2D((0,-1)), Vector2D(( 1, 0)))}

        point1_2d = point1_2d_coordinate_dict[face_point1]
        point2_2d = point2_2d_coordinate_dict[face_point2]

        # The points are on adjacent faces
        if opposite_face_dict[face_point1] != face_point2:
            frame =  combination_dict[(face_point1, face_point2)]
            net_point2 = frame.OldCoordinates(point2_2d)

            # Computes the 3D intersection between the net_line and the edges of the face_point1
            net_line = LineSegment2D(point1_2d, net_point2)
            vertex_points = vertex_2d_coordinate_dict[face_point1]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in zip(vertex_points, vertex_points[1:]+[vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point, a, b = Point2D.LinesIntersection(net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point is not None \
                and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[face_point1]
            disordered_coordinate = [(indice1, edge_intersection_point[0]+offset_dict[indice1]),
                                     (indice2, edge_intersection_point[1]+offset_dict[indice2]),
                                     (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate, key=lambda a: a[0])
            intersection_point_3d = Point3D(tuple([p[1] for p in disordered_coordinate]))

            mesures = [Measure3D(point1_copy, intersection_point_3d),
                       Measure3D(intersection_point_3d, point2_copy)]

            return mesures

        # The points are on opposite faces
        else:
            net_points2_and_frame = []

            faces_number = [1, 2, 3, 4, 5, 6]
            faces_number.remove(face_point1)
            faces_number.remove(face_point2)
            pathes = []
            for face_nb in faces_number:
                path = [(face_point1, face_nb), (face_nb, face_point2)]
                pathes.append(path)

            for path in pathes:
                frame1 = combination_dict[(path[0][0], path[0][1])]
                frame2 = combination_dict[(path[1][0], path[1][1])]
                frame = frame1 + frame2
                net_points2_and_frame.append((Point2D(frame.OldCoordinates(point2_2d).vector), frame))
            net_point2, frame = min(net_points2_and_frame, key=lambda pt: pt[0].point_distance(point1_2d))
            net_line = LineSegment2D(point1_2d, net_point2)

            # Computes the 3D intersection between the net_line and the edges of the face_point1
            vertex_points = vertex_2d_coordinate_dict[face_point1]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in zip(vertex_points, vertex_points[1:]+[vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point1, a, b = Point2D.LinesIntersection(net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point1 is not None \
                and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[face_point1]
            disordered_coordinate = [(indice1, edge_intersection_point1[0]+offset_dict[indice1]),
                                     (indice2, edge_intersection_point1[1]+offset_dict[indice2]),
                                     (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate, key=lambda a: a[0])
            intersection_point1_3d = Point3D(tuple([p[1] for p in disordered_coordinate]))

            # Computes the 3D intersection between the net_line and the edges of the face_point2
            vertex_points = [frame.OldCoordinates(p) for p in vertex_2d_coordinate_dict[face_point2]]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in zip(vertex_points, vertex_points[1:]+[vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point2, a, b = Point2D.LinesIntersection(net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point2 is not None \
                and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            edge_intersection_point2 = Point2D(frame.NewCoordinates(edge_intersection_point2))
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[face_point2]
            disordered_coordinate = [(indice1, edge_intersection_point2[0]+offset_dict[indice1]),
                                     (indice2, edge_intersection_point2[1]+offset_dict[indice2]),
                                     (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate, key=lambda a: a[0])
            intersection_point2_3d = Point3D(tuple([p[1] for p in disordered_coordinate]))

            if point1 == point1_copy:
                mesures = [Measure3D(point1, intersection_point1_3d),
                           Measure3D(intersection_point1_3d, intersection_point2_3d),
                           Measure3D(intersection_point2_3d, point2)]
            else:
                mesures = [Measure3D(point2, intersection_point2_3d),
                           Measure3D(intersection_point2_3d, intersection_point1_3d),
                           Measure3D(intersection_point1_3d, point1)]
            return mesures

    def babylon_script(self):
        height = self.ymax-self.ymin
        width = self.xmax-self.xmin
        depth = self.zmax-self.zmin
        s = 'var box = BABYLON.MeshBuilder.CreateBox("box", {{height: {}, width: {}, depth: {}}}, scene);\n'.format(height, width, depth)
        s += 'box.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n'.format(self.center[0], self.center[1], self.center[2])
        s += 'var bboxmat = new BABYLON.StandardMaterial("bboxmat", scene);\n'
        s += 'bboxmat.alpha = 0.4;\n'
        s += 'var DTWidth = {};\n'.format(width*60)
        s += 'var DTHeight = {};\n'.format(height*60)
        s += 'var font_type = "Arial";\n'
        s += 'var text = "{}";\n'.format(self.name)
        s += 'var dynamicTexture = new BABYLON.DynamicTexture("DynamicTexture", {width:DTWidth, height:DTHeight}, scene);\n'
        s += 'var ctx = dynamicTexture.getContext();\n'
        s += 'var size = 0.8;\n'
        s += 'ctx.font = size + "px " + font_type;\n'
        s += 'var textWidth = ctx.measureText(text).width;\n'
        s += 'var ratio = textWidth/size;\n'
        s += 'var font_size = Math.floor(DTWidth / ratio);\n'
        s += 'var font = font_size + "px " + font_type;\n'
        s += 'dynamicTexture.drawText(text, null, null, font, "#000000", "#ffffff", false);\n'
        s += 'bboxmat.diffuseTexture = dynamicTexture;\n'
        s += 'box.material = bboxmat;\n'
        return s

class Measure2D(LineSegment2D):
    def __init__(self, point1, point2, label='', unit='mm', type_='distance'):
        """
        :param unit: 'mm', 'm' or None. If None, the distance won't be in the label

        """
        # TODO: offset parameter
        LineSegment2D.__init__(self, point1, point2)
        self.label = label
        self.unit = unit
        self.type_ = type_
        
    def MPLPlot(self, ax, ndigits=6):
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]
        xm, ym = 0.5*(self.points[0] + self.points[1])
        distance = self.points[1].point_distance(self.points[0])
        
        if self.label != '':
            label = '{}: '.format(self.label)
        else:
            label = ''
        if self.unit == 'mm':            
            label += '{} mm'.format(round(distance*1000, ndigits))
        else:
            label += '{} m'.format(round(distance, ndigits))
        
        if self.type_ == 'distance':
            arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='<|-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')
        elif self.type_ == 'radius':
            arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')
            

        ax.add_patch(arrow)
        if x2-x1 == 0.:
            theta = 90.
        else:            
            theta = math.degrees(math.atan((y2-y1)/(x2-x1)))
        ax.text(xm, ym, label, va='bottom', ha='center', rotation=theta)
        

class Measure3D(Line3D):
    def __init__(self, point1, point2, color=(1,0,0)):
        self.points = [point1, point2]
        self.color = color
        self.distance = Vector3D(self.points[0]-self.points[1]).Norm()
        self.bounding_box = self._bounding_box()

    # !!! no eq defined!
    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def babylon_script(self):
        s = 'var myPoints = [];\n'
        s += 'var point1 = new BABYLON.Vector3({},{},{});\n'.format(self.points[0][0],self.points[0][1],self.points[0][2])
        s += 'myPoints.push(point1);\n'
        s += 'var point2 = new BABYLON.Vector3({},{},{});\n'.format(self.points[1][0],self.points[1][1],self.points[1][2])
        s += 'myPoints.push(point2);\n'
        s += 'var line = BABYLON.MeshBuilder.CreateLines("lines", {points: myPoints}, scene);\n'
        s += 'line.color = new BABYLON.Color3({}, {}, {});\n'.format(self.color[0], self.color[1], self.color[2])
        return s

class Group:
    def __init__(self, primitives, name):
        self.primitives = primitives
        self.name = name


class StepFunction:
    def __init__(self, function_id, function_name, function_arg):
        self.id = function_id
        self.name = function_name
        self.arg = function_arg

        if self.name == "":
            if self.arg[1][0] == 'B_SPLINE_SURFACE':
                self.simplify('B_SPLINE_SURFACE')
            if self.arg[1][0] == 'B_SPLINE_CURVE':
                self.simplify('B_SPLINE_CURVE')
                
    def simplify(self, new_name):
        # ITERATE ON SUBFUNCTIONS
        args = [subfun[1] for (i, subfun) in enumerate(self.arg) if (len(subfun[1]) != 0 or i == 0)]
        arguments = []
        for arg in args:
            if arg == []:
                arguments.append("''")
            else:
                arguments.extend(arg)
        arguments.pop() # DELETE REPRESENTATION_ITEM('')

        self.name = new_name
        self.arg = arguments
        

class Step:
    def __init__(self, stepfile):
        self.stepfile = stepfile

        self.functions, self.all_connections = self.read_functions()

    def read_functions(self):
        f = open(self.stepfile, "r", encoding = "ISO-8859-1")

        all_connections = []

        previous_line = ""
        functions = {}

        for line in f:

            line = line.replace(" ", "")
            line = line.replace("\n", "")

            # SKIP EMPTY LINE
            if not line:
                continue

            # ASSEMBLE LINES IF THEY ARE SEPARATED
            if line[-1] != ';':
                previous_line = previous_line + line
                continue

            line = previous_line + line

            # SKIP HEADER
            if line[0] != "#":
                previous_line = str()
                continue

            function = line.split("=")
            function_id = int(function[0][1:])
            function_name_arg = function[1].split("(", 1)
            function_name = function_name_arg[0]
            function_arg = function_name_arg[1].split("#")
            function_connections = []
            for connec in function_arg[1:]:
                connec = connec.split(",")
                connec = connec[0].split(")")
                if connec[0][-1] != "'":
                    function_connection = int(connec[0])
                    function_connections.append((function_id, function_connection))

            all_connections.extend(function_connections)

            previous_line = str()

            # FUNCTION ARGUMENTS
            function_arg = function_name_arg[1]
            arguments = step_split_arguments(function_arg)
            if function_name == "":
                arguments = self.step_subfunctions(arguments)

            for i, argument in enumerate(arguments):
                if argument[:2] == '(#' and argument[-1] == ')':
                    arg_list = set_to_list(argument)
                    arguments[i] = arg_list

            function = StepFunction(function_id, function_name, arguments)
            functions[function_id] = function

        f.close()

        return functions, all_connections

    def create_graph(self, draw=False, html=False):

        G = nx.Graph()
        F = nx.DiGraph()
        labels = {}

        for function in self.functions.values():
            if function.name in step_to_volmdlr_primitive:
                G.add_node(function.id)
                F.add_node(function.id)
                labels[function.id] = str(function.id)+' '+function.name

        # Delete connection if node not found
        node_list = list(F.nodes())
        delete_connection = []
        for connection in self.all_connections:
            if connection[0] not in node_list or connection[1] not in node_list:
                delete_connection.append(connection)
        for delete in delete_connection:
            self.all_connections.remove(delete)

        # Create graph connections
        G.add_edges_from(self.all_connections)
        F.add_edges_from(self.all_connections)

        # Remove single nodes
        delete_nodes = []
        for node in F.nodes:
            if F.degree(node) == 0:
                delete_nodes.append(node)
        for node in delete_nodes:
            F.remove_node(node)
            G.remove_node(node)

        if draw:
            # ----------------PLOT----------------
            pos = nx.kamada_kawai_layout(G)
            plt.figure()
            nx.draw_networkx_nodes(F, pos)
            nx.draw_networkx_edges(F, pos)
            nx.draw_networkx_labels(F, pos, labels)
            # ------------------------------------

        if html:

            env = Environment(loader=PackageLoader('powertransmission', 'templates'),
                              autoescape=select_autoescape(['html', 'xml']))
            template = env.get_template('graph_visJS.html')

            nodes = []
            edges = []
            for label in list(labels.values()):
                nodes.append({'name': label, 'shape': 'circular'})

            for edge in G.edges:
                edge_dict = {}
                edge_dict['inode1'] = int(edge[0])-1
                edge_dict['inode2'] = int(edge[1])-1
                edges.append(edge_dict)

            options = {}
            s = template.render(
                name=self.stepfile,
                nodes=nodes,
                edges=edges,
                options=options)

            with open('graph_visJS.html', 'wb') as file:
                file.write(s.encode('utf-8'))

            webbrowser.open('file://' + os.path.realpath('graph_visJS.html'))

        return F

    def draw_graph(self):
        self.graph = self.create_graph()
        labels = {}
        for id_nb, function in self.functions.items():
            labels[id_nb] = str(id_nb)+' '+function.name
        pos = nx.kamada_kawai_layout(self.graph)
        plt.figure()
        nx.draw_networkx_nodes(self.graph, pos)
        nx.draw_networkx_edges(self.graph, pos)
        nx.draw_networkx_labels(self.graph, pos, labels)

    def step_subfunctions(self, subfunctions):
        subfunctions = subfunctions[0]
        parenthesis_count = 0
        subfunction_names = []
        subfunction_args = []
        subfunction_name = ""
        subfunction_arg = ""
        for char in subfunctions:
            
            if char == "(":
                parenthesis_count += 1
                if parenthesis_count == 1:
                    subfunction_names.append(subfunction_name)
                    subfunction_name = ""
                else:
                    subfunction_arg += char
                    
            elif char == ")":
                parenthesis_count -= 1
                if parenthesis_count == 0:
                    subfunction_args.append(subfunction_arg)
                    subfunction_arg = ""
                else:
                    subfunction_arg += char
                
            elif parenthesis_count == 0:
                subfunction_name += char
                
            else:
                subfunction_arg += char

        return [(subfunction_names[i], step_split_arguments(subfunction_args[i])) for i in range(len(subfunction_names))]
            
    def instanciate(self, instanciate_id, object_dict):
        """
        Returns None if the object was instanciate
        """
        
        # print('instanciate_id', instanciate_id)

        name = self.functions[instanciate_id].name
        arguments = self.functions[instanciate_id].arg[:]

        for i, arg in enumerate(arguments):
            if type(arg) == str and arg[0] == '#':
                arguments[i] = int(arg[1:])
            elif type(arg) == str and arg[0:2] == '(#':
                argument = []
                arg_id = ""
                for char in arg[1:-1]:
                    if char == ',':
                        argument.append(arg_id)
                        arg_id = ""
                        continue
                    
                    arg_id += char
                argument.append(arg_id)
                arguments[i] = list(argument)

        if name == 'VERTEX_POINT':
#            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        # elif name == 'LINE':
        #     pass

        elif name == 'ORIENTED_EDGE':
#            object_dict[instanciate_id] = object_dict[arguments[3]]
            volmdlr_object = object_dict[arguments[3]]

        elif name == 'FACE_OUTER_BOUND':
#            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'FACE_BOUND':
#            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'SURFACE_CURVE':
#            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]
        
        elif name == 'SEAM_CURVE':
#            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]
        # elif name == 'EDGE_CURVE':
        #     object_dict[instanciate_id] = object_dict[arguments[3]]

        elif name == 'VERTEX_LOOP' :
#            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'PCURVE' :
        #     # object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]
        
        elif name in step_to_volmdlr_primitive and hasattr(step_to_volmdlr_primitive[name], "from_step"):
            volmdlr_object = step_to_volmdlr_primitive[name].from_step(arguments, object_dict)

#            object_dict[instanciate_id] = volmdlr_object
#            if hasattr(volmdlr_object, "primitive"):
#                primitives.append(volmdlr_object.primitive)primitives
        
        else:
            print('name', name)
            print('arguments', arguments)
            raise NotImplementedError

        return volmdlr_object

    def to_shells3d(self, name):
        self.graph = self.create_graph()
        
        object_dict = {}

        self.graph.add_node("#0")
        for node in self.graph.nodes:
            if node != '#0' and (self.functions[node].name == "CLOSED_SHELL" or self.functions[node].name == "OPEN_SHELL"):
                self.graph.add_edge("#0", node)

        edges = list(nx.algorithms.traversal.breadth_first_search.bfs_edges(self.graph, "#0"))[::-1]
        for edge_nb, edge in enumerate(edges):
            instanciate_id = edge[1]
            volmdlr_object = self.instanciate(instanciate_id, object_dict)
            object_dict[instanciate_id] = volmdlr_object

        shells = []
        for node in list(self.graph.nodes):
            if node != '#0' and (self.functions[node].name == 'CLOSED_SHELL' or self.functions[node].name == "OPEN_SHELL"):
                shells.append(object_dict[node])
        return shells
    
    def to_scatter_volume_model(self, name):
        object_dict = {}
        points3d = []
        for stepfunction in self.functions.values():
            if stepfunction.name == 'CARTESIAN_POINT':
                # INSTANCIATION
                name = self.functions[stepfunction.id].name
                arguments = self.functions[stepfunction.id].arg[:]
                for i, arg in enumerate(arguments):
                    if type(arg) == str and arg[0] == '#':
                        arguments[i] = int(arg[1:])
                volmdlr_object = step_to_volmdlr_primitive[name].from_step(arguments, object_dict)
                points3d.append(volmdlr_object)
        return VolumeModel(points3d)

class VolumeModel(dc.DessiaObject):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes  = ['shells', 'bounding_box']
    _non_eq_attributes = ['name', 'shells', 'bounding_box', 'contours', 'faces']
    _non_hash_attributes = []
    """
    :param groups: A list of two element tuple. The first element is a string naming the group and the second element is a list of primitives of the group
    """
    def __init__(self, primitives, name=''):
        self.primitives = primitives
        self.name = name
        self.shells = []
        if self.primitives:
            self.shells = self._extract_shells()
        if self.shells:
            self.bounding_box = self._bounding_box()
        else : 
            self.bounding_box = BoundingBox(-1, 1, -1, 1, -1, 1)

    # def __hash__(self):
    #     return sum([hash(p) for p in self.primitives])

    def _extract_shells(self):
        shells = []
        for primitive in self.primitives:
            if isinstance(primitive, Shell3D):
                shells.append(primitive)
        return shells

    def Volume(self):
        volume=0
        for primitive in self.primitives:
            volume+=primitive.Volume()
        return volume

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_primitives = [primitive.Rotation(center, axis, angle, copy=True) for primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.Translation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_primitives = [primitive.Translation(offset, copy=True) for primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_primitives = [primitive.frame_mapping(frame, side, copy=True) for primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        new_primitives = [primitive.copy() for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def _bounding_box(self):
        bboxes = []
        points = []
        for primitive in self.primitives:
            if hasattr(primitive, 'bounding_box'):
                bboxes.append(primitive.bounding_box)
            else:
                if primitive.__class__.__name__ == 'Point3D':
                    points.append(primitive)
        if bboxes:
            xmin = min([box.xmin for box in bboxes])
            xmax = max([box.xmax for box in bboxes])
            ymin = min([box.ymin for box in bboxes])
            ymax = max([box.ymax for box in bboxes])
            zmin = min([box.zmin for box in bboxes])
            zmax = max([box.zmax for box in bboxes])
        elif points:
            xmin = min([p[0] for p in points])
            xmax = max([p[0] for p in points])
            ymin = min([p[1] for p in points])
            ymax = max([p[1] for p in points])
            zmin = min([p[2] for p in points])
            zmax = max([p[2] for p in points])
        else:
            # raise ValueError('Bounding box cant be determined')
            return BoundingBox(-1, 1, -1, 1, 1-1, 1)
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)
    
    def plot(self, ax=None, color=None):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        for i, shell in enumerate(self.shells):
            bbox = shell.bbox()
            bbox.plot(ax, color[i])

        return ax

    def MPLPlot(self):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', adjustable='box')
#        ax.set_aspect('equal')
        for primitive in self.primitives:
            primitive.MPLPlot(ax)
        ax.set_aspect('equal')
        ax.margins(0.1)
        return ax

    def FreeCADScript(self, fcstd_filepath,
                      freecad_lib_path='/usr/lib/freecad/lib',
                      export_types=('fcstd',),
                      save_to = '',
                      tolerance=0.0001):
        """
        Generate python a FreeCAD definition of model
        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :type fcstd_filename:str
        """
        fcstd_filepath = os.path.abspath(fcstd_filepath)
        fcstd_filepath = fcstd_filepath.replace('\\','\\\\')
        freecad_lib_path = freecad_lib_path.replace('\\','\\\\')

        s='# -*- coding: utf-8 -*-\n'
        if freecad_lib_path != '':
            s+="import sys\nsys.path.append('"+freecad_lib_path+"')\n"

        s+="import math\nimport FreeCAD as fc\nimport Part\n\ndoc=fc.newDocument('doc')\n\n"
        for ip, primitive in enumerate(self.primitives):
            if primitive.name == '':
                primitive_name = 'Primitive_{}'.format(ip)
            else:
                primitive_name = 'Primitive_{}_{}'.format(ip, primitive.name)
            s += "part = doc.addObject('App::Part','{}')\n".format(primitive_name)
            if hasattr(primitive, 'FreeCADExport'):
                sp = primitive.FreeCADExport(ip)
                if sp != '':
#                        s += (sp+'\n')
                    s += (sp)
                    s += 'shapeobj = doc.addObject("Part::Feature","{}")\n'.format(primitive_name)
                    if isinstance(primitive, BSplineCurve3D) \
                    or isinstance(primitive, BSplineSurface3D) \
                    or isinstance(primitive, Circle3D) \
                    or isinstance(primitive, LineSegment3D) \
                    or isinstance(primitive, Ellipse3D):
#                            print(primitive)
#                            s += 'S = Part.Shape([primitive{}])\n'.format(ip)
#                            s += 'shapeobj.Shape = S\n'
                        s += 'shapeobj.Shape = primitive{}.toShape()\n'.format(ip)
                    else:
                        s += "shapeobj.Shape = primitive{}\n".format(ip)
                    s += 'part.addObject(shapeobj)\n\n'.format(ip, primitive.name)
            # --------------------DEBUG-------------------
#                else:
#                    raise NotImplementedError
            # ---------------------------------------------

        s+='doc.recompute()\n'
        if 'fcstd' in export_types:
            s+="doc.saveAs('"+fcstd_filepath+".fcstd')\n\n"
        if 'stl' in export_types:
            s+="import Mesh\nMesh.export(doc.Objects,'{}.stl', tolerance={})\n".format(fcstd_filepath, tolerance)
        if 'step' in export_types:
            s+="Part.export(doc.Objects,'{}.step')\n".format(fcstd_filepath)


        if save_to != '':
            with open(os.path.abspath(save_to),'w') as file:
                file.write(s)
        return s

    def FreeCADExport(self,fcstd_filepath,
                      python_path='python',
                      freecad_lib_path='/usr/lib/freecad/lib',
                      export_types=('fcstd',),
                      tolerance=0.0001):
        """
        Export model to .fcstd FreeCAD standard

        :param python_path: path of python binded to freecad

            * on windows: something like C:\\\\Program Files\\\\FreeCAD X.XX\\\\bin\\\\python
            * on linux: python if installed by a dstribution package
        :param filepath: path of fcstd file (without extension)
        :param freecad_lib_path: FreeCAD.so lib path (/usr/lib/freecad/lib in general)
        :param tolerance: the tolerance of tesselation for mesh exports

        """
        fcstd_filepath=os.path.abspath(fcstd_filepath)
        s=self.FreeCADScript(fcstd_filepath,
                             freecad_lib_path = freecad_lib_path,
                             export_types = export_types,
                             tolerance=tolerance)
        with tempfile.NamedTemporaryFile(suffix=".py",delete=False) as f:
            f.write(bytes(s,'utf8'))

        arg = f.name
        output = subprocess.call([python_path, arg])

        f.close()
        os.remove(f.name)
        return output

    def babylon_script(self, use_cdn=True, debug=False):

        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))

        template = env.get_template('babylon.html')

        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        primitives_strings=[]
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_script'):
                primitives_strings.append(primitive.babylon_script())
                
        return template.render(name=self.name,
                               center=tuple(center),
                               length=2*max_length,
                               primitives_strings=primitives_strings,
                               use_cdn=use_cdn,
                               debug=debug)

    def babylonjs_from_script(self, page_name=None, use_cdn=True, debug=False):
        script = self.babylon_script(use_cdn=use_cdn, debug=debug)

        if page_name is None:
            with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as file:
                file.write(bytes(script,'utf8'))
            page_name = file.name
        else:
            page_name += '.html'
            with open(page_name,'w')  as file:
                file.write(script)

        webbrowser.open('file://' + os.path.realpath(page_name))

    def babylon_data(self):
        meshes = []
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])
                
        babylon_data = {'meshes': meshes,
                        'max_length': max_length,
                        'center': list(center)}
        return babylon_data

    @classmethod
    def babylonjs_from_babylon_data(cls, babylon_data, page_name=None, use_cdn=True, debug=False):
        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))

        template = env.get_template('babylon_unpacker.html')

        
        script = template.render(babylon_data=babylon_data,
                                 use_cdn=use_cdn,
                                 debug=debug
                                 )
        if page_name is None:
            with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as file:
                file.write(bytes(script,'utf8'))
            page_name = file.name
        else:
            page_name += '.html'
            with open(page_name,'w')  as file:
                file.write(script)

        webbrowser.open('file://' + os.path.realpath(page_name))

        
    def babylonjs(self, page_name=None, use_cdn=True, debug=False):
        # print('self.primitives', self.primitives)
        babylon_data = self.babylon_data()
        self.babylonjs_from_babylon_data(babylon_data, page_name = page_name,
                                         use_cdn = use_cdn, debug = debug)


        


class MovingVolumeModel(VolumeModel):
    def __init__(self, primitives, step_frames, name=''):
        VolumeModel.__init__(self, primitives=primitives, name=name)
        self.step_frames = step_frames
        
        if not self.is_consistent():
            raise dc.ConsistencyError
        
    def is_consistent(self):
        n_primitives = len(self.primitives)
        for frames in self.step_frames:
            if len(frames) != n_primitives:
                return False
        return True

    def step_volume_model(self, istep):
        primitives = []
        for primitive, frame in zip(self.primitives, self.step_frames[istep]):
            primitives.append(primitive.frame_mapping(frame, side='old', copy=True))
        return VolumeModel(primitives)

    def babylon_script(self, use_cdn=True, debug=False):

        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))


        template = env.get_template('babylon.html')

        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        primitives_strings=[]
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_script'):
                primitives_strings.append(primitive.babylon_script())

        positions = []
        orientations = []
        for step in self.step_frames:
            step_positions = []
            step_orientations = []

            for frame in step:
                step_positions.append(list(frame.origin))
                step_orientations.append([list(frame.u),
                                          list(frame.v),
                                          list(frame.w)])

            positions.append(step_positions)
            orientations.append(step_orientations)

        return template.render(name=self.name,
                               center=tuple(center),
                               length=2*max_length,
                               primitives_strings=primitives_strings,
                               positions=positions,
                               orientations=orientations,
                               use_cdn=use_cdn,
                               debug=debug)
    
    def babylon_data(self):
        meshes = []
        primitives_to_meshes = []
        for ip, primitive in enumerate(self.primitives):
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
                primitives_to_meshes.append(ip)
                
        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])
                
        steps = []
        for istep, frames in enumerate(self.step_frames):


            # step_positions = []
            # step_orientations = []
            step = {'time': istep}
            for iframe, frame in enumerate(frames):
                if iframe in primitives_to_meshes:
                    imesh = primitives_to_meshes.index(iframe)
                    step[imesh] = {}
                    step[imesh]['position'] = list(round(frame.origin, 6))
                    step[imesh]['orientations'] = [list(round(frame.u, 6)),
                                                    list(round(frame.v, 6)),
                                                    list(round(frame.w, 6))]

            steps.append(step)
        
        babylon_data = {'meshes': meshes,
                        'max_length': max_length,
                        'center': list(center),
                        'steps': steps}
        return babylon_data


class Routing:
    def __init__(self, point1, point2, volumemodel):
        self.points = [point1, point2]
        self.volumemodel = volumemodel

    def straight_line(self):
        """
        Returns 2 distances :
            - no collision distance
            - collision distance
        """
        line = LineSegment3D(self.points[0], self.points[1])

        intersection_points = []
        abscissea_list = []
        for shell in self.volumemodel.shells:
            for face in shell.faces:
                intersection_point, intersection_abscissea = face.linesegment_intersection(line, abscissea=True)
                if intersection_point is not None and intersection_abscissea != 0 and intersection_abscissea != 1:
                    not_in_abscissea_list = True
                    for abscissea in abscissea_list:
                        if math.isclose(abscissea, intersection_abscissea, abs_tol=1e-8):
                            not_in_abscissea_list = False
                    if not_in_abscissea_list:
                        intersection_points.append((intersection_point, intersection_abscissea))
                        abscissea_list.append(intersection_abscissea)

        if len(intersection_points)%2 != 0:
            raise NotImplementedError

        intersection_points = sorted(intersection_points, key=lambda abscissea: abscissea[1])
        all_points_abscissea = [(self.points[0], 0)] + intersection_points[:] + [(self.points[1], 1)]
        all_points = [p[0] for p in all_points_abscissea]

        no_collision_mesures = []
        collision_mesures = []
        i = 0
        for pt1, pt2 in zip(all_points[:-1], all_points[1:]):
            if i%2 == 0:
                no_collision_mesures.append(Measure3D(pt1, pt2, color=(0,0,1)))
            else:
                collision_mesures.append(Measure3D(pt1, pt2, color=(1,0,0)))
            i += 1

        return no_collision_mesures, collision_mesures

    def straight_line2(self):
        """
        Returns the distance of the line going around each shell's bbox encountered along the path
        """
        line = LineSegment3D(self.points[0], self.points[1])

        all_mesures_abscissea = []
        intersection_points = []
        for shell in self.volumemodel.shells:
            shell_intersection_points = []
            bbox = shell.bounding_box
            for face in shell.faces:
                intersection_point, intersection_abscissea = face.linesegment_intersection(line, abscissea=True)
                if intersection_point is not None and intersection_abscissea != 0 and intersection_abscissea != 1:
                    intersection_points.append((intersection_point, intersection_abscissea))
                    shell_intersection_points.append((intersection_point, intersection_abscissea))

            if len(shell_intersection_points) == 2:
                shell_intersection_points = sorted(shell_intersection_points, key=lambda abscissea: abscissea[1])
                abscissea1 = shell_intersection_points[0][1]
                abscissea2 = shell_intersection_points[1][1]
                shell_intersection_points = [p[0] for p in shell_intersection_points]
                around_bbox_mesures = bbox.distance_between_two_points_on_bbox(shell_intersection_points[0], shell_intersection_points[1])
                all_mesures_abscissea.append((around_bbox_mesures, abscissea1, abscissea2))
            elif len(shell_intersection_points) > 2 or len(shell_intersection_points) == 1:
                raise NotImplementedError

        intersection_points = sorted(intersection_points, key=lambda abscissea: abscissea[1])
        all_mesures_abscissea = sorted(all_mesures_abscissea, key=lambda abscissea: abscissea[1])
        all_points_abscissea = [(self.points[0], 0)] + intersection_points[:] + [(self.points[1], 1)]
        all_points = [p[0] for p in all_points_abscissea]

        no_collision_mesures = []
        i = 0
        for pt1, pt2 in zip(all_points[:-1], all_points[1:]):
            if i%2 == 0:
                no_collision_mesures.append(Measure3D(pt1, pt2, color=(0,0,1)))
            else:
                no_collision_mesures.extend(all_mesures_abscissea[i//2][0])
            i += 1

        return no_collision_mesures


class ViewIso:# TODO: rename this in IsoView
    def __init__(self, component, frame, size):
        self.component = component
        self.frame = frame
        self.size = size
        self.plot_datas = self.plot_data()

    def plot_data(self, detail=True):
        wide = min(self.size)/2
        plot_datas = []
        plot_datas.extend(self.component.plot_data(self.frame, detail=detail))
        plot_datas.extend(self.component.plot_data(Frame3D(self.frame.origin + Point3D((0, self.size[1]/2 + self.size[2]/2 + wide, 0)), self.frame.u, self.frame.w, self.frame.v), detail = detail))
        plot_datas.extend(self.component.plot_data(Frame3D(self.frame.origin + Point3D((self.size[0]/2 + self.size[2]/2 + wide, 0, 0)), self.frame.w, self.frame.v, self.frame.u), detail = detail))
        return plot_datas

    def plot(self):
        plot_data.plot(self.plot_datas)



step_to_volmdlr_primitive = {
        # GEOMETRICAL ENTITIES
        'CARTESIAN_POINT': Point3D,
        'DIRECTION': Vector3D,
        'VECTOR': Vector3D,

        'AXIS1_PLACEMENT': None,
        'AXIS2_PLACEMENT_2D': None, # ??????????????????
        'AXIS2_PLACEMENT_3D': Frame3D,

        'LINE': Line3D, #LineSegment3D,
        'CIRCLE': Circle3D,
        'ELLIPSE': Ellipse3D,
        'PARABOLA': None,
        'HYPERBOLA': None,
#        'PCURVE': None,
        'CURVE_REPLICA': None,
        'OFFSET_CURVE_3D': None,
        'TRIMMED_CURVE': None, # BSplineCurve3D cannot be trimmed on FreeCAD
        'B_SPLINE_CURVE': BSplineCurve3D,
        'B_SPLINE_CURVE_WITH_KNOTS': BSplineCurve3D,
        'BEZIER_CURVE': BSplineCurve3D,
        'RATIONAL_B_SPLINE_CURVE': BSplineCurve3D,
        'UNIFORM_CURVE': BSplineCurve3D,
        'QUASI_UNIFORM_CURVE': BSplineCurve3D,
        'SURFACE_CURVE': None, # TOPOLOGICAL EDGE
        'SEAM_CURVE': None, #LineSegment3D, # TOPOLOGICAL EDGE ############################
        'COMPOSITE_CURVE_SEGMENT': None, # TOPOLOGICAL EDGE
        'COMPOSITE_CURVE': Wire3D, # TOPOLOGICAL WIRE
        'COMPOSITE_CURVE_ON_SURFACE': Wire3D, # TOPOLOGICAL WIRE
        'BOUNDARY_CURVE': Wire3D, # TOPOLOGICAL WIRE

        'PLANE': Plane3D,
        'CYLINDRICAL_SURFACE': CylindricalSurface3D,
        'CONICAL_SURFACE': ConicalSurface3D,
        'SPHERICAL_SURFACE': SphericalSurface3D,
        'TOROIDAL_SURFACE': ToroidalSurface3D,
        'DEGENERATE_TOROIDAL_SURFACE': None,
        'B_SPLINE_SURFACE_WITH_KNOTS': BSplineSurface3D,
        'B_SPLINE_SURFACE': BSplineSurface3D,
        'BEZIER_SURFACE': BSplineSurface3D,
        'OFFSET_SURFACE': None,
        'SURFACE_REPLICA': None,
        'RATIONAL_B_SPLINE_SURFACE': BSplineSurface3D,
        'RECTANGULAR_TRIMMED_SURFACE': None,
        'SURFACE_OF_LINEAR_EXTRUSION': BSplineExtrusion, # CAN BE A BSplineSurface3D
        'SURFACE_OF_REVOLUTION': None,
        'UNIFORM_SURFACE': BSplineSurface3D,
        'QUASI_UNIFORM_SURFACE': BSplineSurface3D,
        'RECTANGULAR_COMPOSITE_SURFACE': PlaneFace3D, # TOPOLOGICAL FACES
        'CURVE_BOUNDED_SURFACE': PlaneFace3D, # TOPOLOGICAL FACE


        # TOPOLOGICAL ENTITIES
        'VERTEX_POINT': None,

        'EDGE_CURVE': Edge3D, # LineSegment3D, # TOPOLOGICAL EDGE
        'ORIENTED_EDGE': None, # TOPOLOGICAL EDGE
        # The one above can influence the direction with their last argument
        # TODO : maybe take them into consideration

        'FACE_BOUND': None, # TOPOLOGICAL WIRE
        'FACE_OUTER_BOUND': None, # TOPOLOGICAL WIRE
        # Both above can influence the direction with their last argument
        # TODO : maybe take them into consideration
        'EDGE_LOOP': Contour3D, # TOPOLOGICAL WIRE
        'POLY_LOOP': Contour3D, # TOPOLOGICAL WIRE
        'VERTEX_LOOP': None, # TOPOLOGICAL WIRE

        'ADVANCED_FACE': Face3D,
        'FACE_SURFACE': Face3D,

        'CLOSED_SHELL': Shell3D,
        'OPEN_SHELL': Shell3D,
#        'ORIENTED_CLOSED_SHELL': None,
        'CONNECTED_FACE_SET': Shell3D,

        }
