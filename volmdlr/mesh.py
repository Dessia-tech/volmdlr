#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Gasmi
"""

import matplotlib.pyplot as plt
import volmdlr.core_compiled 
import volmdlr
import volmdlr.core
from volmdlr.core_compiled import Matrix33

from itertools import combinations
import numpy as npy
import volmdlr.edges 
import volmdlr.wires 
import volmdlr.faces 
from volmdlr.core_compiled import Matrix33
import math
from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple,Dict
import matplotlib
import random
from itertools import product 
from matplotlib.colors import LinearSegmentedColormap
cdict = {'red':  [(0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0)],
         'green': [(0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)],
         'blue':  [(0.0, 1.0, 1.0),
                   (1.0, 0.0, 0.0)]}
blue_red = LinearSegmentedColormap('BLueRed', cdict)

class FlatElementError(Exception):
    pass

def find_duplicate_linear_element(linear_elements1, linear_elements2):
    duplicates = []
    for linear_element in linear_elements1:
        if linear_element in linear_elements2 and linear_element not in duplicates:
            duplicates.append(linear_element)
    return duplicates


class LinearElement(volmdlr.edges.LineSegment2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    def __init__(self, start:volmdlr.Point2D,end:volmdlr.Point2D, interior_normal:volmdlr.Vector2D,name=''):
        self.points=[start, end]
        self.interior_normal = interior_normal
        
        volmdlr.edges.LineSegment2D.__init__(self,start=start,end=end,name=name)

    
    def plot(self, ax=None, color='k', width=None, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        if width is None:
            width=1
        if plot_points:
            ax.plot([self.start.x, self.end.x], [self.start.y, self.end.y], color=color, marker='o', linewidth=width)
        else:
            ax.plot([self.start.x, self.end.x], [self.start.y, self.end.y], color=color, linewidth=width)
        return ax


class TriangularElement(volmdlr.wires.Triangle2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points):
        self.points = points
        self.linear_elements = self._to_linear_elements()
        self.form_functions = self._form_functions()
        
        self.center = (self.points[0]+self.points[1]+self.points[2])/3
        
        self.area = self._area()
        
        volmdlr.wires.Triangle2D.__init__(self,points=points, name='')
        
    def _to_linear_elements(self):
        vec1 = volmdlr.Vector2D(self.points[1].x - self.points[0].x,
                                self.points[1].y - self.points[0].y)
        vec2 = volmdlr.Vector2D(self.points[2].x - self.points[1].x,
                                self.points[2].y - self.points[1].y)
        vec3 = volmdlr.Vector2D(self.points[0].x - self.points[2].x,
                                self.points[0].y - self.points[2].y)
        normal1 = volmdlr.Vector2D(-vec1.y, vec1.x)
        normal2 = volmdlr.Vector2D(-vec2.y, vec2.x)
        normal3 = volmdlr.Vector2D(-vec3.y, vec3.x)
        normal1.normalize()
        normal2.normalize()
        normal3.normalize()
        if normal1.dot(vec2) < 0:
            normal1 = - normal1
        if normal2.dot(vec3) < 0:
            normal2 = - normal2
        if normal3.dot(vec1) < 0:
            normal3 = - normal3
        linear_element_1 = LinearElement(self.points[0], self.points[1],
                                         normal1)
        linear_element_2 = LinearElement(self.points[1], self.points[2],
                                         normal2)
        linear_element_3 = LinearElement(self.points[2], self.points[0],
                                         normal3)
        return [linear_element_1, linear_element_2, linear_element_3]
    
    def _form_functions(self):
        a = Matrix33(1, self.points[0].x, self.points[0].y,
                     1, self.points[1].x, self.points[1].y,
                     1, self.points[2].x, self.points[2].y)
        try:
            inv_a = a.inverse()
        except ValueError:
            self.plot()
            print('buggy element area', self.area)
            raise FlatElementError('form function bug')
        x1 = inv_a.vector_multiplication(volmdlr.X3D)
        x2 = inv_a.vector_multiplication(volmdlr.Y3D)
        x3 = inv_a.vector_multiplication(volmdlr.Z3D)
       
        return x1, x2, x3

    # def _quadratic_form_functions(self):
    #     a = [[1, self.points[0][0], self.points[0][1],self.points[0][0]**2,self.points[0][0]*self.points[0][1],self.points[0][1]**2],
    #           [1, self.points[1][0], self.points[1][1],self.points[1][0]**2,self.points[1][0]*self.points[1][1],self.points[1][1]**2],
    #           [1, self.points[2][0], self.points[2][1],self.points[2][0]**2,self.points[2][0]*self.points[2][1],self.points[2][1]**2],
    #           [1, self.points[3][0], self.points[3][1],self.points[3][0]**2,self.points[3][0]*self.points[3][1],self.points[3][1]**2],
    #           [1, self.points[4][0], self.points[4][1],self.points[4][0]**2,self.points[4][0]*self.points[4][1],self.points[4][1]**2],
    #           [1, self.points[5][0], self.points[5][1],self.points[5][0]**2,self.points[5][0]*self.points[5][1],self.points[5][1]**2]]

    #     try :
    #         inv_a = a.inverse()
    #     except ValueError:
    #         self.plot()
    #         print(self._area())
    #         raise FlatElementError('form function bug')
    #     x1 = inv_a.dot([1,0,0,0,0,0])
    #     x2 = inv_a.dot([1,0,0,0,0,0])
    #     x3 = inv_a.dot([1,0,0,0,0,0])
    #     x4=inv_a.dot([1,0,0,0,0,0])
        
    #     return x1, x2, x3


    def _area(self):
        u = self.points[1] - self.points[0]
        v = self.points[2] - self.points[0]
        return abs(u.cross(v)) / 2

    def point_belongs(self, point):
        polygon = volmdlr.wires.ClosedPolygon2D(self.points)
        point_belongs = polygon.point_belongs(point)
        return point_belongs

    def rotation(self, center, angle, copy=True):
        if copy:
            return TriangularElement([pt.rotation(center, angle, copy=True)
                                      for pt in self.points])
        else:
            for pt in self.points:
                pt.Rotation(center, angle, copy=False)
                
    def translation(self, offset, copy=True):
        if copy:
            return TriangularElement([pt.translation(offset, copy=True)
                                      for pt in self.points])
        else:
            for pt in self.points:
                pt.translation(offset, copy=False)
                
    def axial_symmetry(self, line, copy=True):
        p1, p2 = line.points
        symmetric_points = []
        for point in self.points:
            u = p2 - p1
            t = (point-p1).dot(u) / u.norm()**2
            projection = p1 + t * u
            symmetric_point = volmdlr.Point2D(*(2 * projection - point))
            symmetric_points.append(symmetric_point)
        if copy: 
            return TriangularElement(symmetric_points)
        else:
            for i, point in enumerate(self.points):
                point = symmetric_points[i]
   
    def plot(self, ax=None, color='k', width=None,
             plot_points=False, fill=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
            
        if fill:
            x = [p[0] for p in self.points]
            y = [p[1] for p in self.points]
            plt.fill(x, y, facecolor=color, edgecolor="k")
            return ax
        
        for p1, p2 in zip(self.points, self.points[1:]+[self.points[0]]):
            if width is None:
                width = 1
            if plot_points:
                ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
                        marker='o', linewidth=width)
            else:
                ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
                        linewidth=width)
        return ax


    def triangle_to_polygon(self):
        points = self.points
        return volmdlr.wires.ClosedPolygon2D(points)
       

class ElementsGroup(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, elements:List[TriangularElement], name:str):
        self.elements = elements
        self.name = name

        DessiaObject.__init__(self, name=name)

    def point_to_element(self, point):
        for element in self.elements:
            if element.point_belongs(point):
                return element
        return None

    def rotation(self, center, angle, copy=True):
        if copy:
            return Mesh([elem.rotation(center, angle, copy=True)
                         for elem in self.elements])
        else:
            for elem in self.elements:
                elem.rotation(center, angle, copy=False)
                
    def translation(self, offset, copy=True):
        if copy:
            return Mesh([elem.translation(offset, copy=True)
                         for elem in self.elements])
        else:
            for elem in self.elements:
                elem.translation(offset, copy=False)
                
    def plot(self, ax=None, color='k', fill=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        for element in self.elements:
            element.plot(ax=ax, color=color, fill=fill)
        return ax
        

class Mesh(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['node_to_index']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, elements_groups: List[ElementsGroup]):
        self.elements_groups = elements_groups
        self.nodes = self._set_nodes_number()
        self.node_to_index = {self.nodes[i]: i for i in range(len(self.nodes))}
        
        DessiaObject.__init__(self, name='')

    def __add__(self, other_mesh):
        new_nodes = self.nodes[:]
        new_nodes_index = {p: i for i, p in enumerate(self.points)}
        ip = len(new_nodes)
        for point in other_mesh.nodes:
            if point not in new_nodes_index:
                new_nodes_index[point] = ip
                ip += 1
                new_nodes.append(point)

        new_elements_groups = self.elements_groups[:]
        for i1, i2, i3 in other_mesh.elements_groups:
            p1 = other_mesh.nodes[i1]
            p2 = other_mesh.nodes[i2]
            p3 = other_mesh.nodes[i3]
            new_elements_groups.append((new_nodes_index[p1],
                                        new_nodes_index[p2],
                                        new_nodes_index[p3]))
        return self.__class__(new_elements_groups)

    def _set_nodes_number(self):
        nodes = set()
        for elements_group in self.elements_groups:
            for element in elements_group.elements:
                nodes.add(element.points[0])
                nodes.add(element.points[1])
                nodes.add(element.points[2])
        return tuple(nodes)

    def point_to_element(self, point):
        for element_group in self.elements_groups:
            element = element_group.point_to_element(point)
            if element is not None:
                return element
        return None
    
    def set_node_displacement_index(self):
        indexes = {}
        for node in self.nodes:
            indexes[node] = [2*self.node_to_index[node],
                             2*self.node_to_index[node]+1]
        return indexes
    
    def boundary_dict(self):
        boundary_dict = {}
        for elements_group1, elements_group2 in combinations(
                self.elements_groups, 2):
            linear_elements1 = []
            linear_elements2 = []
            for element in elements_group1.elements:
                linear_elements1.extend(element.linear_elements)
            for element in elements_group2.elements:
                linear_elements2.extend(element.linear_elements)
            duplicate_linear_elements = find_duplicate_linear_element(
                linear_elements1, linear_elements2)
            if duplicate_linear_elements:
                boundary_dict[(elements_group1,
                               elements_group2)] = duplicate_linear_elements
        return boundary_dict
   
    def plot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        for elements_group in self.elements_groups:
            elements_group.plot(ax=ax)
        return ax
    
    def plot_data(self, pos=0, quote=True, constructor=True, direction=1):
        plot_datas = []
        for element_group in self.elements_groups:
            for element in element_group.elements:
                c1 = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(
                    element.points[0], element.points[1])])
                c2 = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(
                    element.points[1], element.points[2])])
                c3 = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(
                    element.points[2], element.points[0])])
                plot_datas.append(c1.plot_data())
                plot_datas.append(c2.plot_data())
                plot_datas.append(c3.plot_data())
                # plot_datas.extend([c1, c2, c3])
        return plot_datas

    def plot_displaced_mesh(self,
                            node_displacement: Dict[volmdlr.Point2D,
                                                    List[float]],
                            ax=None, amplification=0.5):

        deformed_mesh = self.copy()
        nodes = deformed_mesh.nodes
  
        for node in nodes:
            for displaced_node in node_displacement:
                if node == displaced_node:
                    node.x += amplification*node_displacement[
                        displaced_node][0]
                    node.y += amplification*node_displacement[
                        displaced_node][1]
            
        ax = deformed_mesh.plot(ax=ax) 
        ax.set_aspect('equal')           
    
        return ax


class Mesher(DessiaObject):
    
    def __init__(self,interior_contours:List[volmdlr.wires.Contour2D],exterior_contours:List[volmdlr.wires.Contour2D],nodes_len:float):
        self.nodes_len=nodes_len
        self.interior_contours=interior_contours
        self.exterior_contours=exterior_contours
        

    
    def _is_convex(self,p1:volmdlr.Point2D, p2:volmdlr.Point2D, p3:volmdlr.Point2D):
        return self._triangle_sum(p1.x, p1.y, p2.x, p2.y, p3.x,p3.y) < 0
    
    def _is_clockwise(self,polygon:volmdlr.wires.ClosedPolygon2D):
        s = 0
        polygon_count = len(polygon.points)
        for i in range(polygon_count):
            point = polygon.points[i]
            point2 = polygon.points[(i + 1) % polygon_count]
            s += (point2.x - point.x) * (point2.y + point.y)
        return s > 0
    
    def _triangle_sum(self,x1, y1, x2, y2, x3, y3):
        return x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)
    
    def _contains_no_points(self,p1:volmdlr.Point2D,p2:volmdlr.Point2D,p3:volmdlr.Point2D, polygon:volmdlr.wires.ClosedPolygon2D):
       triangle=volmdlr.wires.Triangle2D([p1,p2,p3])
       for pn in polygon.points:
            if pn in [p1, p2, p3]:
                continue
            elif triangle.is_inside_triangle(pn):
                return False
       return True

    def _is_ear(self,p1:volmdlr.Point2D,p2:volmdlr.Point2D,p3:volmdlr.Point2D, polygon:volmdlr.wires.ClosedPolygon2D):
        triangle=volmdlr.wires.Triangle2D([p1,p2,p3])
        ear = self._contains_no_points(p1, p2, p3, polygon) and \
            self._is_convex(p1, p2, p3) and \
            triangle.area > 0
        return ear
    
    def earclip(self,polygon:volmdlr.wires.ClosedPolygon2D):
       
        possible_triangles=[]
        ear_vertex =[]
       
       
        if self._is_clockwise(polygon):
            polygon.points.reverse()
    
        point_count = len(polygon.points)
       
         
        for i in range(point_count)  :
                   
            prev_index = i - 1
            prev_point = polygon.points[prev_index]
          
            point = polygon.points[i]
            next_index = (i + 1) % point_count
            next_point = polygon.points[next_index]
    
            if self._is_ear(prev_point, point, next_point, polygon):
                ear_vertex.append(point)
                    
        while ear_vertex and point_count >=3  :
            ear = ear_vertex.pop(0)
            i = polygon.points.index(ear)
            prev_index = i - 1
            prev_point = polygon.points[prev_index]
            next_index = (i + 1) % point_count
            next_point = polygon.points[next_index]
            prev_prev_point = polygon.points[prev_index - 1]
            next_next_index = (i + 1) % point_count
            next_next_point = polygon.points[next_next_index]
            polygon.points.remove(ear)
            point_count -= 1
            
            p1=volmdlr.Point2D(prev_point[0], prev_point[1])
            p2=volmdlr.Point2D(ear[0], ear[1])
            p3= volmdlr.Point2D(next_point[0], next_point[1])
            triangle=volmdlr.wires.Triangle2D([p1,p2,p3])  
            possible_triangles.append(triangle)
             
            if point_count > 3:
                prev_prev_point = polygon.points[prev_index - 1]
                next_next_index = (i + 1) % point_count
                next_next_point = polygon.points[next_next_index]
    
                groups = [
                    (prev_prev_point, prev_point, next_point, polygon),
                    (prev_point, next_point, next_next_point, polygon),
                ]
                for group in groups:
                    p = group[1]
                    if self._is_ear(*group):
                        if p not in ear_vertex:
                            ear_vertex.append(p)
                    elif p in ear_vertex:
               
                         ear_vertex.remove(p)
                           
        return possible_triangles   
                 
                                 
        
    def basic_triangulation(self,polygon1:volmdlr.wires.ClosedPolygon2D,polygon2:volmdlr.wires.ClosedPolygon2D,
                           segment_to_nodes:Dict[volmdlr.edges.LineSegment2D,List[volmdlr.Point2D]]):
        """
        basic triangulation that simply links the points of the two polygons
        
        """
        
        triangles=[]
        
        
        for segment in polygon1.line_segments:
            index_0=polygon1.line_segments.index(segment)
            
            mid=segment.point_at_abscissa(segment.length()/2)
            
            d=[]
            for in_segment in polygon2.line_segments:
                  
      
                  l=in_segment.point_distance(mid)
                  d.append(l)
            index=d.index(min(d))
            near_segment=polygon2.line_segments[index]
                    
            p=segment_to_nodes[segment]
            q=segment_to_nodes[near_segment]
         
            u=len(p)
            v=len(q)
            if u==2 and v==2 : 
                new_triangle_1=volmdlr.wires.Triangle2D([p[0],p[1],q[0]])
                triangles.append(new_triangle_1)                          
                new_triangle_2=volmdlr.wires.Triangle2D([p[1],q[1],q[0]])
                triangles.append(new_triangle_2)                           
            if u>=v:
                for i in range(v-1):
                    new_triangle_1=volmdlr.wires.Triangle2D([p[i+1],p[i],q[i]])
                   
                    triangles.append(new_triangle_1)
                    new_triangle_2=volmdlr.wires.Triangle2D([q[i+1],p[i+1],q[i]])
                   
                    triangles.append(new_triangle_2)
                              
                for  i in range(v-1,u-1):
                                
                     new_triangle=volmdlr.wires.Triangle2D([p[i],q[v-1],p[i+1]])
                                     
                     triangles.append(new_triangle)  
            else :
                  for i in range(u-1):
                    new_triangle_1=volmdlr.wires.Triangle2D([q[i+1],q[i],p[i]])
                   
                    triangles.append(new_triangle_1)
                    new_triangle_2=volmdlr.wires.Triangle2D([p[i+1],q[i+1],p[i]])
                   
                    triangles.append(new_triangle_2)
                              
                  for  i in range(u-1,v-1):
                                
                     new_triangle=volmdlr.wires.Triangle2D([q[i],p[u-1],q[i+1]])
                                     
                     triangles.append(new_triangle)  
                             
        return triangles  
    def continuous_triangulation(self,polygon1,polygon2,continuity_points,segment_to_nodes): 
        """"
        If  contours/polygons are succesively in one another (offset method), this
        triangulation  ensures the continuity of the discretisation between
        each successive contour/polygon
        
        """
        all_triangles=[]
       
        for projections in continuity_points:
            if len(projections)>1:
                out_segment=volmdlr.edges.LineSegment2D(projections[0],projections[-1])
               
                mid=out_segment.point_at_abscissa(out_segment.length()/2)
                
                d=[]
                for in_segment in polygon2.line_segments:
                      
          
                      l=in_segment.point_distance(mid)
                      d.append(l)
                      
                index=d.index(min(d))
                near_segment=polygon2.line_segments[index]
           
                
                u=len(segment_to_nodes[near_segment])
                v=len(projections)
                
                if u>=v and u>2:
                  
                    for j in range(v-1):
                        
                        new_triangle=volmdlr.wires.Triangle2D([projections[j], segment_to_nodes[near_segment][j+1],
                                            segment_to_nodes[near_segment][j]])   
        
                        all_triangles.append(new_triangle)
                        
                        new_triangle_0=volmdlr.wires.Triangle2D([projections[j+1],projections[j],
                                            segment_to_nodes[near_segment][j+1]])   
                    
                        all_triangles.append(new_triangle_0)
                        
                    for j in  range(v-1,u-1):
                        new_triangle=volmdlr.wires.Triangle2D([segment_to_nodes[near_segment][j], segment_to_nodes[near_segment][j+1],
                                       projections[v-1]])   
                        all_triangles.append(new_triangle) 
               
                if u<v :
                    for j in range(u-1):
                      
                       new_triangle=volmdlr.wires.Triangle2D([segment_to_nodes[near_segment][j],projections[j+1],
                                          projections[j]])   
                                    
                       all_triangles.append(new_triangle)
                       new_triangle_0=volmdlr.wires.Triangle2D([segment_to_nodes[near_segment][j+1],segment_to_nodes[near_segment][j],
                                          projections[j+1]])   
    
                       all_triangles.append(new_triangle_0)
                       
                    for j in  range(u-1,v-1):
                       new_triangle=volmdlr.wires.Triangle2D([projections[j],  projections[j+1],
                                      segment_to_nodes[near_segment][u-1]])   
                       all_triangles.append(new_triangle)
                       
                if u==2 and v==2:
                    if segment_to_nodes[near_segment][1] != segment_to_nodes[near_segment][0] :
                        new_triangle_1=volmdlr.wires.Triangle2D([segment_to_nodes[near_segment][0],segment_to_nodes[near_segment][1],
                                           projections[1]]) 
                        all_triangles.append(new_triangle_1)
                      
                    new_triangle_2=volmdlr.wires.Triangle2D([segment_to_nodes[near_segment][0],projections[0],
                                        projections[1]]) 
                   
                    all_triangles.append(new_triangle_2)
                
                
        return all_triangles 
                           
    def generate_offsets(self,contour,nb_offset):
        # ax=plt.subplot()
        # contour.plot(ax=ax)
        contour_offsets=[]
        all_offsets=[]
        offset_values=[]
        p=1.2
        k=1
        xmin,xmax,ymin,ymax=contour.bounding_points()
        offset_len=max(abs(xmax-xmin),abs(ymax-ymin))
      
        
        while k<10:
           
            repair = False
            contour_offsets=[]
            contour_offsets+=[contour]
            while repair is False :
                
                new_contour=contour.offset(k*offset_len/(100))
                
                if not new_contour.self_intersections() :
                   
                    if new_contour.to_polygon(self.nodes_len).area() > contour.to_polygon(self.nodes_len).area():
                            good_offset=contour.offset(-k*offset_len/100)
                            
                            if not good_offset.self_intersections():
                               contour_offsets.append(good_offset)
                                 
                            else :
                                contour_offsets.append(good_offset.select_reapaired_contour([]))
                            
                                offset_values.append(p*offset_len/(10*k))
                                
                    else :
                
                             offset_values.append(-p*offset_len/(10*k))
                             contour_offsets.append(new_contour)
                           
                 
                else:                
                   
                    contour_offsets.append(new_contour.select_reapaired_contour([]))
            
                    
                if nb_offset==2: 
                        
                    xmin_2,xmax_2,ymin_2,ymax_2=contour_offsets[-1].bounding_points()
                    # offset_len2=min(abs(xmax_2-xmin_2),abs(ymax_2-ymin_2))
                    offset_len2=offset_len
                    offset2=contour_offsets[-1].offset(k*offset_len2/100)
                    if not offset2.self_intersections():
                        
                  
                        
                        if offset2.to_polygon(self.nodes_len).area() > contour_offsets[-1].to_polygon(self.nodes_len).area() :
                            good_offset_2=contour_offsets[-1].offset(-k*offset_len2/100)
                          
                            if not good_offset_2.self_intersections():
                                
                                contour_offsets.append(good_offset_2) 
                                
                                offset_values.append(p/(10*k))
                               
                           
                        else :  
                            
                            contour_offsets.append(offset2)      
                            
                    # else:                                                
                
                    repair=True
                    if len(contour_offsets)>2: 
                        
                        all_offsets.append(contour_offsets)
                else :
                   repair=True
                   if len(contour_offsets)>1:
                       
                       all_offsets.append(contour_offsets)  
                   
                 
            k=k+0.5
            if p<10:
                p=p+0.51
            
        return all_offsets,offset_values
    
    def continuous_earclip(self,contour,continuity_points):
        
        polygon=contour.to_polygon(self.nodes_len)
        all_points=[]
        all_triangles=[]
        for points in continuity_points:
            for point in points :
               if point not in all_points:
                    all_points.append(point)
        
        last_polygon=volmdlr.wires.ClosedPolygon2D(all_points)  
        if polygon.is_convex():
          
            return last_polygon.delaunay_triangulation()
        else :
            ear=self.earclip(last_polygon)
            for triangle in ear :
               
                if triangle.area>10e-9:
                    # if triangle.aspect_ratio()<20:
                        all_triangles.append(triangle)
            return all_triangles
    
    def complete_mesh(self,contour):
        triangles=[]
        segment_to_nodes={}
        polygon=contour.to_polygon(self.nodes_len)
        last_points=[]
        for segment in polygon.line_segments:
            segment_to_nodes[segment]=segment.discretise(self.nodes_len)
            for point in segment_to_nodes[segment]:
                if point not in last_points:
                    last_points.append(point)
        last_polygon=volmdlr.wires.ClosedPolygon2D(last_points)
        if polygon.is_convex():
            
            triangles+=last_polygon.delaunay_triangulation()
        else :
            ear=self.earclip(last_polygon)
        
            for triangle in ear : 
             
                if triangle.area>10e-9:
                    triangles.append(triangle)
        return triangles
    
    
    def alternative_triangulation(self,contour:volmdlr.wires.Contour2D,
      interior_contour:volmdlr.wires.Contour2D,empty:bool,far:bool,nb_offset):
       
        all_aspect_ratios=[]
        segment_to_nodes={}
        all_offsets,offset_values=self.generate_offsets(contour,nb_offset)
        all_meshes=[]
        good_meshes=[]
        good_offsets=[]
        
        if interior_contour is not None :
            interior_polygon=interior_contour.to_polygon(self.nodes_len)
            segment_to_nodes.update(interior_polygon.discretise(self.nodes_len))
                    
        for contour_offsets in all_offsets:
            
            polygon_offsets=[contour.to_polygon(self.nodes_len) for contour in contour_offsets]
            l=len(polygon_offsets)
            all_length=[len(polygon_offsets[k].line_segments)==len(polygon_offsets[k+1].line_segments) for k in range(l-1)]
            offset_triangles=[]  
            
            for polygon in polygon_offsets:
                segment_to_nodes.update(polygon.discretise(self.nodes_len))
            
            if  all(all_length):
                 
                  offset_triangles+=self.basic_triangulation(polygon_offsets[0],polygon_offsets[1],
                                                segment_to_nodes)
            else :
                 
                  offset_triangles+=self.mesh_in_between(polygon_offsets[1],polygon_offsets[0],empty)
                       
                    
            if l >2:
                
                for offset_1,offset_2,offset_0 in zip(polygon_offsets[1:],polygon_offsets[2:]
                                                      ,polygon_offsets[0:]):    
                   
                      
                    if all(all_length):
                       
                        offset_triangles+=self.basic_triangulation(offset_1,
                                                             offset_2,segment_to_nodes)
                    else :
                          
                        continuity_points=self.projection_points(offset_1,offset_0,segment_to_nodes)
                        offset_triangles+=self.continuous_triangulation(offset_1
                                                            ,offset_2,continuity_points,segment_to_nodes)
            
            if not empty :                           
                if all(all_length):
                     
                    contour_last_offset=contour_offsets[-1]
                    offset_triangles+=self.complete_mesh(contour_last_offset)  
                else :
                    # if cut : 
                    #     contour_last_offset=contour_offsets[-1]
                    #     xmin, xmax, ymin, ymax = contour_last_offset.bounding_points()
            
                    #     xi = xmin +  (xmax - xmin) / 2
               
                    #     cut_line = volmdlr.edges.Line2D(volmdlr.Point2D(xi, 0),
                    #                             volmdlr.Point2D(xi, 1))
                    #     cutted = contour_last_offset.cut_by_line(cut_line)
                    #     for c in cutted :
                    #         cutted_polygon=c.to_polygon(self.nodes_len)
                    #         p=self.projection_points(cutted_polygon,polygon_offsets[l-2],segment_to_nodes)  
                    #         offset_triangles+=self.continuous_earclip(c,p)
                    # else :   
                    last_continuity_points=self.projection_points(polygon_offsets[-1],polygon_offsets[l-2],segment_to_nodes)
                    offset_triangles+=self.continuous_earclip(polygon_offsets[-1],last_continuity_points)
            else :
                   if far :
                       interior_polygon=interior_contour.to_polygon(self.nodes_len)
 
                       offset_triangles+=self.mesh_in_between(interior_polygon,polygon_offsets[-1],True)     
       
            all_meshes.append(offset_triangles)   
            
            
        
        for mesh in all_meshes:
            if self.triangulation_max_aspect_ratio(mesh)!=0:
                
                all_aspect_ratios.append(self.triangulation_min_aspect_ratio(mesh))
                p=all_meshes.index(mesh)
                good_meshes.append(mesh)
                good_offsets.append(all_offsets[p])
                
        index=all_aspect_ratios.index(min(all_aspect_ratios))
        
     
       
        return good_meshes[index]
    
    def projection_points(self,in_polygon,out_polygon,segment_to_nodes):
        
        """
        returns the projection of the space of points 
        of the outer polygon on the inner polygon
        
        """

    
        projection_points=[]
        closest_segment={}
        out_point_image={}
      
        for segment in out_polygon.line_segments:
            for point in segment_to_nodes[segment]:
                out_point_image[point]=[]
        for segment in out_polygon.line_segments:
            projection_points.append([]) 
                             
        for out_segment in out_polygon.line_segments:
            
            index_0=out_polygon.line_segments.index(out_segment)
            
            mid=out_segment.point_at_abscissa(out_segment.length()/2)
            
            d=[]
            for in_segment in in_polygon.line_segments:
                  
      
                  l=in_segment.point_distance(mid)
                  d.append(l)
            index=d.index(min(d))
            near_segment=in_polygon.line_segments[index]
            
            closest_segment[out_segment]=near_segment
            for point in segment_to_nodes[out_segment]:

                i=[]
                d_1=[]
                index_point=segment_to_nodes[out_segment].index(point)
                projection, _=near_segment.point_projection(point)
                
                for in_segment in in_polygon.line_segments:
                    
                    if point.point_distance(in_segment.start) < point.point_distance(in_segment.end):
                        
                        d_1.append(point.point_distance(in_segment.start))
                        
                        i.append(0)
                    elif  point.point_distance(in_segment.start) > point.point_distance(in_segment.end):
                        d_1.append(point.point_distance(in_segment.end))
                        i.append(1)
                    elif point.point_distance(in_segment.start) == point.point_distance(in_segment.end):
                         
                          d_1.append(point.point_distance(in_segment.start))
                          i.append(0)
                index_1=d_1.index(min(d_1))
                new_proj=in_polygon.line_segments[index_1].points[i[index_1]]
                
                if len(segment_to_nodes[out_segment])==2  : 
                    
                    if out_point_image[point]==[]:
                        out_point_image[point].append([index_0,new_proj])
                        if new_proj not in projection_points[index_0]:
                          projection_points[index_0].insert(index_point,new_proj)
                    else :
                          if out_point_image[point][0][1] != new_proj :
                            
                 
                        
                            if new_proj not in projection_points[index_0]:
                                    projection_points[index_0].insert(index_point,new_proj)
                            # projection_points[index_0][0].plot(ax=ax,color='r')  
                            # projection_points[out_point_image[point][0][0]].remove(out_point_image[point][0][1])
                            projection_points[out_point_image[point][0][0]].append(new_proj)
                            
                           
                          else :
                              if new_proj not in projection_points[index_0]:
                                    projection_points[index_0].insert(index_point,new_proj)
                                                                                      
                else:                                               
                    if projection is None :
                            
                           
                        if out_point_image[point]==[]:
                          
                            out_point_image[point].append([index_0,new_proj])
                              
                            if new_proj not in projection_points[index_0]:
                               
                                projection_points[index_0].insert(index_point,new_proj)
                            
                        else :
                           
                             # out_point_image[point][0][1].plot(ax=ax,color='b')
                             
                             if out_point_image[point][0][1] != new_proj :
                                
                                  # projection_points[out_point_image[point][0][0]].appebd(out_point_image[point][0][1])
                                  # out_point_image[point][0][1].plot(ax=ax,color='r')
                                  # projection_points[index_0].insert(index_point,out_point_image[point][0][1])
                                   projection_points[out_point_image[point][0][0]].append(new_proj)
                                   if new_proj not in projection_points[index_0]:
                                         projection_points[index_0].insert(index_point,new_proj)
                                 
                             else :
                                 if new_proj not in projection_points[index_0]:
                                    projection_points[index_0].insert(index_point,new_proj)
                    else :
                        x=0
                        
                        if out_point_image[point]==[]:
                              if out_segment.length()/near_segment.length()>2.5:  
                                  i=0
                                  
                                  if projection.point_distance(near_segment.start) > projection.point_distance(near_segment.end) :
                                     i+=1
                                  if projection != near_segment.points[i]:
                                      l = near_segment.length()/projection.point_distance(near_segment.points[i])
                                      
                                      if out_segment.length()/projection.point_distance(near_segment.points[i])<3:
                                         x=10
                                      else :
                                          x=2
                                      if l <= x :
                                       
                                          out_point_image[point].append([index_0,projection]) 
                                          projection_points[index_0].insert(index_point,projection)
                                      else :
                                           
                                            if new_proj not in projection_points[index_0]:     
                                                projection_points[index_0].insert(index_point,new_proj)
                                  else :
                                      out_point_image[point].append([index_0,projection]) 
                                      projection_points[index_0].insert(index_point,projection) 
                              else :
                                      out_point_image[point].append([index_0,projection]) 
                                      projection_points[index_0].insert(index_point,projection) 
                                
                        else :
                            if out_point_image[point][0][1] != projection :
                            
                                if new_proj not in projection_points[index_0]:
                                        projection_points[index_0].insert(index_point,projection)
                                # projection_points[index_0][0].plot(ax=ax,color='r')  
                                # projection_points[out_point_image[point][0][0]].remove(out_point_image[point][0][1])
                                projection_points[out_point_image[point][0][0]].append(projection_points[index_0][0])
                                
                               
                            else :
                                  if projection not in projection_points[index_0]:
                                        projection_points[index_0].insert(index_point,projection)
                                                                                          
                                   
        return projection_points  
    
    def in_between_triangulation(self,out_polygon,projection_points,segment_to_nodes):
        all_triangles=[]
        """
        triangulates between an outer and inner polygon by using the projections points
        of the outer polygon on the inner polygon, and  of points of the outer polygon
        """
     
       
    
        for out_segment in out_polygon.line_segments:
            
            index_0=out_polygon.line_segments.index(out_segment)
            v=len(segment_to_nodes[out_segment])
            u=len(projection_points[index_0])
            
   
            if u>=v and u>2:
             
                    
                       
                for j in range(v-1):
                    
                    new_triangle=volmdlr.wires.Triangle2D([projection_points[index_0][j],projection_points[index_0][j+1],
                                        segment_to_nodes[out_segment][j]])   
    
                    all_triangles.append(new_triangle)
                    
                    new_triangle_0=volmdlr.wires.Triangle2D([projection_points[index_0][j+1],segment_to_nodes[out_segment][j],
                                        segment_to_nodes[out_segment][j+1]])   
                
                    all_triangles.append(new_triangle_0)
                    
                for j in  range(v-1,u-1):
                    # projection_points[index_0][j+1].plot(ax=ax,color='b')
                    new_triangle=volmdlr.wires.Triangle2D([projection_points[index_0][j],projection_points[index_0][j+1],
                                    segment_to_nodes[out_segment][v-1]])  
                    # new_triangle.plot(ax=ax,color='r')
                    all_triangles.append(new_triangle) 
           
            if u<v :
                for j in range(u-1):
                  
                   new_triangle=volmdlr.wires.Triangle2D([projection_points[index_0][j],segment_to_nodes[out_segment][j+1],
                                      segment_to_nodes[out_segment][j]])   
                                
                   all_triangles.append(new_triangle)
                   new_triangle_0=volmdlr.wires.Triangle2D([projection_points[index_0][j+1],projection_points[index_0][j],
                                       segment_to_nodes[out_segment][j+1]])   

                   all_triangles.append(new_triangle_0)
                   
                for j in  range(u-1,v-1):
                   new_triangle=volmdlr.wires.Triangle2D([segment_to_nodes[out_segment][j],  segment_to_nodes[out_segment][j+1],
                                  projection_points[index_0][u-1]])   
                   all_triangles.append(new_triangle)
                   
            if u==2 and v==2:
                
                if projection_points[index_0][1] != projection_points[index_0][0] :
                    new_triangle_1=volmdlr.wires.Triangle2D([projection_points[index_0][0],projection_points[index_0][1],
                                        segment_to_nodes[out_segment][1]]) 
                    all_triangles.append(new_triangle_1)
                  
                new_triangle_2=volmdlr.wires.Triangle2D([projection_points[index_0][0],segment_to_nodes[out_segment][0],
                                    segment_to_nodes[out_segment][1]]) 
               
                all_triangles.append(new_triangle_2)
                 
        return all_triangles    
    
    def mesh_in_between(self,in_polygon:volmdlr.wires.ClosedPolygon2D,
                        out_polygon:volmdlr.wires.ClosedPolygon2D,empty:bool):
    
        segment_to_nodes={}
        all_triangles=[]

        segment_to_nodes.update(out_polygon.discretise(self.nodes_len))
        segment_to_nodes.update(in_polygon.discretise(0))

        projection_points=self.projection_points(in_polygon,out_polygon,
                                                 segment_to_nodes)  
        
        all_triangles+=self.in_between_triangulation(out_polygon,projection_points,
                                                     segment_to_nodes)
        return all_triangles
  
    
    
    def polygon_to_triangles(self,polygons:List[volmdlr.wires.ClosedPolygon2D]):
        triangles=[]
        for polygon in polygons:
            triangles.append(volmdlr.wires.Triangle2D(polygon.points))
        return triangles

    
    def triangulation_max_aspect_ratio(self,triangles:List[volmdlr.wires.Triangle2D]):
        all_aspect_ratios=[]
        for triangle in triangles:
            if triangle.area<10E-9:
                return 0
            else :
                all_aspect_ratios.append(triangle.aspect_ratio())
        index=all_aspect_ratios.index(max(all_aspect_ratios))
        
        return all_aspect_ratios[index]
    
    
    def triangulation_min_aspect_ratio(self,triangles:List[volmdlr.wires.Triangle2D]):
        all_aspect_ratios=[]
        for triangle in triangles:
            if triangle.area<10E-9:
                return 0
            else :
                all_aspect_ratios.append(triangle.aspect_ratio())
        index=all_aspect_ratios.index(min(all_aspect_ratios))
        
        return all_aspect_ratios[index]
    
    def exterior_contours_triangulation(self,polygon,contour,min_aspect_ratio):
        
        triangles=[]
        possible_triangles=self.earclip(polygon)
        
        """For one polygon,earclip created a LineSegment2D as an instance of Triangle2D
        for unknown reasons,hence the area test below """
        
        for triangle in possible_triangles:
            if triangle.area<10E-9:
                possible_triangles.remove(triangle)
                
                
        if min_aspect_ratio == None:
          triangles+=possible_triangles
         
        else: 
            Next=True
            k=0
            while Next is True and k < len(possible_triangles):
             
                triangle=possible_triangles[k]    
               
                if triangle.aspect_ratio() > min_aspect_ratio:
                  Next=False
                  
                  triangles+=self.alternative_triangulation(contour,
                                                            None,False,False,1)
                  
                else :
                
                    k=k+1
                      
            if Next==True :
                triangles+=possible_triangles
                
                
        return triangles   
             
    def mesh_triangulation(self,triangulation):
        
         plot_aspect_ratio_triangles=[]
         mesh_triangles=[]
         segment_to_nodes={}
         all_segments=set()
         all_aspect_ratios={}
         
         for triangle in triangulation:
            
            all_segments= all_segments.union(triangle.line_segments)
            
         for segment in all_segments:
            
            segment_to_nodes[segment]=segment.discretise(self.nodes_len)
            
         for triangle in triangulation :
            
            meshing=triangle.mesh_triangle(segment_to_nodes,self.nodes_len)
            plot_aspect_ratio_triangles+=meshing[0]
            mesh_triangles+=meshing[0]
            all_aspect_ratios.update(meshing[1])
            
      
         return mesh_triangles,all_aspect_ratios,plot_aspect_ratio_triangles

                    
    def generate_mesh(self,min_aspect_ratio:float):
        
        split=False
        ax=plt.subplot()
        interior_polygons=[]
        exterior_polygons=[]
        mesh_triangle_elements=[]
        triangulation=[]
        
        if len(self.interior_contours)>1 : 
            split=True
        
        if self.interior_contours :
            
            for contour in self.interior_contours:
                   
                interior_polygons.append(contour.to_polygon(self.nodes_len))  
       
        for contour in self.exterior_contours:
                 
              exterior_polygons.append(contour.to_polygon(self.nodes_len))   

        if split is True :
         
            Surface=volmdlr.faces.Surface2D(self.exterior_contours[-1],self.interior_contours)
            split_contours = Surface.split_at_centers()
            
            for s in split_contours:
 
                triangulation+=self.alternative_triangulation(s,None,True,False,1)
            
        else :
                      
            """
            The convention adopted is the following : The last polygon of 
            interior_polygons(which stays empty) is inside the last polygon of 
            exterior_polygon 
            
            """ 
         
            if self.interior_contours:
                if exterior_polygons[-1].polygon_distance(interior_polygons[-1]) > exterior_polygons[-1].max_length()/3 :
                    
                    triangulation+=self.alternative_triangulation(self.exterior_contours[-1],
                                                  self.interior_contours[-1],True,True,2)
                else :                
                    
                    triangulation+=self.mesh_in_between(interior_polygons[-1],
                                                    exterior_polygons[-1],True)
                    
                if len(exterior_polygons)>1:
                    for polygon in exterior_polygons:
                        if polygon != exterior_polygons[-1]:
                            p=exterior_polygons.index(polygon)
                            triangulation+=self.exterior_contours_triangulation(polygon,
                                    self.exterior_contours[p],min_aspect_ratio)
                          
            else :
                for polygon in exterior_polygons:
                    p=exterior_polygons.index(polygon) 
                    triangulation+=self.exterior_contours_triangulation(polygon,
                            self.exterior_contours[p],min_aspect_ratio)
                  
        mesh_triangles,all_aspect_ratios,plot_aspect_ratio_triangles=self.mesh_triangulation(triangulation)
      
        for triangle in mesh_triangles:
            # triangle.plot(ax=ax)
            triangle.plot_triangle_points(ax=ax)      
            triangular_element=TriangularElement(triangle.points)
            mesh_triangle_elements.append(triangular_element)
            
        # self.plot_aspect_ratio(plot_aspect_ratio_triangles,all_aspect_ratios,ax)
        ax.set_aspect('equal')
        return mesh_triangle_elements
    
    def plot_aspect_ratio(self,all_triangles:List[volmdlr.wires.Triangle2D],
                          all_aspect_ratios:Dict[volmdlr.wires.Triangle2D,float],
                          ax,min_aspect_ratio=None,max_aspect_ratio=None):
        """
        Plots the mesh with colored triangles representing the 
        the gradient of values of the aspect ratio of each triangle. 
        """
        color_map = ((0,0,1), (1,0,0))
        if ax is None :
            fig, ax = plt.subplots()
        else :
            fig = plt.gcf()
        A= [a for a in list(all_aspect_ratios.values())]
        
        if max_aspect_ratio is None:
           max_aspect = max(A)
           
        if min_aspect_ratio is None:
           min_aspect = min(A)
     
        aspect_ratio_to_color = {}
        for a in A:
            if a > max_aspect:
                x = 1
            else:
                x = (a - min_aspect) / (max_aspect - min_aspect)
            color = (color_map[0][0]-(color_map[0][0]-color_map[1][0])*x, 
                     color_map[0][1]-(color_map[0][1]-color_map[1][1])*x,
                     color_map[0][2]-(color_map[0][2]-color_map[1][2])*x)
            aspect_ratio_to_color[a] = color
        
        for triangle in all_triangles:
            triangle.plot(ax=ax, color=aspect_ratio_to_color[all_aspect_ratios[triangle]], 
                          fill=True) 
            
        norm = matplotlib.colors.Normalize(vmin=min_aspect, vmax=max_aspect)
        sm = plt.cm.ScalarMappable(cmap=blue_red, norm=norm) 
        sm.set_array([])
        cbar = fig.colorbar(sm, ticks=npy.linspace(min_aspect,max_aspect, 10))
        cbar.set_label('Aspect Ratio')
        
        return ax    
      
    
        
               
           
           
       
       
            
        
        
    
  
        
                    
            
       
        
    
    
    
    