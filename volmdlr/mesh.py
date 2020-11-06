#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:05:56 2020

@author: ringhausen
"""

import matplotlib.pyplot as plt
import volmdlr.core as vm
import volmdlr.core_compiled as vmc
from itertools import combinations
import numpy as npy
from volmdlr.core_compiled import Matrix33
import math
from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple,Dict
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import random
from scipy.spatial import Delaunay
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


class LinearElement(vm.LineSegment2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    def __init__(self, point1:vm.Point2D,point2:vm.Point2D, interior_normal:vm.Vector2D,name=''):
        self.points=[point1, point2]
        self.interior_normal = interior_normal
        
        vm.LineSegment2D.__init__(self,point1=point1,point2=point2,name=name)
        
    def __hash__(self):
        return self.point1.__hash__() + self.point2.__hash__()
        
    def __eq__(self, other_linear_element):
        if self.__class__ != other_linear_element.__class__:
            return False
        return (self.point1 == other_linear_element.point1 and self.point2 == other_linear_element.point2) \
            or (self.point1 == other_linear_element.point2 and self.point2 == other_linear_element.point1)
            

    
    # def closest_middle(self,points:List[vm.Point2D]):
        
    #     mid=vm.Point2D([(self.point1[0]+self.point2[0])/2,(self.point1[1]+self.point2[1])/2])
    #     d=[]
    #     u=[]
    #     for point in points :
    #         if d!=0:
    #             d.append(math.sqrt((mid[0]-point[0])**2+(mid[1]-point[1])**2))
    #             u.append(point)
        
    #     index=d.index((min(d)))
    #     return index
    
    
        
    # def common_edge(self,nodes_1:List[vm.Point2D],nodes_2:List[vm.Point2D]):
    #     common_edge=[]
    #     for point_1 in nodes_1:
    #         for point_2 in nodes_2:
    #             if point_1==point_2:
    #                 common_edge.append(point_1)
    #     if len(common_edge)==1:
    #         return True
    #     return False
    
    
    def plot(self, ax=None, color='k', width=None, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        if width is None:
            width=1
        if plot_points:
            ax.plot([self.point1[0], self.point2[0]], [self.point1[1], self.point2[1]], color=color, marker='o', linewidth=width)
        else:
            ax.plot([self.point1[0], self.point2[0]], [self.point1[1], self.point2[1]], color=color, linewidth=width)
        return ax

class TriangularElement(vm.Triangle2D):
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
        
        vm.Triangle2D.__init__(self,points=points, name='')
        
    def _to_linear_elements(self):
        vec1 = vm.Vector2D(self.points[1] - self.points[0])
        vec2 = vm.Vector2D(self.points[2] - self.points[1])
        vec3 = vm.Vector2D(self.points[0] - self.points[2])
        normal1 = vm.Vector2D([-vec1[1], vec1[0]])
        normal2 = vm.Vector2D([-vec2[1], vec2[0]])
        normal3 = vm.Vector2D([-vec3[1], vec3[0]])
        normal1.Normalize()
        normal2.Normalize()
        normal3.Normalize()
        if normal1.Dot(vec2) < 0:
            normal1 = - normal1
        if normal2.Dot(vec3) < 0:
            normal2 = - normal2
        if normal3.Dot(vec1) < 0:
            normal3 = - normal3
        linear_element_1 = LinearElement(self.points[0], self.points[1], normal1)
        linear_element_2 = LinearElement(self.points[1], self.points[2], normal2)
        linear_element_3 = LinearElement(self.points[2], self.points[0], normal3)
        return [linear_element_1, linear_element_2, linear_element_3]
    
    def _form_functions(self):
        a = Matrix33(1, self.points[0][0], self.points[0][1],
                     1, self.points[1][0], self.points[1][1],
                     1, self.points[2][0], self.points[2][1])
        try :
            inv_a = a.inverse()
        except ValueError:
            print('buggy element area', self._area())
            raise FlatElementError('form function bug')
        x1 = inv_a.vector_multiplication(vm.X3D)
        x2 = inv_a.vector_multiplication(vm.Y3D)
        x3 = inv_a.vector_multiplication(vm.Z3D)
       
        return x1, x2, x3
    def _quadratic_form_functions(self):
        a = [[1, self.points[0][0], self.points[0][1],self.points[0][0]**2,self.points[0][0]*self.points[0][1],self.points[0][1]**2],
              [1, self.points[1][0], self.points[1][1],self.points[1][0]**2,self.points[1][0]*self.points[1][1],self.points[1][1]**2],
              [1, self.points[2][0], self.points[2][1],self.points[2][0]**2,self.points[2][0]*self.points[2][1],self.points[2][1]**2],
              [1, self.points[3][0], self.points[3][1],self.points[3][0]**2,self.points[3][0]*self.points[3][1],self.points[3][1]**2],
              [1, self.points[4][0], self.points[4][1],self.points[4][0]**2,self.points[4][0]*self.points[4][1],self.points[4][1]**2],
              [1, self.points[5][0], self.points[5][1],self.points[5][0]**2,self.points[5][0]*self.points[5][1],self.points[5][1]**2]]
                     
    
        try :
            inv_a = a.inverse()
        except ValueError:
            self.plot()
            print(self._area())
            raise FlatElementError('form function bug')
        x1 = inv_a.dot([1,0,0,0,0,0])
        x2 = inv_a.dot([1,0,0,0,0,0])
        x3 = inv_a.dot([1,0,0,0,0,0])
        x4=inv_a.dot([1,0,0,0,0,0])
        
        return x1, x2, x3
    def _area(self):
        u = self.points[1] - self.points[0]
        v = self.points[2] - self.points[0]
        return abs(u.Cross(v)) / 2
        
    def point_belongs(self, point):
        polygon = vm.Polygon2D(self.points)
        point_belongs = polygon.PointBelongs(point)
        return point_belongs
    
    def rotation(self, center, angle, copy=True):
        if copy:
            return TriangularElement([pt.Rotation(center, angle, copy=True) for pt in self.points])
        else:
            for pt in self.points:
                pt.Rotation(center, angle, copy=False)
                
    def translation(self, offset, copy=True):
        if copy:
            return TriangularElement([pt.Translation(offset, copy=True) for pt in self.points])
        else:
            for pt in self.points:
                pt.Translation(offset, copy=False)
                
    def axial_symmetry(self, line, copy=True):
        p1, p2 = line.points
        symmetric_points = []
        for point in self.points:
            u = p2 - p1
            t = (point-p1).Dot(u) / u.Norm()**2
            projection = p1 + t * u
            symmetric_point = vm.Point2D((2 * projection - point).vector)
            symmetric_points.append(symmetric_point)
        if copy: 
            return TriangularElement(symmetric_points)
        else:
            for i, point in enumerate(self.points):
                point = symmetric_points[i]
           
    def line_equation(self,P0:vm.Point2D,P1:vm.Point2D,M:vm.Point2D):
    
        return (P1[0]-P0[0])*(M[1]-P0[1])-(P1[1]-P0[1])*(M[0]-P0[0])  
        
 
    def is_inside_triangle(self,M:vm.Point2D):
        P0=self.points[0]
        P1=self.points[1]
        P2=self.points[2]
        return self.line_equation(P0,P1,M)> 0 and self.line_equation(P1,P2,M) > 0 and self.line_equation(P2,P0,M) > 0
    
    
    def triangle_to_polygon(self):
        points=self.points
        return vm.Polygon2D(points)
    
    def common_edge(self,triangle:'TriangularElement'):
        common_edge=[]
        for point_1 in self.points:
            for point_2 in triangle.points:
                if point_1==point_2:
                    common_edge.append(point_1)
        if len(common_edge)==2:
            return True
        return False
    
    
           
    def common_vertice(self,triangle:'TriangularElement'):
        linear_elements_1=self.linear_elements
        linear_elements_2=triangle.linear_elements
        common_vertice=[]
        for linear_1 in linear_elements_1:
            for linear_2 in linear_elements_2:
                if linear_1==linear_2:
                    common_vertice.append(linear_1)
        if common_vertice!=[]:
            
            return common_vertice[0]
        return None
    
    

           
    def plot(self, ax=None, color='k', width=None, plot_points=False, fill=False):
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
                width=1
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, marker='o', linewidth=width)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=width)
        return ax
    
    
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
            return Mesh([elem.rotation(center, angle, copy=True) for elem in self.elements])
        else:
            for elem in self.elements:
                elem.rotation(center, angle, copy=False)
                
    def translation(self, offset, copy=True):
        if copy:
            return Mesh([elem.translation(offset, copy=True) for elem in self.elements])
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
    def __init__(self, elements_groups:List[ElementsGroup]):
        self.elements_groups = elements_groups
        self.nodes = self._set_nodes_number()
        self.node_to_index = {self.nodes[i]:i for i in range(len(self.nodes))}
        
        DessiaObject.__init__(self, name='')
    
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
        indexes={}
        for node in self.nodes:
           
                indexes[node]=[2*self.node_to_index[node],2*self.node_to_index[node]+1]
        return indexes
    
    def boundary_dict(self):
        boundary_dict = {}
        for elements_group1, elements_group2 in combinations(self.elements_groups, 2):
            linear_elements1 = []
            linear_elements2 = []
            for element in elements_group1.elements:
                linear_elements1.extend(element.linear_elements)
            for element in elements_group2.elements:
                linear_elements2.extend(element.linear_elements)
            duplicate_linear_elements = find_duplicate_linear_element(linear_elements1, linear_elements2)
            if duplicate_linear_elements:
                boundary_dict[(elements_group1, elements_group2)] = duplicate_linear_elements
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
                c1 = vm.Contour2D([vm.LineSegment2D(element.points[0], element.points[1])])
                c2 = vm.Contour2D([vm.LineSegment2D(element.points[1], element.points[2])])
                c3 = vm.Contour2D([vm.LineSegment2D(element.points[2], element.points[0])])
                plot_datas.append(c1.plot_data())
                plot_datas.append(c2.plot_data())
                plot_datas.append(c3.plot_data())
                # plot_datas.extend([c1, c2, c3])
        return plot_datas
    
    
    
    def plot_displaced_mesh(self,node_displacement:Dict[vm.Point2D,List[float]],ax=None,amplification=0.5):
        
        deformed_mesh=self.copy()
        nodes=deformed_mesh.nodes
  
        for node in nodes:
           for displaced_node in node_displacement:
               if node==displaced_node:
                   node[0]+=amplification*node_displacement[displaced_node][0]
                   node[1]+=amplification*node_displacement[displaced_node][1]
            
        ax = deformed_mesh.plot(ax=ax) 
        ax.set_aspect('equal')           
    
        return ax
    
class Mesher(DessiaObject):
    
    def __init__(self,interior_contours:List[vm.Contour2D],exterior_contours:List[vm.Contour2D],triangles:List[TriangularElement],nodes_len:float):
        self.nodes_len=nodes_len
        self.interior_contours=interior_contours
        self.exterior_contours=exterior_contours
        self.triangles=triangles
    
        
    def  neighbour_edge(self,n:int,i:int,di:int):
        return (i+di)%n  
   
    
    def edge_max_distance(self,polygone:vm.Polygon2D,P0:vm.Point2D,P1:vm.Point2D,P2:vm.Point2D,indexes:List[float]):
        n=len(polygone.points)
        distance=0
        j=None
        triangle=vm.Triangle2D([P0,P1,P2])
        for i in range(n):
            if not (i in indexes):
                M=polygone.points[i]
                if triangle.is_inside_triangle(M):
                    d=abs(triangle.line_equation(P1,P2,M))
                    if d > distance:
                        distance=d
                        j=i
       
        return j
    
    
    
    
    def left_edge(self,polygone:vm.Polygon2D):
        n=len(polygone.points)
        x=polygone.points[0][0]
        j=0
        for i in range (1,n):
            if polygone.points[i][0] < x:
                x=polygone.points[i][0]
                j=i
        return j
            
    def new_polygon(self,polygone:vm.Polygon2D,i_beg:int,i_end=int):
        n=len(polygone.points)
        u=[]
        i=i_beg
        while i!=i_end:
            u.append(polygone.points[i])
            i=self.neighbour_edge(n,i,1)
        u.append(polygone.points[i_end])
        p=vm.Polygon2D(u)
        
        return p
    
    def triangulation_polygone_recursive(self,polygone:vm.Polygon2D):
        
        n=len(polygone.points)
        j0=self.left_edge(polygone)
      
       
        j1=self.neighbour_edge(n,j0,1)
       
        j2=self.neighbour_edge(n,j0,-1)
        
        P0=polygone.points[j0]
        P1=polygone.points[j1]
        P2=polygone.points[j2]
        j=self.edge_max_distance(polygone,P0,P1,P2,[j0,j1,j2])
        
        if j==None:
            triangle=vm.Triangle2D([P0,P1,P2])
            self.triangles.append(triangle)
            polygone_1=self.new_polygon(polygone,j1,j2)
                            
            if len(polygone_1.points)==3:
                    
                    new_triangle=vm.Triangle2D([polygone_1.points[0],polygone_1.points[1],
                                                polygone_1.points[2]])    
                    self.triangles.append(new_triangle)  
                                          
            else :
             
                 self.triangulation_polygone_recursive(polygone_1)
     
        else : 
            
            polygone_1=self.new_polygon(polygone,j0,j)
            polygone_2=self.new_polygon(polygone,j,j0)    
            
            if len(polygone_1.points)==3:
                new_triangle=vm.Triangle2D([polygone_1.points[0],polygone_1.points[1],
                                            polygone_1.points[2]])                
                self.triangles.append(new_triangle)
               
            else :
                self.triangulation_polygone_recursive(polygone_1)
                
            if len(polygone_2.points)==3:
                new_triangle=vm.Triangle2D([polygone_2.points[0],polygone_2.points[1],
                                            polygone_2.points[2]])                    
                self.triangles.append(new_triangle)

            else :
                    
                self.triangulation_polygone_recursive(polygone_2)
            
        return self.triangles 
    
    
    
    def _is_convex(self,p1:vm.Point2D, p2:vm.Point2D, p3:vm.Point2D):
        return self._triangle_sum(p1[0], p1[1], p2[0], p2[1], p3[0],p3[1]) < 0
    
    def _is_clockwise(self,polygon:vm.Polygon2D):
        s = 0
        polygon_count = len(polygon.points)
        for i in range(polygon_count):
            point = polygon.points[i]
            point2 = polygon.points[(i + 1) % polygon_count]
            s += (point2[0] - point[0]) * (point2[1] + point[1])
        return s > 0
    
    def _triangle_sum(self,x1, y1, x2, y2, x3, y3):
        return x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)
    
    def _contains_no_points(self,p1:vm.Point2D,p2:vm.Point2D,p3:vm.Point2D, polygon:vm.Polygon2D):
       triangle=vm.Triangle2D([p1,p2,p3])
       for pn in polygon.points:
            if pn in [p1, p2, p3]:
                continue
            elif triangle.is_inside_triangle(pn):
                return False
       return True

    def _is_ear(self,p1:vm.Point2D,p2:vm.Point2D,p3:vm.Point2D, polygon:vm.Polygon2D):
        triangle=vm.Triangle2D([p1,p2,p3])
        ear = self._contains_no_points(p1, p2, p3, polygon) and \
            self._is_convex(p1, p2, p3) and \
            triangle.area > 0
        return ear
    
    def earclip(self,polygon:vm.Polygon2D):
       
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
            
            p1=vm.Point2D([prev_point[0], prev_point[1]])
            p2=vm.Point2D([ear[0], ear[1]])
            p3= vm.Point2D([next_point[0], next_point[1]])
            triangle=vm.Triangle2D([p1,p2,p3])  
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
                 
                                 
        
    def basic_triangulation(self,polygon1:vm.Polygon2D,polygon2:vm.Polygon2D,
                           segment_to_nodes:Dict[vm.LineSegment2D,List[vm.Point2D]]):
        triangles=[]
        
        for j in range(len(polygon1.line_segments)):
                    
            pj=segment_to_nodes[polygon1.line_segments[j]]
            qj=segment_to_nodes[polygon2.line_segments[j]]
            u=len(pj)
            v=len(qj)
            if u==2 and v==2 : 
                new_triangle_1=vm.Triangle2D([pj[0],pj[1],qj[0]])
                triangles.append(new_triangle_1)                          
                new_triangle_2=vm.Triangle2D([pj[1],qj[1],qj[0]])
                triangles.append(new_triangle_2)                           
            if u>=v:
                for i in range(v-1):
                    new_triangle_1=vm.Triangle2D([pj[i+1],pj[i],qj[i]])
                   
                    triangles.append(new_triangle_1)
                    new_triangle_2=vm.Triangle2D([qj[i+1],pj[i+1],qj[i]])
                   
                    triangles.append(new_triangle_2)
                              
                for  i in range(v-1,u-1):
                                
                     new_triangle=vm.Triangle2D([pj[i],qj[v-1],pj[i+1]])
                                     
                     triangles.append(new_triangle)  
            else :
                  for i in range(u-1):
                    new_triangle_1=vm.Triangle2D([qj[i+1],qj[i],pj[i]])
                   
                    triangles.append(new_triangle_1)
                    new_triangle_2=vm.Triangle2D([pj[i+1],qj[i+1],pj[i]])
                   
                    triangles.append(new_triangle_2)
                              
                  for  i in range(u-1,v-1):
                                
                     new_triangle=vm.Triangle2D([qj[i],pj[u-1],qj[i+1]])
                                     
                     triangles.append(new_triangle)  
                             
        return triangles                            
    
    def alternative_triangulation(self,polygon:vm.Polygon2D,interior_polygon:vm.Polygon2D,empty:bool,split:bool,far:bool):
        # ax=plt.subplot()
        all_aspect_ratios=[]
        segment_to_nodes={}
        
        all_offsets=[]
        all_meshes=[]
        offset_values=[]
        good_meshes=[]
        p=9
        k=3
        bound=polygon.bounding_rectangle2()
        offset_len=bound.min_length()
        while  k<8 :
            repair = False
            polygon_offsets=[]
            polygon_offsets+=[polygon]
            while repair is False :
                
                new_polygon=polygon.Offset(-p*offset_len/(10*k))
                
                if not new_polygon.SelfIntersect()[0] :
                    if new_polygon.Area() > polygon.Area():
                        polygon_offsets.append(polygon.Offset(p*offset_len/(10*k)))
                        offset_values.append(p*offset_len/(10*k))
                        # polygon.Offset(9*offset_len/(10*k)).MPLPlot()
                    else :
                         offset_values.append(-p*offset_len/(10*k))
                         polygon_offsets.append(new_polygon)
                         
                    repair=True
                    
                    
                if new_polygon.SelfIntersect()[0] :
                    
                    polygon_offsets.append(new_polygon.select_reapaired_polygon([]))
                    # new_polygon.select_reapaired_polygon([]).MPLPlot()
                    repair=True
                  
            all_offsets.append(polygon_offsets)        
            k=k+1
                    
        for polygon_offsets in all_offsets:
            offset_triangles=[]   
            for polygon in polygon_offsets:
               
                for segment in polygon.line_segments:
                    segment_to_nodes[segment]=segment.discretise(self.nodes_len,None)
                    
                    
            if len(polygon_offsets)>2:      
                for k in range(len(polygon_offsets)-2):
                
                    if len(polygon_offsets[k].line_segments)==len(polygon_offsets[k+1].line_segments):
                        offset_triangles+=self.basic_triangulation(polygon_offsets[k],
                                                            polygon_offsets[k+1],segment_to_nodes)
                  
                    else :
                    
                        offset_triangles+=self.mesh_in_between(polygon_offsets[k+1],
                                                            polygon_offsets[k],True)
                    
                     
            l=len(polygon_offsets)
          
            if len(polygon_offsets[-1].line_segments)==len(polygon_offsets[l-2].line_segments):
               offset_triangles+=self.basic_triangulation(polygon_offsets[l-2],polygon_offsets[-1],segment_to_nodes)
               if  not empty :
                    
                    last_points=[]
                    for segment in polygon_offsets[-1].line_segments:
                        for point in segment_to_nodes[segment]:
                            if point not in last_points:
                                last_points.append(point)
                    last_polygon=vm.Polygon2D(last_points)
                    ear=self.earclip(last_polygon)
                    for triangle in ear : 
                          if triangle.area<10e-9:
                             ear.remove(triangle)
                             
                    offset_triangles+=ear
               else :
                    if far :
                        offset_triangles+=self.mesh_in_between(interior_polygon,polygon_offsets[-1],True)
                               
                  
                                 
                        
                         
            else:
                offset_triangles+=self.mesh_in_between(polygon_offsets[-1],polygon_offsets[l-2],empty)                
              
            all_meshes.append(offset_triangles)                             
         
       
        for mesh in all_meshes:
            if self.triangulation_max_aspect_ratio(mesh)!=0:
                all_aspect_ratios.append(self.triangulation_max_aspect_ratio(mesh))
                good_meshes.append(mesh)
            
            
        index=all_aspect_ratios.index(min(all_aspect_ratios))
        
       
        # plt.plot(offset_values, all_aspect_ratios) 
          
        # plt.xlabel('offset values') 
        # plt.xlim(-1,1) 
        # plt.ylabel(' maximum aspect ratio') 

        
        # plt.show() 
        
        return good_meshes[index]

    def mesh_in_between(self,in_polygon:vm.Polygon2D,out_polygon:vm.Polygon2D,empty:bool):
        
        # ax=plt.subplot()
        projection_points=[]
        segment_to_nodes={}
        closest_segment={}
        all_triangles=[]
        out_point_image={}

        for segment in out_polygon.line_segments:
            segment_to_nodes[segment]=segment.discretise(self.nodes_len,None)
         
                
        for segment in out_polygon.line_segments:
            for point in segment_to_nodes[segment]:
                out_point_image[point]=[]
        for segment in out_polygon.line_segments:
            projection_points.append([]) 
            
        for segment in in_polygon.line_segments:
        
            segment_to_nodes[segment]=segment.discretise(0,None)
       
                                
        for out_segment in out_polygon.line_segments:
            
            index_0=out_polygon.line_segments.index(out_segment)
            
            mid=out_segment.PointAtCurvilinearAbscissa(out_segment.Length()/2)
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
                projection=near_segment.PointProjection2(point)
                for in_segment in in_polygon.line_segments:
                    
                    if point.point_distance(in_segment.point1) < point.point_distance(in_segment.point2):
                        
                        d_1.append(point.point_distance(in_segment.point1))
                        i.append(0)
                    elif  point.point_distance(in_segment.point1) > point.point_distance(in_segment.point2):
                        d_1.append(point.point_distance(in_segment.point2))
                        i.append(1)
                    elif point.point_distance(in_segment.point1) == point.point_distance(in_segment.point2):
                         
                          d_1.append(point.point_distance(in_segment.point1))
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
                            
                            projection_points[out_point_image[point][0][0]].remove(out_point_image[point][0][1])
                            projection_points[out_point_image[point][0][0]].append(new_proj)
                            
                          if new_proj not in projection_points[index_0]:
                              projection_points[index_0].insert(index_point,new_proj)
                  
                                 
                           
                          
                else:                                               
                    if projection is None :
                            
                           
                            if out_point_image[point]==[]:
                              
                                out_point_image[point].append([index_0,new_proj])
                                  
                                if new_proj not in projection_points[index_0]:
                                   
                                    projection_points[index_0].insert(index_point,new_proj)
                                
                            else :
                                  
                                 if out_point_image[point][0][1] != new_proj :
                                      
                                      projection_points[out_point_image[point][0][0]].remove(out_point_image[point][0][1])
                                      projection_points[out_point_image[point][0][0]].append(new_proj)
                                      if new_proj not in projection_points[index_0]:
                                            projection_points[index_0].insert(index_point,new_proj)
                                     
                                 else :
                                     if new_proj not in projection_points[index_0]:
                                        projection_points[index_0].insert(index_point,new_proj)
                    else :
                        if out_point_image[point]==[]:
                          i=0
                          if projection.point_distance(near_segment.point1) > projection.point_distance(near_segment.point2) :
                             i+=1
                          l = near_segment.Length()/10
                          if projection.point_distance(near_segment.points[i]) >= l :
                           
                            out_point_image[point].append([index_0,projection]) 
                            projection_points[index_0].insert(index_point,projection)
                          else :
                               
                               if new_proj not in projection_points[index_0]:     
                                   projection_points[index_0].insert(index_point,new_proj)
                                      
                        else :
                             if new_proj not in projection_points[index_0]:     
                                   projection_points[index_0].insert(index_point,new_proj)
                  
        for out_segment in out_polygon.line_segments:
            
            index_0=out_polygon.line_segments.index(out_segment)
            v=len(segment_to_nodes[out_segment])
            u=len(projection_points[index_0])

            if u>=v and u>2:
              
                for j in range(v-1):
                    
                    new_triangle=vm.Triangle2D([projection_points[index_0][j],projection_points[index_0][j+1],
                                        segment_to_nodes[out_segment][j]])   
                    # new_triangle=vm.Triangle2D([projection_points[index_0][j],segment_to_nodes[out_segment][j+1],
                    #                     segment_to_nodes[out_segment][j]])   
          
                    all_triangles.append(new_triangle)
                    new_triangle_0=vm.Triangle2D([projection_points[index_0][j+1],segment_to_nodes[out_segment][j],
                                        segment_to_nodes[out_segment][j+1]])   
                    # new_triangle_0=vm.Triangle2D([projection_points[index_0][j+1],projection_points[index_0][j],
                    #                     segment_to_nodes[out_segment][j+1]])
                    all_triangles.append(new_triangle_0)
                    
                for j in  range(v-1,u-1):
                      new_triangle=vm.Triangle2D([projection_points[index_0][j],projection_points[index_0][j+1],
                                      segment_to_nodes[out_segment][v-1]])   
                      all_triangles.append(new_triangle)
                      
           
            if u<v :
                for j in range(u-1):
                  
                   new_triangle=vm.Triangle2D([projection_points[index_0][j],segment_to_nodes[out_segment][j+1],
                                      segment_to_nodes[out_segment][j]])   
                                
                   all_triangles.append(new_triangle)
                   new_triangle_0=vm.Triangle2D([projection_points[index_0][j+1],projection_points[index_0][j],
                                       segment_to_nodes[out_segment][j+1]])   
                   # new_triangle_0=vm.Triangle2D([projection_points[index_0][j],segment_to_nodes[out_segment][j],
                   #                    segment_to_nodes[out_segment][j+1]]) 
                   
                   all_triangles.append(new_triangle_0)
                   
                for j in  range(u-1,v-1):
                   new_triangle=vm.Triangle2D([segment_to_nodes[out_segment][j],  segment_to_nodes[out_segment][j+1],
                                  projection_points[index_0][u-1]])   
                   all_triangles.append(new_triangle)
                   
            if u==2 and v==2:
                if projection_points[index_0][1] != projection_points[index_0][0] :
                    new_triangle_1=vm.Triangle2D([projection_points[index_0][0],projection_points[index_0][1],
                                       segment_to_nodes[out_segment][1]]) 
                    all_triangles.append(new_triangle_1)
                  
                new_triangle_2=vm.Triangle2D([projection_points[index_0][0],segment_to_nodes[out_segment][0],
                                    segment_to_nodes[out_segment][1]]) 
               
                all_triangles.append(new_triangle_2)
           
                
           
                                
  
        if empty is False:
            
            last_points=[]
            for k in range(len(projection_points)):
                for point in projection_points[k]:
                    if point not in last_points:
                        last_points.append(point)
            last_polygon=vm.Polygon2D(last_points)
            ear=self.earclip(last_polygon)
            for triangle in ear : 
                if triangle.area<10e-9:
                    ear.remove(triangle)
            all_triangles+=ear
        
        return all_triangles
  
    
    
    def polygon_to_triangles(self,polygons:List[vm.Polygon2D]):
        triangles=[]
        for polygon in polygons:
            triangles.append(vm.Triangle2D(polygon.points))
        return triangles

    
    def triangulation_max_aspect_ratio(self,triangles:List[vm.Triangle2D]):
        all_aspect_ratios=[]
        for triangle in triangles:
            if triangle.area<10E-9:
                return 0
            else :
                all_aspect_ratios.append(triangle.aspect_ratio())
        index=all_aspect_ratios.index(max(all_aspect_ratios))
        
        return all_aspect_ratios[index]
    def triangulation_min_aspect_ratio(self,triangles:List[vm.Triangle2D]):
        all_aspect_ratios=[]
        for triangle in triangles:
            all_aspect_ratios.append(triangle.aspect_ratio())
        index=all_aspect_ratios.index(min(all_aspect_ratios))
        
        return all_aspect_ratios[index]
    def generate_mesh(self,min_aspect_ratio:float,split:bool):
        
        ax=plt.subplot()
        segment_to_nodes={}
        all_segments=set()
        triangles=[]
        interior_polygons=[]
        exterior_polygons=[]
        split_polygons=[]
              
        if self.interior_contours :
            
            for contour in self.interior_contours:
                
                polygon_points=[]
                if contour.primitives[0].__class__.__name__ == 'Circle2D':
                    
                    for point in contour.primitives[0].discretise(self.nodes_len,None):
                        if point not in polygon_points:
                            polygon_points.append(point)
                    polygon=vm.Polygon2D(polygon_points)  
                    interior_polygons.append(polygon)
                    
                else :
                    
                    for primitive in contour.primitives:
                        if isinstance(primitive,vm.LineSegment2D):
                            
                            if primitive.point1 not in polygon_points:
                                polygon_points.append(primitive.point1)
                            
                            # if primitive.point2 not in polygon_points:
                            #     polygon_points.append(primitive.point2)
                        else :
                              print(primitive)
                              for point in primitive.discretise(self.nodes_len,None):
                                  if point not in polygon_points:
                                      polygon_points.append(point)
                    polygon=vm.Polygon2D(polygon_points)          
                    interior_polygons.append(polygon) 
                     
      
       
        for contour in self.exterior_contours:
          polygon_points=[]
          if isinstance(contour,vm.Circle2D):
              for point in contour.discretise(self.nodes_len,None):
                  if point not in polygon_points:
                      polygon_points.append(point)
              polygon=vm.Polygon2D(polygon_points)          
              exterior_polygons.append(polygon)
          else :
              for primitive in contour.primitives:
                  if isinstance(primitive,vm.LineSegment2D):
                      
                      if primitive.point1 not in polygon_points:
                          polygon_points.append(primitive.point1)
                      if primitive.point2 not in polygon_points:
                          polygon_points.append(primitive.point2)   
                  else :
                        for point in primitive.discretise(self.nodes_len,None):
                            if point not in polygon_points:
                                polygon_points.append(point)
              polygon=vm.Polygon2D(polygon_points)          
              exterior_polygons.append(polygon)   
              
    
         
    
        if split is True :
          
            Surface=vm.Surface2D(self.exterior_contours[-1],self.interior_contours)
            split_contours = Surface.split_regularly(2)
             
            if self.interior_contours : 
                for contour in split_contours:
                    polygon_points=[]
                    if isinstance(contour,vm.Circle2D):
                        for point in contour.discretise(self.nodes_len,None):
                            if point not in polygon_points:
                                polygon_points.append(point)
                        polygon=vm.Polygon2D(polygon_points)          
                        split_polygons.append(polygon)
                    else :
                        for primitive in contour.primitives:
                            if isinstance(primitive,vm.LineSegment2D):
                                
                                if primitive.point1 not in polygon_points:
                                    polygon_points.append(primitive.point1)
                                
                                if primitive.point2 not in polygon_points:
                                    polygon_points.append(primitive.point2)
                            else :
                                  for point in primitive.discretise(self.nodes_len,None):
                                      if point not in polygon_points:
                                          polygon_points.append(point)
                        polygon=vm.Polygon2D(polygon_points)          
                        split_polygons.append(polygon) 
            # split_polygons[0].MPLPlot(ax=ax)
            # rec=split_polygons[0].bounding_rectangle2()
            # offset_len=rec.min_length()
        
            # offset=split_polygons[0].Offset(offset_len/15)
            # print(offset.SelfIntersect()[0])
            # offset.select_reapaired_polygon([]).MPLPlot()            
            for s in split_polygons :
                triangles+=self.alternative_triangulation(split_polygons[0],None,True,True,False)
                
            
        else :
                      
            """
            
            The convention adopted is the following : The last polygon of 
            interior_polygons(which stays empty) is inside the last polygon of 
            exterior_polygon 
            
            """ 
             
            
            # rec=exterior_polygons[0].bounding_rectangle2()
            # offset_len=rec.min_length()
            # offset=exterior_polygons[0].Offset(offset_len/40)
            # offset.MPLPlot()
            # w=[exterior_polygons[0],exterior_polygons[1]]
            # for p in w:
         
            #     triangles+=self.alternative_triangulation(p,None,False,False,False)
            if self.interior_contours:
                if exterior_polygons[-1].polygon_distance(interior_polygons[-1]) > exterior_polygons[-1].max_length()/3 :
                    
                    triangles+=self.alternative_triangulation(exterior_polygons[-1],interior_polygons[-1],True,False,True)
                else :                
                    
                    triangles+=self.mesh_in_between(interior_polygons[-1],exterior_polygons[-1],True)
            
                if len(exterior_polygons)>1:
                    for polygon in exterior_polygons:
                      
                        if polygon != exterior_polygons[-1]:
                            
                            polygon_1=polygon.copy()
                           
                            possible_triangles=self.earclip(polygon)
                            
                            """For one polygon,earclip created a LineSegment2D as an instance of Triangle2D
                            for unknown reasons,hence the area test below """
                            
                            for triangle in possible_triangles:
                                if triangle.area==0:
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
                                        triangles+=self.alternative_triangulation(polygon_1,None,False,False,False)
                                        
                                      else :
                                      
                                          k=k+1
                                          
                                if Next==True :
                                    triangles+=possible_triangles
                          
            else :
                for polygon in exterior_polygons:
                      
               
                    
                    polygon_1=polygon.copy()
                   
                    possible_triangles=self.earclip(polygon)
                    
                    """For one polygon,earclip created a LineSegment2D as an instance of Triangle2D
                    for unknown reasons,hence the area test below """
                    
                    for triangle in possible_triangles:
                        if triangle.area==0:
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
                                triangles+=self.alternative_triangulation(polygon_1,None,False,False,False)
                                
                              else :
                              
                                  k=k+1
                                  
                        if Next==True :
                            triangles+=possible_triangles

        all_triangles=[]
        all_triangles+=triangles
        all_triangle_elements=[]
        plot_aspect_ratio_triangles=[]
        all_aspect_ratios={}
        
        # for triangle in all_triangles:
        #       triangle.MPLPlot(ax=ax)

        for triangle in triangles:
            
            all_segments= all_segments.union(triangle.line_segments)
            
        for segment in all_segments:
            
            segment_to_nodes[segment]=segment.discretise(self.nodes_len,None)
                
        for triangle in triangles :
            
            meshing=triangle.mesh_triangle(segment_to_nodes,self.nodes_len,None)
            plot_aspect_ratio_triangles+=meshing[0]
            all_triangles+=meshing[0]
            all_aspect_ratios.update(meshing[1])
            
        for triangle in all_triangles:
            
            # triangle.MPLPlot(ax=ax)
            triangular_element=TriangularElement(triangle.points)
            all_triangle_elements.append(triangular_element)
            
        self.plot_aspect_ratio(plot_aspect_ratio_triangles,all_aspect_ratios,ax)
        ax.set_aspect('equal')
        return all_triangle_elements
    
    def plot_aspect_ratio(self,all_triangles:List[vm.Triangle2D],
                          all_aspect_ratios:Dict[vm.Triangle2D,float],
                          ax,min_aspect_ratio=None,max_aspect_ratio=None):
        
        """
        Plots the mesh with colored triangles representing the \
        value of the aspect ratio. 
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
  

                
        
                

          
      
    
        
               
           
           
       
       
            
        
        
    
  
        
                    
            
       
        
    
    
    
    