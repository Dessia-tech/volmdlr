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
class FlatElementError(Exception):
    pass

def find_duplicate_linear_element(linear_elements1, linear_elements2):
    duplicates = []
    for linear_element in linear_elements1:
        if linear_element in linear_elements2 and linear_element not in duplicates:
            duplicates.append(linear_element)
    return duplicates


class LinearElement(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    def __init__(self, points:Tuple[vm.Point2D, vm.Point2D],
                 interior_normal:vm.Vector2D):
        self.points = points
        self.interior_normal = interior_normal
        
        DessiaObject.__init__(self, name='')
        
    def __hash__(self):
        return self.points[0].__hash__() + self.points[1].__hash__()
        
    def __eq__(self, other_linear_element):
        if self.__class__ != other_linear_element.__class__:
            return False
        return (self.points[0] == other_linear_element.points[0]
                and self.points[1] == other_linear_element.points[1]) \
            or (self.points[0] == other_linear_element.points[1]
                and self.points[1] == other_linear_element.points[0])
            
    def length(self):
        return self.points[1].point_distance(self.points[0])
    
    def plot(self, ax=None, color='k', width=None, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        if width is None:
            width=1
        if plot_points:
            ax.plot([self.points[0][0], self.points[1][0]],
                    [self.points[0][1], self.points[1][1]],
                    color=color, marker='o', linewidth=width)
        else:
            ax.plot([self.points[0][0], self.points[1][0]],
                    [self.points[0][1], self.points[1][1]],
                    color=color, linewidth=width)
        return ax


class TriangularElement(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points:Tuple[vm.Point2D, vm.Point2D, vm.Point2D]):
        self.points = points
        
        self.linear_elements = self._to_linear_elements()
        self.form_functions = self._form_functions()
        
        self.center = (self.points[0]+self.points[1]+self.points[2])/3
        
        self.area = self._area()
        
        DessiaObject.__init__(self, name='')
        
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
        linear_element_1 = LinearElement([self.points[0], self.points[1]],
                                         normal1)
        linear_element_2 = LinearElement([self.points[1], self.points[2]],
                                         normal2)
        linear_element_3 = LinearElement([self.points[2], self.points[0]],
                                         normal3)
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

    def quadratic_form_functions(self):
        a = [[1, self.points[0][0], self.points[0][1], self.points[0][0]**2,
              self.points[0][0]*self.points[0][1], self.points[0][1]**2],
             [1, self.points[1][0], self.points[1][1], self.points[1][0]**2,
              self.points[1][0]*self.points[1][1], self.points[1][1]**2],
             [1, self.points[2][0], self.points[2][1], self.points[2][0]**2,
              self.points[2][0]*self.points[2][1], self.points[2][1]**2],
             [1, self.points[3][0], self.points[3][1], self.points[3][0]**2,
              self.points[3][0]*self.points[3][1], self.points[3][1]**2],
             [1, self.points[4][0], self.points[4][1], self.points[4][0]**2,
              self.points[4][0]*self.points[4][1], self.points[4][1]**2],
             [1, self.points[5][0], self.points[5][1], self.points[5][0]**2,
              self.points[5][0]*self.points[5][1], self.points[5][1]**2]]

        try:
            inv_a = a.inverse()
        except ValueError:
            self.plot()
            print(self._area())
            raise FlatElementError('form function bug')
        x1 = inv_a.dot([1, 0, 0, 0, 0, 0])
        x2 = inv_a.dot([1, 0, 0, 0, 0, 0])
        x3 = inv_a.dot([1, 0, 0, 0, 0, 0])
        x4 = inv_a.dot([1, 0, 0, 0, 0, 0])
        
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
           
    def line_equation(self, P0: vm.Point2D, P1: vm.Point2D, M: vm.Point2D):
        return (P1[0]-P0[0])*(M[1]-P0[1])-(P1[1]-P0[1])*(M[0]-P0[0])

    def is_inside_triangle(self, M: vm.Point2D):
        P0 = self.points[0]
        P1 = self.points[1]
        P2 = self.points[2]
        return self.line_equation(P0, P1, M) > 0 and \
               self.line_equation(P1, P2, M) > 0 and \
               self.line_equation(P2, P0, M) > 0
                
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
        
        for p1, p2 in zip(self.points, list(self.points[1:])+[self.points[0]]):
            if width is None:
                width = 1
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                        color=color, marker='o', linewidth=width)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                        color=color, linewidth=width)
        return ax
    
    
class ElementsGroup(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, elements: List[TriangularElement], name: str):
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
        indexes = {}
        for node in self.nodes:
            indexes[node] = [2*self.node_to_index[node],
                             2*self.node_to_index[node]+1]
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
   
    def plot(self, ax=None, fill=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        for elements_group in self.elements_groups:
            elements_group.plot(ax=ax, fill=fill)
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

    def plot_displaced_mesh(self, node_displacement: Dict[vm.Point2D, List[float]],ax=None,amplification=0.5):
        deformed_mesh=self.copy()
        nodes=deformed_mesh.nodes
        for node in nodes:
            for displaced_node in node_displacement:
                if node == displaced_node:
                    node[0] += amplification*node_displacement[displaced_node][0]
                    node[1] += amplification*node_displacement[displaced_node][1]
        ax = deformed_mesh.plot(ax=ax)
        return ax
    
class Mesher(DessiaObject):
    
    def __init__(self,contour:vm.Contour2D,triangles:List[TriangularElement],nodes_len:float):
        self.nodes_len=nodes_len
        self.contour=contour
        self.polygon=self.contour._get_polygon()
        self.triangles=triangles
                
        
        
    def neighbour_edge(self,n:int,i:int,di:int):
        return (i+di)%n  
   
    
    def edge_max_distance(self,polygone:vm.Polygon2D,P0:vm.Point2D,P1:vm.Point2D,P2:vm.Point2D,indexes:List[float]):
        n=len(polygone.points)
        distance=0
        j=None
        triangle=TriangularElement([P0,P1,P2])
        for i in range(n):
            if not (i in indexes):
                M=polygone.points[i]
                if triangle.is_inside_triangle(M):
                    d=math.abs(triangle.line_equation(P1,P2,M))
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
            self.triangles.append(TriangularElement([P0,P1,P2]))
           
            polygone_1=self.new_polygon(polygone,j1,j2)
            
            if len(polygone_1.points)==3:
                self.triangles.append(TriangularElement([polygone_1.points[0],polygone_1.points[1],polygone_1.points[2]]))
                
            else :
                
                self.triangulation_polygone_recursive(polygone_1)
  
        else : 
            
            
            polygone_1=self.new_polygon(polygone,j0,j)
            polygone_2=self.new_polygon(polygone,j,j0)    
            
            if len(polygone_1.points)==3:
                self.triangles.append(polygone_1)
            else :
                self.triangulate_polygone_recursive(polygone_1)
                
            if len(polygone_2.points)==3:
                self.triangles.append(polygone_2)
            else :
                self.triangulate_polygone_recursive(polygone_2)
                
        return self.triangles 
                
    def plot_triangulised_contour(self):
        fig,ax=plt.subplots()
        patches=[]
        for triangle in self.triangles:
            triangle.plot(ax=ax)
        return ax
        
    def nodes_on_segments(self,triangles:List[TriangularElement]):
        n=self.nodes_len
        nodes=[]
        node_linear_element={}
        
        for triangle in triangles:
            
            for linear_element in triangle.linear_elements:
               
                l0= int(n*linear_element.length())
                
                
                for k in range(1,l0):

                    node=vm.Point2D([linear_element.points[1][0]*k/l0+(1-k/l0)*linear_element.points[0][0],linear_element.points[1][1]*k/l0+(1-k/l0)*linear_element.points[0][1]])
                    if node not in nodes:
                        nodes.append(node)
                    
                    node_linear_element[node]=linear_element
        random.shuffle(nodes)
            
        return [nodes, node_linear_element]
    
    def closest_neightbours(self,triangles:List[TriangularElement],node:vm.Point2D):
        
        neightbours=[]
        d=[]
        possible_nodes=[]
        nodes=self.nodes_on_segments(triangles)[0]
        node_linear_element=self.nodes_on_segments(triangles)[1]
        k=True
        # while k < len(nodes):
        #     if node_linear_element[node]==node_linear_element[nodes[k]] and node_linear_element[node]==node_linear_element[nodes[k+1]]:
        #         k=k+2
        #     else :
        for p in nodes :
            d_0=node.point_distance(p)
            if d_0!=0:
                    d.append(d_0)
                # d.append(math.sqrt((node[0]-nodes[k+1])**2+(node[1]-nodes[k+1])**2))
                
                    possible_nodes.append(p)
        print(d)
       
            
        while k :
           
                
                index_min_1=d.index(min(d))
                d.pop(index_min_1)
                nodes.pop(index_min_1)
                print(index_min_1)
                index_min_2=d.index(min(d))
                print(index_min_2)
                n1=possible_nodes[index_min_1]
                n2=possible_nodes[index_min_2]
                if n1!=n2:
                    if (n1[0]!= node[0] and n2[0]!=node[0]) or (n1[1]!= node[1] and n2[1]!=node[1]) :
                        
                           
                           
                    
                    # if node_linear_element[possible_nodes[index_min_1]]!=node_linear_element[node] and node_linear_element[possible_nodes[index_min_2]]!=node_linear_element[node]:
                           
                           
                           neightbours.append(possible_nodes[index_min_1])
                           neightbours.append(possible_nodes[index_min_2])
                           k=False
                           print('ok')
                    
                else : 
                    
                    k=True
                    
        return neightbours
                
        
                
            
    def assemble_mesh(self,triangles:List[TriangularElement],trigger):
        memo=[]
        for triangle in triangles: 
           l0=triangle.linear_elements[0].length()/3
           if l0 > trigger :
               
               nodes=self.nodes_on_segments(triangles)[0]
               
               
               node_linear_element=self.nodes_on_segments(triangles)[1]
               node_counter=[]
               segment_counter=[]
               new_triangles=[]
               
               for node in nodes :
                   print(node)
                  
                   n1=self.closest_neightbours(triangles,node)[0]
                   n2=self.closest_neightbours(triangles,node)[1]
                   
                   triangle = TriangularElement([node,n1,n2])
                   
                   new_triangles.append(triangle)
               memo.extend(self.assemble_mesh(new_triangles,trigger)) 
            
          
          
        return memo 
                       
    
          
      
           
               
           
           
       
       
            
        
        
    
  
        
                    
            
       
        
    
    
    
    