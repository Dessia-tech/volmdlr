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
class FlatElementError(Exception):
    pass

def find_duplicate_linear_element(linear_elements1, linear_elements2):
    duplicates = []
    for linear_element in linear_elements1:
        if linear_element in linear_elements2 and linear_element not in duplicates:
            duplicates.append(linear_element)
    return duplicates


class LinearElement(vm.edges.LineSegment2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    def __init__(self, point1:vm.Point2D,point2:vm.Point2D, interior_normal:vm.Vector2D,name=''):

        self.interior_normal = interior_normal
        
        vm.edges.LineSegment2D.__init__(self,point1=point1,point2=point2,name=name)
        
    def __hash__(self):
        return self.point1.__hash__() + self.point2.__hash__()
        
    def __eq__(self, other_linear_element):
        if self.__class__ != other_linear_element.__class__:
            return False
        return (self.point1 == other_linear_element.point1 and self.point2 == other_linear_element.point2) \
            or (self.point1 == other_linear_element.point2 and self.point2 == other_linear_element.point1)
            
    def length(self):
        
        return self.point2.point_distance(self.point1)
    
    def closest_middle(self,points:List[vm.Point2D]):
        
        mid=vm.Point2D([(self.point1[0]+self.point2[0])/2,(self.point1[1]+self.point2[1])/2])
        d=[]
        u=[]
        for point in points :
            if d!=0:
                d.append(math.sqrt((mid[0]-point[0])**2+(mid[1]-point[1])**2))
                u.append(point)
        
        index=d.index((min(d)))
        return index
    
    # def mesh_segment(self,n:float):
        
        
    #     segment_to_nodes={}
    
         
    #     nodes=[]
    #     if n*self.length() < 1 :
    #         segment_to_nodes[self]=[]
    #     else :
    #          l0= int(math.ceil(n*self.length()))
        
                    
   
    #          for k in range(l0):
                      
                             
    #              node=self.PointAtCurvilinearAbscissa(k/n)
                       
    #              nodes.append(node)
    #          nodes.insert(len(nodes),self.point2)
                   
    #          segment_to_nodes[self]=nodes
         
                
               
            
    #     return segment_to_nodes
        
               
        
    
    
    
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

class TriangularElement(vm.Polygon2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    def __init__(self, points:Tuple[vm.Point2D, vm.Point2D, vm.Point2D],subdivsion=[]):
        self.points = points
        self.subdivsion=subdivsion
        self.linear_elements = self._to_linear_elements()
        self.form_functions = self._form_functions()
        
        self.center = (self.points[0]+self.points[1]+self.points[2])/3
        
        self.area = self._area()
        
        vm.Polygon2D.__init__(self,points=points, name='')
        
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
            self.plot()
            print('buggy element area', self._area())
            raise FlatElementError('form function bug')
        x1 = inv_a.vector_multiplication(vm.X3D)
        x2 = inv_a.vector_multiplication(vm.Y3D)
        x3 = inv_a.vector_multiplication(vm.Z3D)
       
        return x1, x2, x3
    def quadratic_form_functions(self):
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
    
    
    def min_length(self):
         L=[]
        
         for k in range(len(self.linear_elements)):
             L.append(self.linear_elements[k].length())
       
         return min(L)     
              
    
        
        
                           
    


       
    def closest_neightbours(self,node:vm.Point2D,nodes:List[vm.Point2D]):
        
        
      
        
        neightbours=[]
        d=[]
        possible_nodes=[]
        k=1
      
        for p in nodes :
            d_0=node.point_distance(p)
            if d_0!=0:
                    d.append(d_0)
                
                    possible_nodes.append(p)
        
                
        while k==1 :
          
                if len(d)==2:
                
                     k=0
                 
                index_min_1=d.index(min(d))
                d.pop(index_min_1)
                
                index_min_2=d.index(min(d))
             
                index_3=d.index(max(d))
                
          
                    
                if index_min_1==index_min_2:
                   index_min_2=index_min_2+1
                   
                n1=possible_nodes[index_min_1]
                n2=possible_nodes[index_min_2]
                n3=possible_nodes[index_3]
                
                line=vm.Line2D(n1,n2)
                
                if line.PointProjection(node)!=node:   
                               
                       
                           
                           neightbours.append(n1)
                           neightbours.append(n2)
                      
                           
                           k=0
                    
                else : 
                    if n3!=n1:
                        if vm.Line2D(n3,n1).PointProjection(node)!=node:
                            neightbours.append(n3)
                            neightbours.append(n1)
                           
                           
                            k=0
                    else : 
                        if n3!=n2:
                            if vm.Line2D(n3,n2).PointProjection(node)!=node:
                                neightbours.append(n3)
                                neightbours.append(n2)
                               
                                k=0
                         
                            
                        
                         
                            else : 
                                   k=0
        return neightbours
    
    
    
    
    def triangle_to_polygon(self):
        points=self.points
        return vm.Polygon2D(points)
    
    
    
    def mesh_segments(self,n:float):
        
        
        segment_to_nodes={}
        for linear_element in self.linear_elements:
            
            nodes=[]
            if n*linear_element.length() < 1:
                
                 segment_to_nodes[linear_element]=[linear_element.point1,linear_element.point2]
                
            else :
                   
                   l0= int(math.ceil(n*linear_element.length()))
                    
           
                   for k in range(l0):
                  
                         
                       node=linear_element.PointAtCurvilinearAbscissa(k/n)
                       
                       nodes.append(node)
                   nodes.insert(len(nodes),linear_element.point2)
                   segment_to_nodes[linear_element]=nodes

                
               
            
        return segment_to_nodes
        
      
        
   
    
        
        
    
    
    
    
    
    
    

                            
    #                          # H=[]
    #                          # for k in range(len(new_triangle.linear_elements)):
    #                          #    H.append(2*new_triangle.area/new_triangle.linear_elements[k].length())
                                
    #                          #    E=new_triangle.min_length()
    #                          #    h=min(H)
                               
    #                          #    if E/h<20:
    #                              subdivsion.append(new_triangle)
    #                              subdivsion.extend(new_triangle.mesh_triangle(n=n,trigger=trigger))
                                
                              
                        
    #                 else :
    #                     return [self] 
               
                
                    
    #     else:
    #         return [self]
        
    #     return subdivsion
       
    # def mesh_triangle(self,n:float,trigger:float):
       
    #    all_segments=[]
    #    subdivision=[self]
       

                   
               
                
    # def mesh_triangle(self,n:float):
    #     set_1=[]
    #     set_2=[]
    #     sets=[]
    #     all_polygons=[]
    #     new_points=[]
    #     segment_to_node={}
    #     linear_elements=self.linear_elements
        
        
    #     segment_to_node=self.mesh_segment(n)
       
    #     if segment_to_node[linear_elements[0]]!=[] and segment_to_node[linear_elements[1]]!=[] and segment_to_node[linear_elements[2]]!=[]:
    #         common_edge=linear_elements[0].line_intersection(linear_elements[1])
    #         nodes_1=segment_to_node[linear_elements[0]]
           
    #         nodes_2=segment_to_node[linear_elements[1]]
           
            
    #         best_point1=linear_elements[0].closest_middle(nodes_1)  
            
    #         best_point2=linear_elements[1].closest_middle(nodes_2) 
      
    #         for k in range(best_point1+1):
                
    #                 set_1.append(nodes_1[k])
    #         for  j in range(best_point2,len(nodes_2)):
                    
    #                 set_1.append(nodes_2[j])
                    
    #         sets.append(set_1)       
                    
    #         for k in range(best_point1,len(nodes_1)):
                
    #                 set_2.append(nodes_1[k])
    #         for  j in range(best_point2+1):
    #             if nodes_2[j] not in set_2:
    #                 set_2.append(nodes_2[j])
    #         sets.append(set_2)
           
    #         for set_0 in sets:
    #             if common_edge not in set_0:
    #                 new_polygon=vm.Polygon2D(set_0)
                   
    #                 all_polygons.extend(new_polygon.delaunay_triangulation())
    #             else : 
                
    #                 triangle=TriangularElement([nodes_1[best_point1],nodes_2[best_point2],common_edge])
                   
                 
                   
                   
    #                 all_polygons.extend(triangle.mesh_triangle(n))
                    
    #     else : 
    #             return [self.triangle_to_polygon()]
    #     print(len(all_polygons))
    #     return all_polygons
    
    # def complete_mesh(self,n:float):
      
    #      all_polygons = self.mesh_triangle(n)
    #      all_points=[p.points for p in all_polygons]
    #      for points in all_points:
    #          for k in range(len(points)-1):
    #              all_segments.append(vm.LineSegment2D(points[k],points[k+1]))
    #              for seg in all_segments:
                  
    #                  for polygon in all_polygons:
                      
    #                      if polygon.is_intersecting(seg)==False:
    #                          all_polygons.append(vm.Polygon2D([seg.point1,seg.point2,))
                 
    #    return all_polygons
         
           
           
    def mesh_triangle(self,n:float):
       segment_to_nodes=self.mesh_segments(n)
       # print(segment_to_nodes)

       
     
       linear_elements=self.linear_elements
       
       min_segment=[]
       interior_segments=[]
       interior_segment_nodes={}
       all_triangles=[]
       
       k=0
       nodes_0=[]
       nodes_1=[]
       if len(segment_to_nodes[linear_elements[1]])> len(segment_to_nodes[linear_elements[0]]) and len(segment_to_nodes[linear_elements[2]])>= len(segment_to_nodes[linear_elements[0]]) :
           
           nodes_0=segment_to_nodes[linear_elements[1]]
           nodes_1=segment_to_nodes[linear_elements[2]]
           min_segment.append(linear_elements[0])
           
       if len(segment_to_nodes[linear_elements[1]])< len(segment_to_nodes[linear_elements[0]]) and len(segment_to_nodes[linear_elements[2]])>= len(segment_to_nodes[linear_elements[0]]) :
           
           nodes_0=segment_to_nodes[linear_elements[0]]
           nodes_1=segment_to_nodes[linear_elements[2]]
           min_segment.append(linear_elements[1])
             
       if len(segment_to_nodes[linear_elements[0]])> len(segment_to_nodes[linear_elements[2]]) and len(segment_to_nodes[linear_elements[1]])>= len(segment_to_nodes[linear_elements[0]]) :
           
           nodes_0=segment_to_nodes[linear_elements[0]]
           nodes_1=segment_to_nodes[linear_elements[1]]
           min_segment.append(linear_elements[2])
       if len(segment_to_nodes[linear_elements[0]])> len(segment_to_nodes[linear_elements[1]]) and len(segment_to_nodes[linear_elements[1]])>= len(segment_to_nodes[linear_elements[2]]) :
           
           nodes_0=segment_to_nodes[linear_elements[0]]
           nodes_1=segment_to_nodes[linear_elements[1]]
           min_segment.append(linear_elements[2])
           
       if len(segment_to_nodes[linear_elements[0]])== len(segment_to_nodes[linear_elements[2]]) and len(segment_to_nodes[linear_elements[1]])==len(segment_to_nodes[linear_elements[0]]):
         
           nodes_0=segment_to_nodes[linear_elements[0]]
           nodes_1=segment_to_nodes[linear_elements[1]]
           min_segment.append(linear_elements[2])
       # print(nodes_1)
       # print(nodes_0)
       
       
       l0=min(len(nodes_0),len(nodes_1))
       l1=max(len(nodes_0),len(nodes_1))
      
       if len(nodes_0)>len(nodes_1):
           while k < l0-2:
              
               interior_segment=vm.LineSegment2D(nodes_0[k+1],nodes_1[len(nodes_1)-k-2])
              
               interior_segments.append(interior_segment)
               k=k+1
       else :
            while k < l0-2:
              
               interior_segment=vm.LineSegment2D(nodes_1[k+1],nodes_0[len(nodes_0)-k-2])
              
               interior_segments.append(interior_segment)
           
               k=k+1
    
       # if len(nodes_0)>len(nodes_1):
       #       for k in range(l0,len(nodes_0)-1):
       #           new_triangle=TriangularElement([nodes_0[k],nodes_1[len(nodes_1)],nodes_1[len(nodes_1)-1]])
               
       #           all_triangles.append(new_triangle)
          
              
       # if len(nodes_1)>len(nodes_0):
       #        for k in range(l0,len(nodes_1)-1):
       #          new_triangle=vm.TriangularElement([nodes_1[k],nodes_0[len(nodes_0)],nodes_0[len(nodes_0)-1]])
             
       #          all_triangles.append(new_triangle)
       
       
       
       for seg in interior_segments:
         

               interior_segment_nodes[seg]=seg.mesh_segment(n)[seg]
               
       if min_segment[0].point_distance(vm.Point2D(interior_segments[0].points[0])) < min_segment[0].point_distance(interior_segments[len(interior_segments)-1].points[0]):
           
           interior_segments.insert(0,min_segment[0])
           
           interior_segment_nodes[interior_segments[0]]=segment_to_nodes[interior_segments[0]]
       
       else :
           
           interior_segments.insert(len(interior_segments),min_segment[0])
           
       
           interior_segment_nodes[interior_segments[len(interior_segments)-1]]=segment_to_nodes[interior_segments[len(interior_segments)-1]]
       for k in range(len(interior_segments)):
           
            print(interior_segments[k].points )
          
       
       for k in range(len(interior_segments)-1):
           
           u=len(interior_segment_nodes[interior_segments[k]])
           v=len(interior_segment_nodes[interior_segments[k+1]])
           print(u)
           print(v)
           if interior_segment_nodes[interior_segments[k]][0].point_distance(interior_segment_nodes[interior_segments[k+1]][0])> interior_segment_nodes[interior_segments[k]][0].point_distance(interior_segment_nodes[interior_segments[k+1]][1]):
               interior_segment_nodes[interior_segments[k+1]].reverse()
               
           if u>2 and v>2:
               print('ok')
               for j in range(min(u,v)-1):
                       
                       
                       new_triangle_1=TriangularElement([interior_segment_nodes[interior_segments[k]][j+1],interior_segment_nodes[interior_segments[k]][j],interior_segment_nodes[interior_segments[k+1]][j]])
                       new_triangle_2=TriangularElement([interior_segment_nodes[interior_segments[k]][j+1],interior_segment_nodes[interior_segments[k+1]][j+1],interior_segment_nodes[interior_segments[k+1]][j]])
                       if new_triangle_1 not in all_triangles:
                           
                           all_triangles.append(new_triangle_1)
                           
                       if new_triangle_2 not in all_triangles:
                              # print(new_triangle_2.points)
                           
                              all_triangles.append(new_triangle_2)
               # for j in range(min(u,v),u-1):
               #             new_triangle=TriangularElement([interior_segment_nodes[interior_segments[k]][j],interior_segment_nodes[interior_segments[k]][j+1],interior_segment_nodes[interior_segments[k+1]][v-1]])
                       
               #             if new_triangle not in all_triangles:
                           
                #                all_triangles.append(new_triangle)       
           if v==2 and u>2:
                 for j in range(math.ceil(u/2)):
                   
                   new_triangle_1=TriangularElement([interior_segment_nodes[interior_segments[k]][j],interior_segment_nodes[interior_segments[k]][j+1],interior_segment_nodes[interior_segments[k+1]][0]])
                  
                   if new_triangle_1 not in all_triangles:
                     
                       all_triangles.append(new_triangle_1) 
                 for j in range((math.ceil(u/2)-1),u-1):
                   new_triangle_2=TriangularElement([interior_segment_nodes[interior_segments[k]][j],interior_segment_nodes[interior_segments[k]][j+1],interior_segment_nodes[interior_segments[k+1]][1]])
           
          
                   new_triangle_2=TriangularElement([interior_segment_nodes[interior_segments[k]][0],interior_segment_nodes[interior_segments[k]][1],interior_segment_nodes[interior_segments[k+1]][0]]) 
                   if new_triangle_2 not in all_triangles:
                           
                                  all_triangles.append(new_triangle_2)
                  # for j in range(u-1):
                  #    # new_triangle=TriangularElement([interior_segment_nodes[interior_segments[k]][j],interior_segment_nodes[interior_segments[k+1]][j+1],interior_segment_nodes[interior_segments[k+1]][0]])
                  #    new_triangle=TriangularElement([interior_segment_nodes[interior_segments[k]][j],interior_segment_nodes[interior_segments[k+1]][0],interior_segment_nodes[interior_segments[k+1]][1]])
                  #    if new_triangle not in all_triangles:
                           
                  #            all_triangles.append(new_triangle)
                         
                  # new_triangle_1=TriangularElement([interior_segment_nodes[interior_segments[k]][u-1],interior_segment_nodes[interior_segments[k]][u-2],interior_segment_nodes[interior_segments[k+1]][1]])
             
                  # if new_triangle_1 not in all_triangles:
                           
                  #               all_triangles.append(new_triangle_1)
                  
                  # new_triangle_2=TriangularElement([interior_segment_nodes[interior_segments[k]][0],interior_segment_nodes[interior_segments[k]][1],interior_segment_nodes[interior_segments[k+1]][0]])
               
                  # if new_triangle_2 not in all_triangles:
                           
                  #              all_triangles.append(new_triangle_2)
             
           if u==2 and v>2:
                for j in range(math.ceil(v/2)):
                   
                   new_triangle_1=TriangularElement([interior_segment_nodes[interior_segments[k+1]][j],interior_segment_nodes[interior_segments[k+1]][j+1],interior_segment_nodes[interior_segments[k]][0]])
                  
                   if new_triangle_1 not in all_triangles:
                     
                       all_triangles.append(new_triangle_1) 
                for j in range((math.ceil(v/2)-1),v-1):
                   new_triangle_2=TriangularElement([interior_segment_nodes[interior_segments[k+1]][j],interior_segment_nodes[interior_segments[k+1]][j+1],interior_segment_nodes[interior_segments[k]][1]])
           
          
                   new_triangle_2=TriangularElement([interior_segment_nodes[interior_segments[k+1]][0],interior_segment_nodes[interior_segments[k+1]][1],interior_segment_nodes[interior_segments[k]][0]]) 
                   if new_triangle_2 not in all_triangles:
                           
                                  all_triangles.append(new_triangle_2)
       return all_triangles  
    
      
    
          
           
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
                     
    
        return ax
    
class Mesher(DessiaObject):
    
    def __init__(self,contour:vm.Contour2D,triangles:List[TriangularElement],nodes_len:float):
        self.nodes_len=nodes_len
        self.contour=contour
        self.polygon=self.contour._get_polygon()
        self.triangles=triangles
        
                
        
        
    def  neighbour_edge(self,n:int,i:int,di:int):
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
    
       
    def polygon_to_triangles(self,polygons:List[vm.Polygon2D]):
        triangles=[]
        for polygon in polygons:
            triangles.append(TriangularElement(polygon.points))
        return triangles
    # def delaunay_triangulation(self,points:List[vm.Point2D]):
    #     new_points=[]
    #     delaunay_triangles=[]
        
    #     for point in points : 
    #         new_points.append([point[0],point[1]])
        
            
    #     delaunay=npy.array(new_points)  
        
    #     tri=Delaunay(delaunay)
        
       
    #     fig,ax=plt.subplots()
        
      
    #     for simplice in delaunay[tri.simplices]:
        
    #         triangle=TriangularElement([vm.Point2D(simplice[0]),vm.Point2D(simplice[1]),vm.Point2D(simplice[2])])
    #         delaunay_triangles.append(triangle)
    #         triangle.plot(ax=ax)
                
        
        
    #     return delaunay_triangles
         
             
  
        
        
        
    def plot_triangulised_contour(self):
        fig,ax=plt.subplots()
        patches=[]
        for triangle in self.triangles:
            triangle.plot(ax=ax)
        return ax
    
    
    # def assemble_mesh(self,triangles:List[TriangularElement]):
    #     all_triangles=[]
    #     all_triangles+=self.delaunay_triangulation(self.polygon.points)   
        
      
    #     for triangle in triangles :
                
    #             all_polygons=[]
    #             linear_elements=triangle.linear_elements
    #             segment_to_nodes=triangle.mesh_segment(self.nodes_len)
    #             nodes_1=segment_to_nodes[linear_elements[0]]
    #             nodes_2=segment_to_nodes[linear_elements[1]]
    #             nodes_3=nodes_2=segment_to_nodes[linear_elements[2]]
    #             nodes=list(product(nodes_1,nodes_2))
             
    #             all_polygons.append(vm.Polygon2D([nodes[0][0],nodes[1][0],nodes[0][1],nodes[1][1]]))
    #             length=len(all_polygons)
    #             for k in range(1,len(nodes)):
    #                 for j in range(2,len(nodes)-1):
                    
    #                           if abs(nodes[k][0][0]nodes[j][0][])<10e-10 and abs(nodes[k][1]!=nodes[j][1])<10e-10:
                                     
                                    
    #                                   i=0
                                     
    #                                   while i<length:
    #                                       print(vm.Polygon2D([nodes[k][0],nodes[j][0],nodes[k][1],nodes[j][1]]).is_intersecting(all_polygons[i]))
                                         
    #                                       if vm.Polygon2D([nodes[k][0],nodes[j][0],nodes[k][1],nodes[j][1]]).is_intersecting(all_polygons[i])==False:
                                           
                                           
    #                                         all_polygons.append(vm.Polygon2D([nodes[k][0],nodes[j][0],nodes[k][1],nodes[j][1]]))
    #                                         i+=1
                                            
                                           
    #                                       else : 
    #                                           i+=1
                                             
                                                
                                    
              
    #             print(all_polygons)
    #             for polygon in all_polygons:
    #                 polygon.points()
    #                 all_triangles+=self.delaunay_triangulation(polygon.points)
                        
                
                
                
    #     return all_triangles
                    
               
                      
    def assemble_mesh(self,triangles:List[TriangularElement]):
        
        all_triangles=[]
        all_triangles+=triangles
        for triangle in triangles:
            print('next')
            all_triangles+=triangle.mesh_triangle(self.nodes_len)
        
        return all_triangles
   
                
            
        
        
      
        
               
          
          
                                                
                                    
              
                
                        
                
                
                

                    
                          
                
          
                
     
                
                
                
                    
        
        
        
        
       
        
        
        # for triangle in triangles :
        #     all_triangles+=(triangle.mesh_triangle(self.nodes_len,trigger))
            
        # return all_triangles
            
        
  
    
   
                        
                        
                        
                      
      

                
        
                
            
    # def assemble_mesh(self,triangles:List[TriangularElement],trigger):
    #     memo=[]
    #     nodes=self.nodes_on_segments(triangles)[0]
    #     for triangle in triangles: 
    #         l0=triangle.linear_elements[0].length()/3
    #         if l0 > trigger and len(nodes) >= 3 :
               
    #             nodes=self.nodes_on_segments(triangles)[0]
                
             
    #             new_triangles=[]
               
    #             for node in nodes :
                   
    #                 u=self.closest_neightbours(triangles,node)
    #                 if u!=[]:
                       
    #                     n1=u[0]
    #                     n2=u[1]
                       
                    
    #                     new_triangle = TriangularElement([node,n1,n2])
    #                     new_triangles.append(new_triangle)
    #                     if new_triangle not in memo:
    #                         memo.append(new_triangle)
                       
                        
    #             memo.extend(self.assemble_mesh(triangles=new_triangles,trigger=trigger))
               
               
          
      
    #     return memo
     
          
      
    
        
               
           
           
       
       
            
        
        
    
  
        
                    
            
       
        
    
    
    
    