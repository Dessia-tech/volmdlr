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
import triangle
from volmdlr.core_compiled import Matrix33

from dessia_common import DessiaObject
from typing import TypeVar, List, Tuple


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
    def __init__(self, points:Tuple[vm.Point2D, vm.Point2D], interior_normal:vm.Vector2D):
        self.points = points
        self.interior_normal = interior_normal
        
        DessiaObject.__init__(self, name='')
        
    def __hash__(self):
        return self.points[0].__hash__() + self.points[1].__hash__()
        
    def __eq__(self, other_linear_element):
        if self.__class__ != other_linear_element.__class__:
            return False
        return (self.points[0] == other_linear_element.points[0] and self.points[1] == other_linear_element.points[1]) \
            or (self.points[0] == other_linear_element.points[1] and self.points[1] == other_linear_element.points[0])
            
    def length(self):
        return self.points[1].point_distance(self.points[0])
    
    def plot(self, ax=None, color='k', width=None, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        if width is None:
            width=1
        if plot_points:
            ax.plot([self.points[0][0], self.points[1][0]], [self.points[0][1], self.points[1][1]], color=color, marker='o', linewidth=width)
        else:
            ax.plot([self.points[0][0], self.points[1][0]], [self.points[0][1], self.points[1][1]], color=color, linewidth=width)
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
        linear_element_1 = LinearElement([self.points[0], self.points[1]], normal1)
        linear_element_2 = LinearElement([self.points[1], self.points[2]], normal2)
        linear_element_3 = LinearElement([self.points[2], self.points[0]], normal3)
        return [linear_element_1, linear_element_2, linear_element_3]
    
    def _form_functions(self):
        a = Matrix33(1, self.points[0][0], self.points[0][1],
                     1, self.points[1][0], self.points[1][1],
                     1, self.points[2][0], self.points[2][1])
        try :
            inv_a = a.inverse()
        except ValueError:
            self.plot()
            print(self._area())
            raise FlatElementError('form function bug')
        x1 = inv_a.vector_multiplication(vm.X3D)
        x2 = inv_a.vector_multiplication(vm.Y3D)
        x3 = inv_a.vector_multiplication(vm.Z3D)
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
    def __init__(self, elements:List[TriangularElement], mu_total:float, name:str):
        self.elements = elements
        self.mu_total = mu_total
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



    