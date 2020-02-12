#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:05:56 2020

@author: ringhausen
"""

import matplotlib.pyplot as plt
import volmdlr.core as vm

class TriangularElement:
    def __init__(self, points):
        if len(points) != 3:
            raise ValueError
            
        self.points = points
        
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
                
    def plot(self, ax=None, color='k', width=None, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        
        for p1, p2 in zip(self.points, self.points[1:]+[self.points[0]]):
            if width is None:
                width=1
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, marker='o', linewidth=width)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color, linewidth=width)
            
        return ax
    
class ElementsGroup:
    def __init__(self, elements, name):
        self.elements = elements 
        self.name = name
        
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
                
    def plot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
            
        for element in self.elements:
            element.plot(ax=ax)
            
        return ax
        
class Mesh:
    def __init__(self, elements_groups):
        self.elements_groups = elements_groups

    