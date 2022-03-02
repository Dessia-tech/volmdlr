#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module containing grid and relative objects
"""

import volmdlr as vm
import volmdlr.wires
import numpy as npy
from dessia_common import DessiaObject



class Grid2D(DessiaObject):
    
    def __init__(self, list_points: list[list[volmdlr.vm.Point2D]],
                 direction: list[str] = ['+x','+y'],
                 name: str = ''):
        self.list_points = list_points
        direction = direction
        DessiaObject.__init__(self, name=name)

    @classmethod
    def from_properties(x_limits, y_limits, points_nbr, direction):
        """
        Define Grid2d based on the given properties

        Parameters
        ----------
        x_limits : tuple(float, float)
            (x_min, x_max)
        y_limits : tuple(float, float)
            (y_min, y_max)
        points_nbr : tuple(int, int)
            number of points along x-axis and y-axis
        direction : list[str]
            it is used to order the generated points

        Returns
        -------
        Grid2d

        """

        xmin, xmax = x_limits
        ymin, ymax = y_limits
        points_x, points_y = points_nbr

        # points_2d = []
        grid2d = []
        points = []

        if direction == ['+x','+y']:
            x = npy.linspace(xmin, xmax, points_x)
            y = npy.linspace(ymin, ymax, points_y)

            for yi in y:
                for xi in x:
                    # points_2d.append(vm.Point2D(xi, yi))
                    points.append(vm.Point2D(xi, yi))
                
                grid2d.append(points)
                points = []
        
        elif direction == ['-x','+y']:
            x = npy.linspace(xmax, xmin, points_x) 
            y = npy.linspace(ymin, ymax, points_y) 

            for yi in y:
                for xi in x:
                    # points_2d.append(vm.Point2D(xi, yi))
                    points.append(vm.Point2D(xi, yi))

                grid2d.append(points)
                points = []
                
        elif direction == ['+y','+x']:
            x = npy.linspace(xmin, xmax, points_x) 
            y = npy.linspace(ymin, ymax, points_y) 
            
            for xi in x:
                for yi in y:
                    # points_2d.append(vm.Point2D(xi, yi))
                    points.append(vm.Point2D(xi, yi))
        
                grid2d.append(points)
                points = []

        elif direction == ['-y','+x']:
            x = npy.linspace(xmin, xmax, points_x) 
            y = npy.linspace(ymax, ymin, points_y) 

            for xi in x:
                for yi in y:
                    # points_2d.append(vm.Point2D(xi, yi))
                    points.append(vm.Point2D(xi, yi))
                    
                grid2d.append(points)
                points = []
           
        elif direction == ['+x','-y']:
            x = npy.linspace(xmin, xmax, points_x) 
            y = npy.linspace(ymax, ymin, points_y) 

            for yi in y:
                for xi in x:
                    # points_2d.append(vm.Point2D(xi, yi))
                    points.append(vm.Point2D(xi, yi))
                
                grid2d.append(points)
                points = []

        elif direction == ['-x','-y']:
            x = npy.linspace(xmax, xmin, points_x) 
            y = npy.linspace(ymax, ymin, points_y) 

            for yi in y:
                for xi in x:
                    # points_2d.append(vm.Point2D(xi, yi)) 
                    points.append(vm.Point2D(xi, yi))
                    
                grid2d.append(points)
                points = []    
                    
        elif direction == ['+y','-x']:
            x = npy.linspace(xmax, xmin, points_x) 
            y = npy.linspace(ymin, ymax, points_y) 

            for xi in x:
                for yi in y:
                    # points_2d.append(vm.Point2D(xi, yi)) 
                    points.append(vm.Point2D(xi, yi))
                    
                grid2d.append(points)
                points = []

        elif direction == ['-y','-x']:
            x = npy.linspace(xmax, xmin, points_x) 
            y = npy.linspace(ymax, ymin, points_y) 

            for xi in x:
                for yi in y:
                    # points_2d.append(vm.Point2D(xi, yi))
                    points.append(vm.Point2D(xi, yi))
                    
                grid2d.append(points)
                points = []

        return cls(list_points=grid2d, direction=direction)

    # @property
    # def x_limits(self):
    #     """
    #     find the limits (min, max) of points along x_direction_axis
    #     """

    #     x_limits = (self.list_points[0][0].x, self.list_points[-1][-1].x)
    #     x_min, x_max = min(x_limits), max(x_limits)

    #     return (x_min, x_max)

    # @property
    # def y_limits(self):
    #     """
    #     find the limits (min, max) of points along y_direction_axis
    #     """

    #     y_limits = (self.list_points[0][0].y, self.list_points[-1][-1].y)
    #     y_min, y_max = min(y_limits), max(y_limits)

    #     return (y_min, y_max)

    @property
    def limits_xy(self):
        """
        find the limits (min, max) of points along x & y direction_axis
        """

        x_limits = (self.list_points[0][0].x, self.list_points[-1][-1].x)
        x_min, x_max = min(x_limits), max(x_limits)
        y_limits = (self.list_points[0][0].y, self.list_points[-1][-1].y)
        y_min, y_max = min(y_limits), max(y_limits)

        return ((x_min, x_max), (y_min, y_max))

    # @property
    # def points_x(self):
    #     """
    #     find how many points there are along x_direction_axis
    #     """

    #     index = self.find_direction_index(direction_axis = 'x')        
    #     if index == 0:
    #         points_x = len(self.list_points[0])
    #     else:
    #         points_x = len(self.list_points)

    #     return points_x

    # @property        
    # def points_y(self):
    #     """
    #     find how many points there are along y_direction_axis 
    #     """

    #     index = self.find_direction_index(direction_axis = 'y')      
    #     if index == 0:
    #         points_y = len(self.list_points[0])
    #     else:
    #         points_y = len(self.list_points)

    #     return points_y

    @property
    def points_xy(self):
        """
        find how many points there are along x & y direction_axis
        """
        
        index = self.find_direction_index(direction_axis = 'x')        
        if index == 0:
            points_x = len(self.list_points[0])
            points_y = len(self.list_points)
        else:
            points_x = len(self.list_points)
            points_y = len(self.list_points[0])
        
        return (points_x, points_y)

    def find_direction_index(self, direction_axis: str):
        """
        find the index of a given direction_axis

        Parameters
        ----------
        direction_axis : str
            'x' OR 'y'

        Returns
        -------
        index : int

        """

        try:
            index = self.direction.index('+'+direction_axis)
        except ValueError:
            index = self.direction.index('-'+direction_axis)

        return index

