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
