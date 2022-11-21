#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module containing grid and relative objects
"""

from typing import List
import numpy as npy
from dessia_common import DessiaObject
import volmdlr
import volmdlr.wires


class Grid2D(DessiaObject):

    def __init__(self, lists_points: List[List[volmdlr.Point2D]],
                 direction: List[str],
                 name: str = ''):
        self.lists_points = lists_points
        self.direction = direction
        DessiaObject.__init__(self, name=name)

    def displacement_compared_to(self, initial_grid2d):
        """
        compute the deformation/displacement (dx,dy) of a grid2d based on an another grid2d

        Parameters
        ----------
        grid2d : Grid2D

        """

        points_2d = initial_grid2d.points
        points_2d_deformed = self.points

        displacement = npy.ones(shape=(len(points_2d), 2))  # Grid2D points displacement
        for i in range(0, len(displacement)):
            displacement[i][0] = points_2d_deformed[i][0] - points_2d[i][0]
            displacement[i][1] = points_2d_deformed[i][1] - points_2d[i][1]

        return displacement

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
            index = self.direction.index('+' + direction_axis)
        except ValueError:
            index = self.direction.index('-' + direction_axis)

        return index

    @classmethod
    def from_points(cls, points, points_dim_1, direction):
        """
        define a Grid2D given a list of points, number of points along the 1st dimension, and a direction

        Parameters
        ----------
        points : list[volmdlr.Point2D]
        points_dim1 : int
        direction : List[str]

        """

        lists_points = [points[i:i + points_dim_1] for i in range(0, len(points), points_dim_1)]

        return cls(lists_points, direction)

    @classmethod
    def from_properties(cls, x_limits, y_limits, points_nbr,
                        direction=None):
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
        if direction is None:
            direction = ['+x', '+y']

        xmin, xmax = x_limits
        ymin, ymax = y_limits
        points_x, points_y = points_nbr

        # points_2d = []
        grid2d = []
        points = []

        if direction == ['+x', '+y']:
            x = npy.linspace(xmin, xmax, points_x)
            y = npy.linspace(ymin, ymax, points_y)

            for yi in y:
                for xi in x:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction == ['-x', '+y']:
            x = npy.linspace(xmax, xmin, points_x)
            y = npy.linspace(ymin, ymax, points_y)

            for yi in y:
                for xi in x:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction == ['+y', '+x']:
            x = npy.linspace(xmin, xmax, points_x)
            y = npy.linspace(ymin, ymax, points_y)

            for xi in x:
                for yi in y:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction == ['-y', '+x']:
            x = npy.linspace(xmin, xmax, points_x)
            y = npy.linspace(ymax, ymin, points_y)

            for xi in x:
                for yi in y:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction == ['+x', '-y']:
            x = npy.linspace(xmin, xmax, points_x)
            y = npy.linspace(ymax, ymin, points_y)

            for yi in y:
                for xi in x:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction == ['-x', '-y']:
            x = npy.linspace(xmax, xmin, points_x)
            y = npy.linspace(ymax, ymin, points_y)

            for yi in y:
                for xi in x:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction == ['+y', '-x']:
            x = npy.linspace(xmax, xmin, points_x)
            y = npy.linspace(ymin, ymax, points_y)

            for xi in x:
                for yi in y:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction == ['-y', '-x']:
            x = npy.linspace(xmax, xmin, points_x)
            y = npy.linspace(ymax, ymin, points_y)

            for xi in x:
                for yi in y:
                    # points_2d.append(volmdlr.Point2D(xi, yi))
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        return cls(lists_points=grid2d, direction=direction)

    def grid_pattern(self):
        """
        define a list of quadrilateral polygons defined based on Grid2d points

        Returns
        -------
        quadrilateral_polygons : list[volmdlr.wires.ClosedPolygon2D]

        """

        quadrilateral_polygons = []
        length_1, length_2 = len(self.lists_points[0]), len(self.lists_points)
        for i in range(0, length_1 - 1):
            for j in range(0, length_2 - 1):
                quadrilateral_polygons.append(volmdlr.wires.ClosedPolygon2D(
                    [self.lists_points[i][j],
                     self.lists_points[i + 1][j],
                     self.lists_points[i + 1][j + 1],
                     self.lists_points[i][j + 1]]))

        return quadrilateral_polygons

    @property
    def limits_xy(self):
        """
        find the limits (min, max) of points along x & y direction_axis
        """

        x_limits = (self.lists_points[0][0].x, self.lists_points[-1][-1].x)
        x_min, x_max = min(x_limits), max(x_limits)
        y_limits = (self.lists_points[0][0].y, self.lists_points[-1][-1].y)
        y_min, y_max = min(y_limits), max(y_limits)

        return ((x_min, x_max), (y_min, y_max))

    # @property
    # def x_limits(self):
    #     """
    #     find the limits (min, max) of points along x_direction_axis
    #     """

    #     x_limits = (self.lists_points[0][0].x, self.lists_points[-1][-1].x)
    #     x_min, x_max = min(x_limits), max(x_limits)

    #     return (x_min, x_max)

    # @property
    # def y_limits(self):
    #     """
    #     find the limits (min, max) of points along y_direction_axis
    #     """

    #     y_limits = (self.lists_points[0][0].y, self.lists_points[-1][-1].y)
    #     y_min, y_max = min(y_limits), max(y_limits)

    #     return (y_min, y_max)

    @property
    def points(self):
        """
        return all the points in lists_points in just one list

        Returns
        -------
        points: list[volmdlr.Point2D]

        """

        points = []
        for list_point in self.lists_points:
            points.extend(list_point)
        return points

    @property
    def points_xy(self):
        """
        find how many points there are along x & y direction_axis
        """

        index = self.find_direction_index(direction_axis='x')
        if index == 0:
            points_x = len(self.lists_points[0])
            points_y = len(self.lists_points)
        else:
            points_x = len(self.lists_points)
            points_y = len(self.lists_points[0])

        return (points_x, points_y)

    # @property
    # def points_x(self):
    #     """
    #     find how many points there are along x_direction_axis
    #     """

    #     index = self.find_direction_index(direction_axis = 'x')
    #     if index == 0:
    #         points_x = len(self.lists_points[0])
    #     else:
    #         points_x = len(self.lists_points)

    #     return points_x

    # @property
    # def points_y(self):
    #     """
    #     find how many points there are along y_direction_axis
    #     """

    #     index = self.find_direction_index(direction_axis = 'y')
    #     if index == 0:
    #         points_y = len(self.lists_points[0])
    #     else:
    #         points_y = len(self.lists_points)

    #     return points_y
