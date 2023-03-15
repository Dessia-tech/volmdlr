#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module containing grid and relative objects.

"""

from typing import List

import numpy as npy

from dessia_common.core import DessiaObject  # isort: skip

import volmdlr
import volmdlr.wires


class Grid2D(DessiaObject):
    """
    A class defined with a list of points and characterized with a chosen direction.

    :param lists_points: A list of a list of points
    :type lists_points: List[List[volmdlr.Point2D]]
    :param direction: A direction
    :type direction: List[str]
    """

    def __init__(self, lists_points: List[List[volmdlr.Point2D]],
                 direction: List[str],
                 name: str = ''):

        self.lists_points = lists_points
        self.direction = direction
        DessiaObject.__init__(self, name=name)

    def displacement_compared_to(self, initial_grid2d):
        """
        Computes the deformation/displacement (dx,dy) of a grid2d based on another grid2d.

        :param initial_grid2d: A 2 dimensional grid
        :type initial_grid2d: :class:`volmdlr.grid.Grid2D`
        :return: The displacement of the 2 dimensional grid
        :rtype:
        """
        points_2d = initial_grid2d.points
        points_2d_deformed = self.points

        # Grid2D points displacement
        displacement = npy.ones(shape=(len(points_2d), 2))
        for i, displacement_i in enumerate(displacement):
            displacement_i[0] = points_2d_deformed[i][0] - points_2d[i][0]
            displacement_i[1] = points_2d_deformed[i][1] - points_2d[i][1]

        return displacement

    def find_direction_index(self, direction_axis: str):
        """
        Finds the index of a given direction_axis.

        :param direction_axis: 'x' or 'y'
        :type direction_axis: str
        :return: The direction index
        :rtype: int
        """

        try:
            index = self.direction.index('+' + direction_axis)
        except ValueError:
            index = self.direction.index('-' + direction_axis)

        return index

    @classmethod
    def from_points(cls, points, points_dim_1, direction):
        """
        Defines a Grid2D given a list of points, number of points along the 1st dimension, and a direction.

        :param points:
        :type points: List[:class:`volmdlr.Point2D`]
        :param points_dim_1:
        :type points_dim_1: int
        :param direction:
        :type direction: List[str]
        :return:
        :rtype:
        """

        lists_points = [points[i:i + points_dim_1] for i in range(0, len(points), points_dim_1)]

        return cls(lists_points, direction)

    @classmethod
    def from_properties(cls, x_limits, y_limits, points_nbr,
                        direction=None):
        """
        Defines Grid2d based on the given properties.

        :param x_limits: x_min and x_max
        :type x_limits: Tuple[float, float]
        :param y_limits: y_min and y_max
        :type y_limits: Tuple[float, float]
        :param points_nbr: Number of points along the x-axis and the y-axis
        :type points_nbr: Tuple[int, int]
        :param direction: Used for ordering the generated points
        :type direction: List[str]
        :return: The 2 dimensional grid
        :rtype: :class:`volmdlr.grid.Grid2D`
        """
        if direction is None:
            direction = ['+x', '+y']

        xmin, xmax = x_limits
        ymin, ymax = y_limits
        points_x, points_y = points_nbr

        directions_properties = {
            ('+x', '+y'): (xmin, xmax, ymin, ymax),
            ('-x', '+y'): (xmax, xmin, ymin, ymax),
            ('+y', '+x'): (xmin, xmax, ymin, ymax),
            ('-y', '+x'): (xmin, xmax, ymax, ymin),
            ('+x', '-y'): (xmin, xmax, ymax, ymin),
            ('-x', '-y'): (xmax, xmin, ymax, ymin),
            ('+y', '-x'): (xmax, xmin, ymin, ymax),
            ('-y', '-x'): (xmax, xmin, ymax, ymin)
        }

        grid2d = []
        points = []

        xmin, xmax, ymin, ymax = directions_properties[tuple(direction)]
        x = npy.linspace(xmin, xmax, points_x)
        y = npy.linspace(ymin, ymax, points_y)

        if direction in [['+x', '+y'], ['-x', '+y'], ['+x', '-y'], ['-x', '-y']]:
            for yi in y:
                for xi in x:
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        elif direction in [['+y', '+x'], ['-y', '+x'], ['+y', '-x'], ['-y', '-x']]:
            for xi in x:
                for yi in y:
                    points.append(volmdlr.Point2D(xi, yi))

                grid2d.append(points)
                points = []

        return cls(lists_points=grid2d, direction=direction)

    def grid_pattern(self):
        """
        Defines a list of quadrilateral polygons defined based on Grid2d points.

        :return: The list of quadrilateral polygons
        :rtype: List[volmdlr.wires.ClosedPolygon2D]
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
        Finds the limits (min, max) of points along x & y direction_axis.

        :return:
        :rtype:
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
        Returns all the points in lists_points in just one list.

        :return: The flattened list of points
        :rtype: List[:class:`volmdlr.Point2D`]
        """

        points = []
        for list_point in self.lists_points:
            points.extend(list_point)
        return points

    @property
    def points_xy(self):
        """
        Finds how many points there are along x & y direction_axis.

        :return: Two counts, one for the x direction_axis and another one for
            the y direction_axis
        :rtype: Tuple[int, int]
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
