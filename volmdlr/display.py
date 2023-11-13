#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes to define mesh for display use. Display mesh do not require good aspect ratios on elements.
"""

import math
import warnings
from typing import List, Tuple

import dessia_common.core as dc

import volmdlr.edges


class Node2D(volmdlr.Point2D):
    """
    A node is a point with some hash capabilities for performance.
    """

    def __init__(self, x: float, y: float, name: str = ""):
        self.x = x
        self.y = y
        volmdlr.Point2D.__init__(self, x, y, name)

    def __hash__(self):
        return int(1e6 * (self.x + self.y))

    def __eq__(self, other_node: 'Node2D'):
        if other_node.__class__.__name__ not in ['Vector2D', 'Point2D',
                                                 'Node2D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-06) \
            and math.isclose(self.y, other_node.y, abs_tol=1e-06)

    @classmethod
    def from_point(cls, point2d, name: str = ''):
        """
        Creates a Node 2D from a Point 2D.

        :param point2d: 2d point.:
        :param name: point's name.
        """
        return cls(point2d.x, point2d.y, name=name)


class Node3D(volmdlr.Point3D):
    """
    A node is a point with some hash capabilities for performance.
    """

    def __init__(self, x: float, y: float, z: float, name: str = ""):
        self.x = x
        self.y = y
        self.z = z
        volmdlr.Point3D.__init__(self, x, y, z, name)

    def __hash__(self):
        return int(1e6 * (self.x + self.y + self.z))

    def __eq__(self, other_node: 'Node3D'):
        if other_node.__class__.__name__ not in ['Vector3D', 'Point3D',
                                                 'Node3D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-06) \
            and math.isclose(self.y, other_node.y, abs_tol=1e-06) \
            and math.isclose(self.z, other_node.z, abs_tol=1e-06)

    @classmethod
    def from_point(cls, point3d):
        """
        Creates a Node 3D from a Point 3D.

        :param point3d: 3d point.
        """
        return cls(point3d.x, point3d.y, point3d.z)


class DisplayMesh(dc.DessiaObject):
    """
    A DisplayMesh is a list of points linked by triangles.

    This is an abstract class for 2D & 3D.
    """
    _linesegment_class = volmdlr.edges.LineSegment

    def __init__(self, points, triangles, name: str = ''):

        self.points = points
        self.triangles = triangles
        # Avoiding calling dessia object init because its inefficiency
        # dc.DessiaObject.__init__(self, name=name)
        self.name = name
        self._point_index = None

    def check(self):
        npoints = len(self.points)
        for triangle in self.triangles:
            if max(triangle) >= npoints:
                return False
        return True

    @property
    def point_index(self):
        if self._point_index is None:
            self._point_index = {point: index for index, point in enumerate(self.points)}
        return self._point_index

    @classmethod
    def merge_meshes(cls, meshes: List['DisplayMesh'], name: str = ''):
        """
        Merge several meshes into one.
        """
        # Collect points
        i_points = 0
        point_index = {}
        points = []
        if len(meshes) == 1:
            return cls(meshes[0].points, meshes[0].triangles, name=name)
        for mesh in meshes:
            if not mesh:
                continue
            for point in mesh.points:
                if point not in point_index:
                    point_index[point] = i_points
                    i_points += 1
                    points.append(point)

        triangles = []
        for mesh in meshes:
            if not mesh:
                continue
            for vertex1, vertex2, vertex3 in mesh.triangles:
                point1 = mesh.points[vertex1]
                point2 = mesh.points[vertex2]
                point3 = mesh.points[vertex3]
                triangles.append((point_index[point1],
                                  point_index[point2],
                                  point_index[point3]))
        return cls(points, triangles, name=name)

    def merge_mesh(self, other_mesh):
        """
        Merge two meshes.

        :param other_mesh: other mesh.
        :return:
        """
        i_points = len(self.points)
        for point in other_mesh.points:
            if point not in self.point_index:
                self.point_index[point] = i_points
                i_points += 1
                self.points.append(point)

        for vertex1, vertex2, vertex3 in other_mesh.triangles:
            point1 = other_mesh.points[vertex1]
            point2 = other_mesh.points[vertex2]
            point3 = other_mesh.points[vertex3]
            self.triangles.append((self._point_index[point1],
                                   self._point_index[point2],
                                   self._point_index[point3]))

    def __add__(self, other_mesh):
        """
        Defines how to add two meshes.
        """
        new_points = self.points[:]
        new_point_index = self.point_index.copy()
        i_points = len(new_points)
        for point in other_mesh.points:
            if point not in new_point_index:
                new_point_index[point] = i_points
                i_points += 1
                new_points.append(point)

        new_triangles = self.triangles[:]
        for vertex1, vertex2, vertex3 in other_mesh.triangles:
            point1 = other_mesh.points[vertex1]
            point2 = other_mesh.points[vertex2]
            point3 = other_mesh.points[vertex3]
            new_triangles.append((new_point_index[point1],
                                  new_point_index[point2],
                                  new_point_index[point3]))

        return self.__class__(new_points, new_triangles)

    def plot(self, ax=None, numbering=False):
        """Plots the mesh with Matplotlib."""
        for i_points, point in enumerate(self.points):
            ax = point.plot(ax=ax)
            if numbering:
                ax.text(*point, f'node {i_points + 1}',
                        ha='center', va='center')

        for vertex1, vertex2, vertex3 in self.triangles:
            point1 = self.points[vertex1]
            point2 = self.points[vertex2]
            point3 = self.points[vertex3]
            if not point1.is_close(point2):
                self._linesegment_class(point1, point2).plot(ax=ax)
            if not point2.is_close(point3):
                self._linesegment_class(point2, point3).plot(ax=ax)
            if not point1.is_close(point3):
                self._linesegment_class(point1, point3).plot(ax=ax)

        return ax


class DisplayMesh2D(DisplayMesh):
    """
    A mesh for display purposes in 2D.

    """

    _linesegment_class = volmdlr.edges.LineSegment2D
    _point_class = volmdlr.Point2D

    def __init__(self, points: List[volmdlr.Point2D],
                 triangles: List[Tuple[int, int, int]],
                 name: str = ''):
        DisplayMesh.__init__(self, points, triangles, name=name)

    def area(self):
        """
        Return the area as the sum of areas of triangles.
        """
        area = 0.
        for (vertex1, vertex2, vertex3) in self.triangles:
            point1 = self.points[vertex1]
            point2 = self.points[vertex2]
            point3 = self.points[vertex3]
            area += 0.5 * abs((point2 - point1).cross(point3 - point1))
        return area


class DisplayMesh3D(DisplayMesh):
    """
    A mesh for display purposes in 3D.

    """

    _linesegment_class = volmdlr.edges.LineSegment3D
    _point_class = volmdlr.Point3D

    def __init__(self, points: List[volmdlr.Point3D],
                 triangles: List[Tuple[int, int, int]], name=''):
        self._faces = None
        DisplayMesh.__init__(self, points, triangles, name=name)

    def to_babylon(self):
        """
        Returns mesh in babylonjs format.

        https://doc.babylonjs.com/how_to/custom
        """
        positions = []
        for point in self.points:
            # positions.extend(list(round(p, 6)))
            # Not using round for performance
            positions.extend([int(1e6 * point.x) / 1e6, int(1e6 * point.y) / 1e6, int(1e6 * point.z) / 1e6])

        flatten_indices = []
        for vertex in self.triangles:
            flatten_indices.extend(vertex)
        return positions, flatten_indices

    @property
    def faces(self):
        """
        Gets the mesh's triangular faces.

        """
        if not self._faces:
            self._faces = self.triangular_faces()
        return self._faces

    def triangular_faces(self):
        """
        Calculates the mesh's triangular faces.

        """
        triangular_faces = []
        for (vertex1, vertex2, vertex3) in self.triangles:
            point1 = self.points[vertex1]
            point2 = self.points[vertex2]
            point3 = self.points[vertex3]
            if not point1.is_close(point2) and not point2.is_close(point3) and not point1.is_close(point3):
                face = volmdlr.faces.Triangle3D(point1, point2, point3)
                if face.area() >= 1e-11:
                    triangular_faces.append(face)
        return triangular_faces

    def to_stl(self):
        """
        Exports to STL.

        """
        warnings.warn('Please use the Stl.from_display_mesh method instead')
