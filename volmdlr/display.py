#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes to define mesh for display use. Display mesh do not require good aspect ratios on elements.
"""

import math
import warnings
from typing import List, Tuple, Union, Dict, Any
import numpy as np
from numpy.typing import NDArray

import dessia_common.core as dc
from dessia_common.typings import JsonSerializable

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


class MeshMixin:
    """
    A DisplayMesh is a list of points linked by triangles.

    This is an abstract class for 2D & 3D.
    """

    # def __add__(self, other_mesh):
    #     """
    #     Defines how to add two meshes.
    #     """
    #     new_points = self.vertices[:]
    #     new_point_index = self.point_index.copy()
    #     i_points = len(new_points)
    #     for point in other_mesh.points:
    #         if point not in new_point_index:
    #             new_point_index[point] = i_points
    #             i_points += 1
    #             new_points.append(point)
    #
    #     new_triangles = self.triangles[:]
    #     for vertex1, vertex2, vertex3 in other_mesh.triangles:
    #         point1 = other_mesh.points[vertex1]
    #         point2 = other_mesh.points[vertex2]
    #         point3 = other_mesh.points[vertex3]
    #         new_triangles.append((new_point_index[point1],
    #                               new_point_index[point2],
    #                               new_point_index[point3]))
    #
    #     return self.__class__(new_points, new_triangles)

    def to_dict(self, *args, **kwargs):
        """Overload of 'to_dict' for performance."""
        dict_ = self.base_dict()

        dict_["vertices"] = self.vertices.tolist()
        dict_["triangles"] = self.triangles.tolist()

        return dict_

    @classmethod
    def dict_to_object(cls, dict_: JsonSerializable, force_generic: bool = False,
                       global_dict=None, pointers_memo: Dict[str, Any] = None,
                       path: str = "#", name: str = "") -> "Union[Mesh2D, Mesh3D]":
        """Overload of 'dict_to_object' for performance."""

        vertices = np.array(dict_["vertices"])
        triangles = np.array(dict_["triangles"])
        name = dict_["name"]

        display_triangle_shell = cls(vertices, triangles, name)

        return display_triangle_shell

    def concatenate(self, other: "Union[Mesh2D, Mesh3D]") -> "Union[Mesh2D, Mesh3D]":
        """
        Concatenates two Mesh instances into a single instance.

        This method merges the vertices and indices of both mesh. If the same vertex exists in both mesh,
        it is only included once in the merged shell to optimize memory usage. It also ensures that each face is
        represented uniquely by sorting the vertices of each triangle.

        :param other: Another Mesh instance to concatenate with this instance.
        :return: A new Mesh instance representing the concatenated shells.
        """
        if len(self.vertices) == 0 or len(self.triangles) == 0:
            return other
        if len(other.vertices) == 0 or len(other.triangles) == 0:
            return self

        # Merge and remove duplicate vertices
        merged_vertices = np.vstack((self.vertices, other.vertices))
        unique_vertices, indices_map = np.unique(merged_vertices, axis=0, return_inverse=True)

        # Adjust indices to account for duplicates and offset from concatenation
        self_indices_adjusted = self.triangles
        other_indices_adjusted = other.triangles + self.vertices.shape[0]

        # Re-map indices to unique vertices
        all_indices = np.vstack((self_indices_adjusted, other_indices_adjusted))
        final_indices = indices_map[all_indices]

        # Use np.unique to find unique subarrays
        _, unique_indices = np.unique(np.sort(final_indices, axis=1), axis=0, return_index=True)

        # Get the unique subarrays
        merged_indices = final_indices[unique_indices]

        # Create a new DisplayTriangleShell3D with merged data
        return self.__class__(unique_vertices, merged_indices, name=self.name + "+" + other.name)

    def __add__(self, other: "Union[Mesh2D, Mesh3D]") -> "Union[Mesh2D, Mesh3D]":
        """
        Overloads the + operator to concatenate two Mesh instances.

        :param other: Another Mesh instance to concatenate with this instance.
        :type other: Mesh

        :return: A new Mesh instance representing the concatenated shells.
        :rtype: Union[Mesh2D, Mesh3D]
        """
        return self.concatenate(other)

    def __hash__(self):
        return hash(
            (
                self.__class__.__name__,
                (tuple(self.triangles[0]), tuple(self.triangles[-1]), len(self.triangles)),
                (tuple(self.vertices[0]), tuple(self.vertices[-1]), len(self.vertices)),
            )
        )

    def __eq__(self, other):
        return hash(self) == hash(other)

    def _data_hash(self):
        return hash(
            (
                self.__class__.__name__,
                (tuple(self.triangles[0]), tuple(self.triangles[-1]), len(self.triangles)),
                (tuple(self.vertices[0]), tuple(self.vertices[-1]), len(self.vertices)),
            )
        )

    def _data_eq(self, other_object):
        if other_object.__class__.__name__ != self.__class__.__name__:
            return False
        return self._data_hash() == other_object._data_hash()


    @property
    def triangles_vertices(self):
        """
        Actual triangles of the mesh (points, not indexes)

        :return: Points of triangle vertices
        :rtype: (n, 3, d) float
        """
        # use of advanced indexing on our tracked arrays will
        # trigger a change flag which means the hash will have to be
        # recomputed. We can escape this check by viewing the array.
        triangles = self.vertices.view(np.ndarray)[self.triangles]

        return triangles

    @property
    def point_index(self):
        if self._point_index is None:
            self._point_index = {point: index for index, point in enumerate(self.vertices)}
        return self._point_index

    def check(self):
        npoints = len(self.vertices)
        for triangle in self.triangles:
            if max(triangle) >= npoints:
                return False
        return True

    def triangles_crosses(self):
        vectors = np.diff(self.triangles_vertices, axis=1)
        crosses = np.cross(vectors[:, 0], vectors[:, 1])
        return crosses

    @classmethod
    def merge_meshes(cls, meshes: List[Union['Mesh2D', 'Mesh3D']], name: str = ''):
        """
        Merge several meshes into one.
        """
        if len(meshes) == 1:
            return cls(meshes[0].vertices, meshes[0].triangles, name=name)

        points_list = []
        triangles_list = []
        i_points = 0

        for mesh in meshes:
            if not mesh:
                continue
            points_list.append(mesh.vertices)
            triangles_list.append(mesh.triangles + i_points)
            i_points += mesh.vertices.shape[0]
        if points_list:
            points = np.concatenate(points_list, axis=0)
            triangles = np.concatenate(triangles_list, axis=0)
            return cls(points, triangles, name=name)
        return None

    def merge_mesh(self, other_mesh):
        """
        Merge two meshes.

        :param other_mesh: other mesh.
        :return:
        """
        result = self + other_mesh
        self.vertices = result.vertices
        self.triangles = result.triangles

    def plot(self, ax=None, numbering=False):
        """Plots the mesh with Matplotlib."""
        for i_points, point in enumerate(self.vertices):
            ax = self._point_class(*point).plot(ax=ax)
            if numbering:
                ax.text(*point, f'node {i_points + 1}',
                        ha='center', va='center')

        for vertex1, vertex2, vertex3 in self.triangles_vertices:
            point1 = self._point_class(*vertex1)
            point2 = self._point_class(*vertex2)
            point3 = self._point_class(*vertex3)
            if not point1.is_close(point2):
                self._linesegment_class(point1, point2).plot(ax=ax)
            if not point2.is_close(point3):
                self._linesegment_class(point2, point3).plot(ax=ax)
            if not point1.is_close(point3):
                self._linesegment_class(point1, point3).plot(ax=ax)

        return ax


class Mesh2D(MeshMixin, dc.PhysicalObject):
    """
    A mesh for display purposes in 2D.

    """

    _linesegment_class = volmdlr.edges.LineSegment2D
    _point_class = volmdlr.Point2D

    def __init__(self, vertices: NDArray[float], triangles: NDArray[int], name: str = ''):
        self.vertices = vertices
        self.triangles = triangles
        # Avoiding calling dessia object init because its inefficiency
        # dc.DessiaObject.__init__(self, name=name)
        self.name = name
        self._point_index = None

    def area(self):
        """
        Return the area as the sum of areas of triangles.
        """
        areas = np.sqrt((self.triangles_crosses() ** 2)) / 2.0
        return areas.sum()


class Mesh3D(MeshMixin, dc.PhysicalObject):
    """
    A mesh for display purposes in 3D.

    """

    _linesegment_class = volmdlr.edges.LineSegment3D
    _point_class = volmdlr.Point3D

    def __init__(self, vertices: NDArray[float], triangles: NDArray[int], name: str = ''):
        self.vertices = vertices
        self.triangles = triangles
        # Avoiding calling dessia object init because its inefficiency
        # dc.DessiaObject.__init__(self, name=name)
        self.name = name
        self._point_index = None
        self._faces = None

    def area(self):
        """
        Return the area as the sum of areas of triangles.
        """
        areas = np.sqrt((self.triangles_crosses() ** 2).sum(axis=1)) / 2.0
        return areas.sum()

    def to_babylon(self):
        """
        Returns mesh in babylonjs format.

        https://doc.babylonjs.com/how_to/custom
        """
        babylon_mesh = {"positions": self.vertices.flatten().tolist(), "indices": self.triangles.flatten().tolist()}
        return babylon_mesh

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
        for vertex1, vertex2, vertex3 in self.triangles_vertices:
            try:
                point1 = volmdlr.Point3D(*vertex1)
                point2 = volmdlr.Point3D(*vertex2)
                point3 = volmdlr.Point3D(*vertex3)
            except TypeError:
                print(True)
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
