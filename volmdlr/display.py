#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes to define mesh for display use. Display mesh do not require good aspect ratios on elements.
"""

import math
import warnings
from typing import List, Union, TypeVar

import numpy as np
from dessia_common.typings import JsonSerializable
from dessia_common.core import DessiaObject, PhysicalObject
from numpy.typing import NDArray

import volmdlr.edges


class Node2D(volmdlr.Point2D):
    # TODO: remove this class when no longer used
    """
    A node is a point with some hash capabilities for performance.
    """

    def __init__(self, x: float, y: float, name: str = ""):
        self.x = x
        self.y = y
        volmdlr.Point2D.__init__(self, x, y, name)

    def __hash__(self):
        return int(1e6 * (self.x + self.y))

    def __eq__(self, other_node: "Node2D"):
        if other_node.__class__.__name__ not in ["Vector2D", "Point2D", "Node2D"]:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-06) and math.isclose(self.y, other_node.y, abs_tol=1e-06)

    @classmethod
    def from_point(cls, point2d, name: str = ""):
        """
        Creates a Node 2D from a Point 2D.

        :param point2d: 2d point.:
        :param name: point's name.
        """
        return cls(point2d.x, point2d.y, name=name)


class Node3D(volmdlr.Point3D):
    # TODO: remove this class when no longer used
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

    def __eq__(self, other_node: "Node3D"):
        if other_node.__class__.__name__ not in ["Vector3D", "Point3D", "Node3D"]:
            return False
        return (
            math.isclose(self.x, other_node.x, abs_tol=1e-06)
            and math.isclose(self.y, other_node.y, abs_tol=1e-06)
            and math.isclose(self.z, other_node.z, abs_tol=1e-06)
        )

    @classmethod
    def from_point(cls, point3d):
        """
        Creates a Node 3D from a Point 3D.

        :param point3d: 3d point.
        """
        return cls(point3d.x, point3d.y, point3d.z)


class MeshMixin:
    """
    Mixin class for 2D and 3D meshes.

    This is an abstract class.
    """

    MeshType = TypeVar("MeshType", bound="MeshMixin")

    _standalone_in_db = True
    _non_serializable_attributes = ["vertices", "triangles"]

    # MANIPULATION
    def merge(
        self, other: "MeshType", mutualize_vertices: bool = False, mutualize_triangles: bool = False
    ) -> "MeshType":
        """
        Merge two meshes.

        :param other:
        :param mutualize_vertices:
        :param mutualize_triangles:
        :return:
        """
        if self.__class__.__name__ != other.__class__.__name__:
            raise ValueError("Meshes should have same dimension.")

        if len(self.vertices) == 0 or len(self.triangles) == 0:
            return other
        if len(other.vertices) == 0 or len(other.triangles) == 0:
            return self

        merged_vertices = np.concatenate((self.vertices, other.vertices))
        merged_triangles = np.concatenate((self.triangles, other.triangles + len(self.vertices)))

        mesh = self.__class__(merged_vertices, merged_triangles)

        if mutualize_vertices:
            mesh = mesh.mutualize_vertices()
        if mutualize_triangles:
            mesh = mesh.mutualize_triangles()

        return mesh

    def round_vertices(self, decimals: int = 9) -> "MeshType":
        """Round the mesh vertices to a given number of decimals."""

        rounded_vertices = np.round(self.vertices, decimals)

        return self.__class__(rounded_vertices, self.triangles, self.name)

    def remove_degenerate_triangles(self, tol: float = 0.0) -> "MeshType":
        """Remove degenerate triangles from the mesh."""
        # Get vertices for each corner of the triangles
        v0, v1, v2 = (
            self.vertices[self.triangles[:, 0]],
            self.vertices[self.triangles[:, 1]],
            self.vertices[self.triangles[:, 2]],
        )

        # Calculate the squared distance between each pair of vertices
        dist_sq_v0_v1 = np.sum((v0 - v1) ** 2, axis=1)
        dist_sq_v1_v2 = np.sum((v1 - v2) ** 2, axis=1)
        dist_sq_v0_v2 = np.sum((v0 - v2) ** 2, axis=1)

        # Find triangles where all three vertices are distinct (no zero distances)
        valid_triangles_mask = (dist_sq_v0_v1 > tol) & (dist_sq_v1_v2 > tol) & (dist_sq_v0_v2 > tol)

        # Filter out invalid triangles
        valid_triangles = self.triangles[valid_triangles_mask]

        # Create a new Mesh3D instance with non-flat triangles
        return self.__class__(self.vertices, valid_triangles, self.name)

    def mutualize_vertices(self) -> "MeshType":
        """Remove duplicated vertices and remap triangles."""

        unique_vertices, indices_map = np.unique(self.vertices, axis=0, return_inverse=True)
        remapped_triangles = indices_map[self.triangles]

        return self.__class__(unique_vertices, remapped_triangles, self.name)

    def mutualize_triangles(self) -> "MeshType":
        """Remove duplicated triangles from a mesh with unique vertices."""

        sorted_triangles = np.sort(self.triangles, axis=1)
        _, unique_triangle_indices = np.unique(sorted_triangles, axis=0, return_index=True)
        unique_triangles = self.triangles[unique_triangle_indices]

        return self.__class__(self.vertices, unique_triangles, self.name)

    def __add__(self, other: "MeshType") -> "MeshType":
        """
        Overload the "+" operator to merge two Mesh instances, without mutualization of vertices and triangles.

        :param other: Another Mesh instance to concatenate with this instance.
        :type other: MeshType

        :return: A new Mesh instance representing the merged shells.
        :rtype: MeshType
        """
        return self.merge(other, mutualize_vertices=False, mutualize_triangles=False)

    def __or__(self, other: "MeshType") -> "MeshType":
        """
        Overload the "|" operator to merge two Mesh instances, with mutualization of vertices and triangles.

        :param other: Another Mesh instance to concatenate with this instance.
        :type other: MeshType

        :return: A new Mesh instance representing the concatenated shells.
        :rtype: MeshType
        """
        return self.merge(other, mutualize_vertices=True, mutualize_triangles=True)

    @classmethod
    def merge_meshes(cls, meshes: List[Union["Mesh2D", "Mesh3D"]], name: str = ""):
        """
        Merge several meshes into one.
        """
        # TODO: refactor this method to "from_meshes"

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

        points = np.concatenate(points_list, axis=0)
        triangles = np.concatenate(triangles_list, axis=0)

        return cls(points, triangles, name=name)

    def merge_mesh(self, other_mesh):
        """
        Merge two meshes.

        :param other_mesh: other mesh.
        :return:
        """
        # TODO: remove this method and use "merge"
        result = self + other_mesh
        self.vertices = result.vertices
        self.triangles = result.triangles

    # CHECK
    def check_concistency(self):
        """Check mesh concistency."""

        n_points = len(self.vertices)

        for triangle in self.triangles:
            if max(triangle) >= n_points:
                return False
        return True

    # COMPUTATION
    def triangles_vertices(self):
        """
        Actual triangles of the mesh (points, not indexes).

        :return: Points of triangle vertices.
        :rtype: np.ndarray[float]
        """
        triangles = self.vertices.view(np.ndarray)[self.triangles]
        return triangles

    def triangles_cross_products(self):
        """
        Compute the cross products of edges for each triangle in the mesh.
        """
        vectors = np.diff(self.triangles_vertices(), axis=1)
        return np.cross(vectors[:, 0], vectors[:, 1])

    def plot(self, ax=None, numbering: bool = False):
        """Plot the mesh with Matplotlib."""
        for i_points, point in enumerate(self.vertices):
            ax = self._point_class(*point).plot(ax=ax)
            if numbering:
                ax.text(*point, f"node {i_points + 1}", ha="center", va="center")

        for vertex1, vertex2, vertex3 in self.triangles_vertices():
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

    # SERIALIZATION
    def to_dict(self, *args, **kwargs):
        """Overload of 'to_dict' for numpy usage."""

        dict_ = self.base_dict()
        dict_["vertices"] = self.vertices.tolist()
        dict_["triangles"] = self.triangles.tolist()

        return dict_

    @classmethod
    def dict_to_object(cls, dict_: JsonSerializable, *args, **kwargs) -> "MeshType":
        """Overload of 'dict_to_object' for numpy usage."""

        vertices = np.array(dict_["vertices"])
        triangles = np.array(dict_["triangles"])
        name = dict_["name"]

        return cls(vertices, triangles, name)

    # HASH AND EQUALITY
    def __hash__(self):
        """Computation of hash."""
        return hash((self.__class__.__name__, self.vertices.tobytes(), self.triangles.tobytes()))

    def __eq__(self, other):
        """Equality."""
        return hash(self) == hash(other)

    def _data_hash(self):
        """Computation of hash for Dessia platform usage."""
        return hash(self)

    def _data_eq(self, other_object) -> bool:
        """Equality for Dessia platform usage."""
        if other_object.__class__.__name__ != self.__class__.__name__:
            return False
        return self == other_object


class Mesh2D(MeshMixin, DessiaObject):
    """
    2D mesh.
    """

    _linesegment_class = volmdlr.edges.LineSegment2D
    _point_class = volmdlr.Point2D

    def __init__(self, vertices: NDArray[float], triangles: NDArray[int], name: str = ""):
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
        areas = np.sqrt((self.triangles_cross_products() ** 2)) / 2.0
        return areas.sum()


class Mesh3D(MeshMixin, PhysicalObject):
    """
    3D mesh.
    """

    _linesegment_class = volmdlr.edges.LineSegment3D
    _point_class = volmdlr.Point3D

    def __init__(self, vertices: NDArray[float], triangles: NDArray[int], name: str = ""):
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
        areas = np.sqrt((self.triangles_cross_products() ** 2).sum(axis=1)) / 2.0
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
        warnings.warn("Please use the Stl.from_display_mesh method instead")
