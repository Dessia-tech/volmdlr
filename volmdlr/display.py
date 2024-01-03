#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes to define mesh for display use. Display mesh do not require good aspect ratios on elements.
"""

import math
import warnings
from typing import List, TypeVar, Union

import numpy as np
from dessia_common.core import DessiaObject, PhysicalObject
from dessia_common.typings import JsonSerializable
from dessia_common.serialization import BinaryFile
from numpy.typing import NDArray

import trimesh
from trimesh import Trimesh

import volmdlr.edges


# TODO: make this module "mesh" as it is not useful only for display


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
    def resize(self, scale_factor: float) -> "MeshType":
        return self.__class__(self.vertices * scale_factor, self.triangles, self.name)

    def merge(self, other: "MeshType", merge_vertices: bool = False, merge_triangles: bool = False) -> "MeshType":
        """
        Merge two meshes.

        :param other:
        :param merge_vertices:
        :param merge_triangles:
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

        mesh = self.__class__(merged_vertices, merged_triangles, self.name)

        if merge_vertices:
            mesh = mesh.merge_vertices()
        if merge_triangles:
            mesh = mesh.merge_triangles()

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

    def merge_vertices(self) -> "MeshType":
        """Merge duplicated vertices and remap triangles."""

        unique_vertices, indices_map = np.unique(self.vertices, axis=0, return_inverse=True)
        remapped_triangles = indices_map[self.triangles]

        return self.__class__(unique_vertices, remapped_triangles, self.name)

    def unmerge_vertices(self) -> "MeshType":
        """Unmerge shared vertices between triangles."""

        unmerged_vertices = self.vertices[self.triangles.ravel()]
        unmerged_triangles = np.arange(len(self.triangles) * 3).reshape(-1, 3)

        return self.__class__(unmerged_vertices, unmerged_triangles, self.name)

    def merge_triangles(self) -> "MeshType":
        """Merge duplicated triangles."""

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
        return self.merge(other, merge_vertices=False, merge_triangles=False)

    def __or__(self, other: "MeshType") -> "MeshType":
        """
        Overload the "|" operator to merge two Mesh instances, with mutualization of vertices and triangles.

        :param other: Another Mesh instance to concatenate with this instance.
        :type other: MeshType

        :return: A new Mesh instance representing the concatenated shells.
        :rtype: MeshType
        """
        return self.merge(other, merge_vertices=True, merge_triangles=True)

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

        # Plot vertices
        for i_point, point in enumerate(self.vertices):
            ax = self._point_class(*point).plot(ax=ax)
            if numbering:
                ax.text(*point, f"node {i_point}", ha="center", va="center")

        # Plot line segments
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
        """Overload of 'to_dict' for numpy usage and memory perf."""

        dict_ = self.base_dict()
        dict_["vertices"] = self.vertices.flatten().tolist()
        dict_["triangles"] = self.triangles.flatten().tolist()

        return dict_

    @classmethod
    def dict_to_object(cls, dict_: JsonSerializable, *args, **kwargs) -> "MeshType":
        """Overload of 'dict_to_object' for numpy usage and memory perf."""

        vertices = np.array(dict_["vertices"]).reshape(-1, 3)
        triangles = np.array(dict_["triangles"]).reshape(-1, 3)
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

        DessiaObject.__init__(self, name=name)

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

        self._faces = None

        PhysicalObject.__init__(self, name=name)

    def volmdlr_primitives(self, **kwargs):
        return [self]

    def area(self):
        """
        Return the area as the sum of areas of triangles.
        """
        areas = np.sqrt((self.triangles_cross_products() ** 2).sum(axis=1)) / 2.0
        return areas.sum()

    @property
    def faces(self):
        """
        Get the mesh faces as Triangle3D objects.

        :return: The triangles comosing the mesh.
        :rtype: list[Triangle3D]
        """
        if not self._faces:
            self._faces = self.to_triangles3d()
        return self._faces

    # IMPORT
    @classmethod
    def from_trimesh(cls, trimesh_: Trimesh) -> "Mesh3D":
        return cls(trimesh_.vertices, trimesh_.faces)

    @classmethod
    def from_trimesh_scene(cls, trimesh_scene: trimesh.Scene) -> "Mesh3D":
        mesh = cls(np.array([]), np.array([]))
        for trimesh_ in trimesh_scene.geometry.values():
            mesh += cls.from_trimesh(trimesh_)

        return mesh

    @classmethod
    def from_stl_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        return cls.from_trimesh(trimesh.load(filepath, "stl")).resize(scale_factor)

    @classmethod
    def from_stl_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "stl")).resize(scale_factor)

    @classmethod
    def from_obj_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        return cls.from_trimesh(trimesh.load(filepath, "obj")).resize(scale_factor)

    @classmethod
    def from_obj_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "obj")).resize(scale_factor)

    @classmethod
    def from_ply_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        return cls.from_trimesh(trimesh.load(filepath, "ply")).resize(scale_factor)

    @classmethod
    def from_ply_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "ply")).resize(scale_factor)

    @classmethod
    def from_off_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        return cls.from_trimesh(trimesh.load(filepath, "off")).resize(scale_factor)

    @classmethod
    def from_off_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "off")).resize(scale_factor)

    @classmethod
    def from_3mf_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        return cls.from_trimesh_scene(trimesh.load(filepath, "3mf")).resize(scale_factor)

    @classmethod
    def from_3mf_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        stream.seek(0)
        return cls.from_trimesh_scene(trimesh.load(stream, "3mf")).resize(scale_factor)

    # EXPORT
    def triangular_faces(self):
        """
        Export the mesh faces as Triangle3D objects.

        :return: The triangles comosing the mesh.
        :rtype: list[Triangle3D]
        """
        warnings.warn("Deprecated: use to_triangles3d instead.", DeprecationWarning)
        return self.to_triangles3d()

    def to_triangles3d(self):
        """
        Export the mesh faces as Triangle3D objects.

        :return: The triangles comosing the mesh.
        :rtype: list[Triangle3D]
        """
        triangles3d = []
        for vertex1, vertex2, vertex3 in self.remove_degenerate_triangles(tol=1e-6).triangles_vertices():
            point1 = volmdlr.Point3D(*vertex1)
            point2 = volmdlr.Point3D(*vertex2)
            point3 = volmdlr.Point3D(*vertex3)

            triangles3d.append(volmdlr.faces.Triangle3D(point1, point2, point3))

        return triangles3d

    def to_closed_shell(self):
        """
        Convert the Mesh3D object to a closed triangle shell.

        :return: A closed triangle shell representation of the Mesh3D object.
        :rtype: ClosedTriangleShell3D
        """
        warnings.warn(
            """
            ClosedTriangleShell3D is not an efficient object to deal with mesh data.
            Try to stick to Mesh3D or Trimesh object if you can.
            """
        )
        return volmdlr.shells.ClosedTriangleShell3D(faces=self.to_triangles3d(), name=self.name)

    def to_open_shell(self):
        """
        Convert the Mesh3D object to an open triangle shell.

        :return: An open triangle shell representation of the Mesh3D object.
        :rtype: OpenTriangleShell3D
        """
        warnings.warn(
            """
            OpenTriangleShell3D is not an efficient object to deal with mesh data.
            Try to stick to Mesh3D or Trimesh object if you can.
            """
        )
        return volmdlr.shells.OpenTriangleShell3D(faces=self.to_triangles3d(), name=self.name)

    def to_trimesh(self):
        return Trimesh(self.vertices, self.triangles)

    def to_babylon(self):
        """
        Convert the mesh in babylonjs format.

        https://doc.babylonjs.com/how_to/custom
        """
        mesh = self.round_vertices(decimals=6)
        babylon_mesh = {"positions": mesh.vertices.flatten().tolist(), "indices": mesh.triangles.flatten().tolist()}

        return babylon_mesh

    def babylon_meshes(self, merge_meshes=True):
        babylon_param = {"alpha": 1.0, "name": self.name, "color": [0.8, 0.8, 0.8]}
        babylon_mesh = self.to_babylon()
        babylon_mesh.update(babylon_param)

        return [babylon_mesh]

    # SAVING
    def save_to_stl_file(self, filepath: str, scale_factor: float = 1000.0):
        if not filepath.lower().endswith(".stl"):
            filepath += ".stl"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_stl_stream(file, scale_factor=scale_factor)

    def save_to_stl_stream(self, stream, scale_factor: float = 1000.0):
        self.resize(scale_factor).to_trimesh().export(stream, "stl")

    def save_to_obj_file(self, filepath: str, scale_factor: float = 1000.0):
        if not filepath.lower().endswith(".obj"):
            filepath += ".obj"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_obj_stream(file, scale_factor=scale_factor)

    def save_to_obj_stream(self, stream, scale_factor: float = 1000.0):
        self.resize(scale_factor).to_trimesh().export(stream, "obj")

    def save_to_ply_file(self, filepath: str, scale_factor: float = 1000.0):
        if not filepath.lower().endswith(".ply"):
            filepath += ".ply"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_ply_stream(file, scale_factor=scale_factor)

    def save_to_ply_stream(self, stream, scale_factor: float = 1000.0):
        self.resize(scale_factor).to_trimesh().export(stream, "ply")

    def save_to_off_file(self, filepath: str, scale_factor: float = 1000.0):
        if not filepath.lower().endswith(".off"):
            filepath += ".off"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_off_stream(file, scale_factor=scale_factor)

    def save_to_off_stream(self, stream, scale_factor: float = 1000.0):
        self.resize(scale_factor).to_trimesh().export(stream, "off")

    def save_to_3mf_file(self, filepath: str, scale_factor: float = 1000.0):
        if not filepath.lower().endswith(".3mf"):
            filepath += ".3mf"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_3mf_stream(file, scale_factor=scale_factor)

    def save_to_3mf_stream(self, stream, scale_factor: float = 1000.0):
        self.resize(scale_factor).to_trimesh().export(stream, "3mf")
