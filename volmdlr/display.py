#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes to define mesh for display use. Display mesh do not require good aspect ratios on elements.
"""

import warnings
from typing import List, Tuple, TypeVar, Union

import numpy as np
import pyfqmr
import trimesh
from dessia_common.core import DessiaObject
from dessia_common.serialization import BinaryFile
from dessia_common.typings import JsonSerializable
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from trimesh import Trimesh

import volmdlr.edges
from volmdlr.core import Primitive3D

# TODO: make this module "mesh" as it is not useful only for display


class MeshMixin:
    """
    Mixin class for 2D and 3D meshes.

    This is an abstract class.
    """

    MeshType = TypeVar("MeshType", bound="MeshMixin")

    _standalone_in_db = True
    _non_serializable_attributes = ["vertices", "triangles"]

    @property
    def n_vertices(self) -> int:
        """Get number of vertices in the mesh."""
        return len(self.vertices)

    @property
    def n_triangles(self) -> int:
        """Get number of triangles in the mesh."""
        return len(self.triangles)

    @property
    def dimension(self) -> int:
        """Get the dimension of the mesh ("2" for 2D mesh or "3" for 3D mesh)."""
        return self.vertices.shape[1]

    # MANIPULATION
    def resize(self, scale_factor: float) -> MeshType:
        """
        Resize the Mesh instance by scaling its vertices.

        :param scale_factor: The factor by which to scale the mesh.
        :type scale_factor: float

        :return: A new Mesh instance representing the scaled mesh.
        :rtype: MeshType
        """
        return self.__class__(self.vertices * scale_factor, self.triangles, name=self.name)

    def round_vertices(self, decimals: int = 9) -> MeshType:
        """
        Round the vertices of the Mesh instance to a given number of decimals.

        :param decimals: The number of decimal places to round the vertices to (default is 9).
        :type decimals: int, optional

        :return: A new Mesh instance with rounded vertices.
        :rtype: MeshType
        """
        rounded_vertices = np.round(self.vertices, decimals)

        return self.__class__(rounded_vertices, self.triangles, name=self.name)

    def remove_degenerate_triangles(self, tol: float = 0.0) -> MeshType:
        """
        Remove degenerate triangles from the Mesh instance.

        Degenerate triangles are triangles with vertices that are too close to each other.
        This method checks for degenerate triangles based on a tolerance value and removes them.

        :param tol: The tolerance value to determine whether a triangle is degenerate or not (default is 0.0).
        :type tol: float, optional

        :return: A new Mesh instance with degenerate triangles removed.
        :rtype: MeshType
        """
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
        return self.__class__(self.vertices, valid_triangles, name=self.name)

    def merge(self, other: MeshType, merge_vertices: bool = False, merge_triangles: bool = False) -> MeshType:
        """
        Merge two meshes.

        :param other: Another Mesh instance to merge with this instance.
        :type other: MeshType
        :param merge_vertices: Flag to indicate whether to merge vertices.
        :type merge_vertices: bool, optional
        :param merge_triangles: Flag to indicate whether to merge triangles.
        :type merge_triangles: bool, optional

        :return: A new Mesh instance representing the merged meshes.
        :rtype: MeshType
        """
        if self.__class__.__name__ != other.__class__.__name__:
            raise ValueError("Meshes should have same dimension.")

        if len(self.vertices) == 0 or len(self.triangles) == 0:
            return other
        if len(other.vertices) == 0 or len(other.triangles) == 0:
            return self

        merged_vertices = np.concatenate((self.vertices, other.vertices))
        merged_triangles = np.concatenate((self.triangles, other.triangles + len(self.vertices)))

        mesh = self.__class__(merged_vertices, merged_triangles, name=self.name)

        if merge_vertices:
            mesh = mesh.merge_vertices()
        if merge_triangles:
            mesh = mesh.merge_triangles()

        return mesh

    def merge_vertices(self) -> MeshType:
        """
        Merge duplicated vertices in the Mesh instance and remap triangles accordingly.

        This method identifies duplicate vertices and combines them into a single unique set of vertices,
        updating the triangles to use the unique vertices.

        :return: A new Mesh instance with merged vertices and updated triangles.
        :rtype: MeshType
        """
        unique_vertices, indices_map = np.unique(self.vertices, axis=0, return_inverse=True)
        remapped_triangles = indices_map[self.triangles]

        return self.__class__(unique_vertices, remapped_triangles, name=self.name)

    def merge_triangles(self) -> MeshType:
        """
        Merge duplicated triangles in the Mesh instance.

        This method identifies and removes duplicate triangles, resulting in a Mesh with unique triangles.

        :return: A new Mesh instance with merged triangles.
        :rtype: MeshType
        """

        sorted_triangles = np.sort(self.triangles, axis=1)
        _, unique_triangle_indices = np.unique(sorted_triangles, axis=0, return_index=True)
        unique_triangles = self.triangles[unique_triangle_indices]

        return self.__class__(self.vertices, unique_triangles, name=self.name)

    def split_shared_vertices(self) -> MeshType:
        """
        Split the shared vertices between triangles in the Mesh instance.

        This method recreates distinct vertices for each triangle, effectively unmerging shared vertices.
        The resulting mesh will have three times the number of vertices as the number of triangles.

        :return: A new Mesh instance with unmerged vertices and original triangles.
        :rtype: MeshType
        """
        unmerged_vertices = self.vertices[self.triangles.ravel()]
        unmerged_triangles = np.arange(len(self.triangles) * 3).reshape(-1, 3)

        return self.__class__(unmerged_vertices, unmerged_triangles, name=self.name)

    def __add__(self, other: MeshType) -> MeshType:
        """
        Overload the "+" operator to merge two Mesh instances, without mutualization of vertices and triangles.

        :param other: Another Mesh instance to concatenate with this instance.
        :type other: MeshType

        :return: A new Mesh instance representing the merged shells.
        :rtype: MeshType
        """
        return self.merge(other, merge_vertices=False, merge_triangles=False)

    def __or__(self, other: MeshType) -> MeshType:
        """
        Overload the "|" operator to merge two Mesh instances, with mutualization of vertices and triangles.

        :param other: Another Mesh instance to concatenate with this instance.
        :type other: MeshType

        :return: A new Mesh instance representing the concatenated shells.
        :rtype: MeshType
        """
        return self.merge(other, merge_vertices=True, merge_triangles=True)

    @classmethod
    def from_meshes(
        cls, meshes: List[MeshType], merge_vertices: bool = False, merge_triangles: bool = False
    ) -> MeshType:
        """
        Merge two meshes.

        :param meshes: A list of Mesh instance to merge all together.
        :type meshes: MeshType
        :param merge_vertices: Flag to indicate whether to merge vertices.
        :type merge_vertices: bool, optional
        :param merge_triangles: Flag to indicate whether to merge triangles.
        :type merge_triangles: bool, optional

        :return: A new Mesh instance representing the merged meshes.
        :rtype: MeshType
        """
        merged_mesh = cls(np.array([]), np.array([]))
        for mesh in meshes:
            merged_mesh = merged_mesh.merge(mesh, merge_vertices=merge_vertices, merge_triangles=merge_triangles)

        return merged_mesh

    # CHECK
    def check_consistency(self) -> bool:
        """
        Check the consistency of the Mesh instance.

        This method verifies that all vertices referenced by triangles are within the valid range of vertex indices.

        :return: True if the mesh is consistent, False otherwise.
        :rtype: bool
        """
        max_vertex_indices = np.max(self.triangles, axis=None)
        return np.all(max_vertex_indices < len(self.vertices))

    # COMPUTATION
    def triangles_vertices(self):
        """
        Get the actual triangles of the mesh represented by their vertices (not indices).

        :return: An array containing the vertices of the triangles.
        :rtype: np.ndarray[float]
        """
        triangles = self.vertices.view(np.ndarray)[self.triangles]
        return triangles

    def triangles_cross_products(self):
        """
        Compute the cross products of edges for each triangle in the mesh.

        :return: An array containing the cross products of edges for each triangle.
        :rtype: np.ndarray[float]
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
    def dict_to_object(cls, dict_: JsonSerializable, *args, **kwargs) -> MeshType:
        """Overload of 'dict_to_object' for numpy usage and memory perf."""

        vertices = np.array(dict_["vertices"]).reshape(-1, 3)
        triangles = np.array(dict_["triangles"]).reshape(-1, 3)
        name = dict_["name"]

        return cls(vertices=vertices, triangles=triangles, name=name)

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
    2D triangle mesh.
    """

    _linesegment_class = volmdlr.edges.LineSegment2D
    _point_class = volmdlr.Point2D

    def __init__(self, vertices: NDArray[float], triangles: NDArray[int], name: str = ""):
        """
        Initialize a 2D mesh.

        :param vertices: An array of 2D vertices specifying the 2D mesh.
        :type vertices: ndarray[float]
        :param triangles: An array of triangles representing the connectivity of the 2D mesh.
        :type triangles: ndarray[int]
        :param name: A name for the mesh (default is an empty string).
        :type name: str, optional
        """
        self.vertices = vertices
        self.triangles = triangles

        DessiaObject.__init__(self, name=name)

    def area(self) -> float:
        """
        Calculate the total area of the 2D mesh as the sum of areas of triangles.

        :return: The total area of the mesh.
        :rtype: float
        """
        areas = np.sqrt((self.triangles_cross_products() ** 2)) / 2.0
        return areas.sum()


class Mesh3D(MeshMixin, Primitive3D):
    """
    3D triangle mesh.
    """

    # pylint: disable=too-many-public-methods

    _linesegment_class = volmdlr.edges.LineSegment3D
    _point_class = volmdlr.Point3D

    def __init__(
        self,
        vertices: NDArray[float],
        triangles: NDArray[int],
        color: Tuple[float, float, float] = None,
        alpha: float = 1.0,
        name: str = "",
    ):
        """
        Initialize a 3D mesh.

        :param vertices: An array of 3D vertices specifying the 3D mesh.
        :param triangles: An array of triangles representing the connectivity of the 3D mesh.
        :param color: A color for the mesh, optional.
        :param alpha: An alpha value for the mesh, optional.
        :param name: A name for the mesh, optional (default is an empty string).
        """
        self.vertices = vertices
        self.triangles = triangles

        self._faces = None
        self._bounding_box = None

        Primitive3D.__init__(self, color=color, alpha=alpha, name=name)

    def triangulation(self) -> "Mesh3D":
        """Return self as triangulation to enable VolumeModel usage."""
        return self

    def to_babylon(self):
        """
        Convert the mesh to the Babylon.js format.

        This method rounds the vertices to 6 decimal places and returns the mesh in a Babylon.js compatible format.
        https://doc.babylonjs.com/how_to/custom

        :return: A dictionary representing the mesh in Babylon.js format with 'positions' and 'indices' keys.
        :rtype: dict
        """
        mesh = self.round_vertices(decimals=6)
        babylon_mesh = {"positions": mesh.vertices.flatten().tolist(), "indices": mesh.triangles.flatten().tolist()}

        return babylon_mesh

    @property
    def bounding_box(self):
        """Bounding box of current mesh."""
        if self._bounding_box is None:
            maximums = np.max(self.vertices, axis=0)
            minimums = np.min(self.vertices, axis=0)

            self._bounding_box = volmdlr.core.BoundingBox(
                minimums[0], maximums[0], minimums[1], maximums[1], minimums[2], maximums[2]
            )

        return self._bounding_box

    def area(self) -> float:
        """
        Calculate the total surface area of the 3D mesh as the sum of areas of triangles.

        :return: The total surface area of the 3D mesh.
        :rtype: float
        """
        areas = np.sqrt((self.triangles_cross_products() ** 2).sum(axis=1)) / 2.0
        return areas.sum()

    def minimum_distance(self, other_mesh: "Mesh3D", return_points: bool = False):
        """
        Compute the minimum distance between this 3D mesh and another 3D mesh.

        This is an approximation: only vertices are taken in account for minimum distance computation.

        :param other_mesh: The other 3D mesh to compare against.
        :type other_mesh: Mesh3D
        :param return_points: Whether to return the closest points.
        :type return_points: bool, optional

        :return: The minimum distance between the two meshes, and optionally, the closest points.
        :rtype: float or (float, ndarray[float], ndarray[float])
        """
        # Create KD-Trees for both meshes (cKDTree for improved performance)
        self_tree = cKDTree(self.vertices)
        other_tree = cKDTree(other_mesh.vertices)

        # Query the KD-Tree to find the nearest neighbors for all vertices in one go
        _, self_to_other_indices = other_tree.query(self.vertices, k=1)
        _, other_to_self_indices = self_tree.query(other_mesh.vertices, k=1)

        # Calculate the minimum distance between vertices using vectorized operations
        self_to_other_distances = np.linalg.norm(self.vertices - other_mesh.vertices[self_to_other_indices], axis=1)
        other_to_self_distances = np.linalg.norm(other_mesh.vertices - self.vertices[other_to_self_indices], axis=1)
        min_distances = [self_to_other_distances.min(), other_to_self_distances.min()]
        if min_distances[0] < min_distances[1]:
            min_distance = min_distances[0]
            # Get the points corresponding to the minimum distance
            min_distance_index = np.argmin(self_to_other_distances)
            closest_point_self = self.vertices[min_distance_index]
            closest_point_other = other_mesh.vertices[self_to_other_indices[min_distance_index]]
        else:
            min_distance = min_distances[1]
            # Get the points corresponding to the minimum distance
            min_distance_index = np.argmin(other_to_self_distances)
            closest_point_self = self.vertices[other_to_self_indices[min_distance_index]]
            closest_point_other = other_mesh.vertices[min_distance_index]

        if return_points:

            closest_point_self = volmdlr.Point3D(closest_point_self[0],
                                                 closest_point_self[1], closest_point_self[2])
            closest_point_other = volmdlr.Point3D(closest_point_other[0],
                                                  closest_point_other[1], closest_point_other[2])
            return min_distance, closest_point_self, closest_point_other

        return min_distance

    def get_edges_triangles(self):
        """
        Compute lengths edges of triangles.

        :return: A 3D numpy array representing edges of triangles. The dimensions are n_triangles x 3 x 2,
             where each entry contains the start and end points of an edge.
        :rtype: np.ndarray
        """
        edges = np.stack([self.triangles[:, [0, 1]], self.triangles[:, [0, 2]], self.triangles[:, [1, 2]]], axis=1)

        return edges

    def compute_len_edges(self):
        """
        Compute the lengths of edges for each triangle in the mesh.

        :return: Lengths of edges (3 edges per triangles) of dimensions n_simplices x 3 and edges of dimensions
            n_simplices x 3 x 2
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        edges = self.get_edges_triangles()

        indexed_points = self.vertices[edges]
        vectors = indexed_points[..., 0, :] - indexed_points[..., 1, :]
        return np.linalg.norm(vectors, axis=-1), edges

    def get_mesh_border(self):
        """
        Retrieve the topological border of a triangle mesh.

        This function identifies and returns the edges that belong to only one triangle,
        effectively representing the border of the mesh.

        :return: A tuple of two numpy arrays. The first array contains the unique border edges,
            and the second array includes all edges of the mesh.
        :rtype: tuple[np.ndarray, np.ndarray]
        """
        edges = self.get_edges_triangles().reshape((-1, 2))
        unique_edges, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
        border_edges = unique_edges[counts == 1]

        return border_edges, edges

    def remove_large_triangles(self, threshold_edge_length: float) -> "Mesh3D":
        """
        Remove triangles from the mesh whose edge lengths exceed the specified threshold.

        :param threshold_edge_length: The maximum allowed edge length for a triangle to remain in the mesh.
        :type threshold_edge_length: float

        :return: A new Mesh3D instance with large triangles removed.
        :rtype: Mesh3D
        """

        # Compute the lengths of all edges in the mesh
        edge_lengths, _ = self.compute_len_edges()

        # Find triangles where all edges are below the threshold
        valid_triangles = np.all(edge_lengths < threshold_edge_length, axis=1)

        # Keep only the triangles that are valid
        return Mesh3D(self.vertices, self.triangles[valid_triangles])

    def decimate(
        self,
        target_count: int,
        update_rate: int = 5,
        aggressiveness: float = 7.0,
        max_iterations: int = 100,
        verbose: bool = False,
        lossless: bool = False,
        threshold_lossless: float = 1e-3,
        alpha: float = 1e-9,
        k: int = 3,
        preserve_border: bool = True,
    ) -> "Mesh3D":
        """
        Decimate the Mesh3D, and return it as a new instance.

        Vertices of the mesh should be merged (and maybe rounded) for efficient decimation.

        Note: threshold = alpha * pow(iteration + k, aggressiveness)

        :param target_count: Target number of triangles. Not used if `lossless` is True.
        :type target_count: int
        :param update_rate: Number of iterations between each update. If `lossless` flag is set to True, rate is 1.
        :type update_rate: int
        :param aggressiveness: Parameter controlling the growth rate of the threshold at each iteration when `lossless`
            is False.
        :type aggressiveness: float
        :param max_iterations: Maximal number of iterations.
        :type max_iterations: int
        :param verbose: Control verbosity.
        :type verbose: bool
        :param lossless: Use the lossless simplification method.
        :type lossless: bool
        :param threshold_lossless: Maximal error after which a vertex is not deleted. Only for `lossless` method.
        :type threshold_lossless: float
        :param alpha: Parameter for controlling the threshold growth.
        :type alpha: float
        :param k: Parameter for controlling the threshold growth.
        :type k: int
        :param preserve_border: Flag for preserving vertices on open border.
        :type preserve_border: bool

        :return: The decimated mesh.
        :rtype: Mesh3D
        """
        # pylint: disable=too-many-arguments

        simplifier = pyfqmr.Simplify()
        simplifier.setMesh(self.vertices, self.triangles)
        simplifier.simplify_mesh(
            target_count=target_count,
            update_rate=update_rate,
            aggressiveness=aggressiveness,
            max_iterations=max_iterations,
            verbose=verbose,
            lossless=lossless,
            threshold_lossless=threshold_lossless,
            alpha=alpha,
            K=k,
            preserve_border=preserve_border,
        )

        vertices, triangles, _ = simplifier.getMesh()

        return self.__class__(vertices, triangles)

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
        """
        Create a 3D mesh from a Trimesh object.

        :param trimesh_: A Trimesh object representing the 3D mesh.
        :type trimesh_: Trimesh

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        return cls(trimesh_.vertices, trimesh_.faces)

    @classmethod
    def from_trimesh_scene(cls, trimesh_scene: trimesh.Scene, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from a Trimesh Scene.

        :param trimesh_scene: A Trimesh Scene containing multiple geometry objects.
        :type trimesh_scene: trimesh.Scene
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        return cls.from_meshes(cls.trimesh_scene_to_meshes(trimesh_scene, scale_factor))

    @classmethod
    def trimesh_scene_to_meshes(cls, trimesh_scene: trimesh.Scene, scale_factor: float = 0.001) -> List["Mesh3D"]:
        """
        Create a 3D mesh from a Trimesh Scene.

        :param trimesh_scene: A Trimesh Scene containing multiple geometry objects.
        :type trimesh_scene: trimesh.Scene
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A list of new 3D mesh instance.
        :rtype: list[Mesh3D]
        """
        meshes = []
        for trimesh_ in trimesh_scene.geometry.values():
            meshes.append(cls.from_trimesh(trimesh_).resize(scale_factor))

        return meshes

    @classmethod
    def from_stl_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an STL file.

        :param filepath: The path to the STL file.
        :type filepath: str
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        return cls.from_trimesh(trimesh.load(filepath, "stl")).resize(scale_factor)

    @classmethod
    def from_stl_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an STL stream.

        :param stream: A binary stream containing STL data.
        :type stream: BinaryFile
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "stl")).resize(scale_factor)

    @classmethod
    def from_obj_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an OBJ file.

        :param filepath: The path to the OBJ file.
        :type filepath: str
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        return cls.from_trimesh(trimesh.load(filepath, "obj")).resize(scale_factor)

    @classmethod
    def from_obj_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an OBJ stream.

        :param stream: A binary stream containing OBJ data.
        :type stream: BinaryFile
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "obj")).resize(scale_factor)

    @classmethod
    def from_ply_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an PLY file.

        :param filepath: The path to the PLY file.
        :type filepath: str
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        return cls.from_trimesh(trimesh.load(filepath, "ply")).resize(scale_factor)

    @classmethod
    def from_ply_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an PLY stream.

        :param stream: A binary stream containing PLY data.
        :type stream: BinaryFile
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "ply")).resize(scale_factor)

    @classmethod
    def from_off_file(cls, filepath: str, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an OFF file.

        :param filepath: The path to the OFF file.
        :type filepath: str
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        return cls.from_trimesh(trimesh.load(filepath, "off")).resize(scale_factor)

    @classmethod
    def from_off_stream(cls, stream: BinaryFile, scale_factor: float = 0.001) -> "Mesh3D":
        """
        Create a 3D mesh from an OFF stream.

        :param stream: A binary stream containing OFF data.
        :type stream: BinaryFile
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        stream.seek(0)
        return cls.from_trimesh(trimesh.load(stream, "off")).resize(scale_factor)

    @classmethod
    def from_3mf_file(
        cls, filepath: str, scale_factor: float = 0.001, merge_meshes: bool = True
    ) -> Union["Mesh3D", List["Mesh3D"]]:
        """
        Create a 3D mesh from an 3MF file.

        :param filepath: The path to the 3MF file.
        :type filepath: str
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional
        :param merge_meshes: A flag to choose to merge all the 3mf meshes in one, or return a list of meshes.
        :type merge_meshes: bool

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        if merge_meshes:
            return cls.from_trimesh_scene(trimesh.load(filepath, "3mf"), scale_factor)
        return cls.trimesh_scene_to_meshes(trimesh.load(filepath, "3mf"), scale_factor)

    @classmethod
    def from_3mf_stream(
        cls, stream: BinaryFile, scale_factor: float = 0.001, merge_meshes: bool = True
    ) -> Union["Mesh3D", List["Mesh3D"]]:
        """
        Create a 3D mesh from an 3MF stream.

        :param stream: A binary stream containing 3MF data.
        :type stream: BinaryFile
        :param scale_factor: The scale factor to apply to the mesh (default is 0.001).
        :type scale_factor: float, optional
        :param merge_meshes: A flag to choose to merge all the 3mf meshes in one, or return a list of meshes.
        :type merge_meshes: bool

        :return: A new 3D mesh instance.
        :rtype: Mesh3D
        """
        stream.seek(0)

        if merge_meshes:
            return cls.from_trimesh_scene(trimesh.load(stream, "3mf"), scale_factor)
        return cls.trimesh_scene_to_meshes(trimesh.load(stream, "3mf"), scale_factor)

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
        # pylint: disable=import-outside-toplevel, cyclic-import
        from volmdlr.faces import Triangle3D

        triangles3d = []
        for vertex1, vertex2, vertex3 in self.remove_degenerate_triangles(tol=1e-6).triangles_vertices():
            point1 = volmdlr.Point3D(*vertex1)
            point2 = volmdlr.Point3D(*vertex2)
            point3 = volmdlr.Point3D(*vertex3)

            triangles3d.append(Triangle3D(point1, point2, point3))

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
            Try to stick to Mesh3D or Trimesh object if possible.
            """
        )

        # pylint: disable=import-outside-toplevel, cyclic-import
        from volmdlr.shells import ClosedTriangleShell3D

        return ClosedTriangleShell3D(faces=self.to_triangles3d(), name=self.name)

    def to_open_shell(self):
        """
        Convert the Mesh3D object to an open triangle shell.

        :return: An open triangle shell representation of the Mesh3D object.
        :rtype: OpenTriangleShell3D
        """
        warnings.warn(
            """
            OpenTriangleShell3D is not an efficient object to deal with mesh data.
            Try to stick to Mesh3D or Trimesh object if possible.
            """
        )

        # pylint: disable=import-outside-toplevel, cyclic-import
        from volmdlr.shells import OpenTriangleShell3D

        return OpenTriangleShell3D(faces=self.to_triangles3d(), name=self.name)

    def to_trimesh(self):
        """
        Convert the Mesh3D instance to a Trimesh object.

        :return: A Trimesh object representing the 3D mesh.
        :rtype: Trimesh
        """
        return Trimesh(self.vertices, self.triangles)

    # SAVING
    def save_to_stl_file(self, filepath: str, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to an STL file.

        :param filepath: The path to the STL file.
        :type filepath: str

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        if not filepath.lower().endswith(".stl"):
            filepath += ".stl"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_stl_stream(file, scale_factor=scale_factor)

    def save_to_stl_stream(self, stream, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to an STL stream.

        :param stream: A binary stream to write the STL data.
        :type stream: BinaryFile

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        self.resize(scale_factor).to_trimesh().export(stream, "stl")

    def save_to_obj_file(self, filepath: str, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to an OBJ file.

        :param filepath: The path to the OBJ file.
        :type filepath: str

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        if not filepath.lower().endswith(".obj"):
            filepath += ".obj"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_obj_stream(file, scale_factor=scale_factor)

    def save_to_obj_stream(self, stream, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to an OBJ stream.

        :param stream: A binary stream to write the OBJ data.
        :type stream: BinaryFile

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        self.resize(scale_factor).to_trimesh().export(stream, "obj")

    def save_to_ply_file(self, filepath: str, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to a PLY file.

        :param filepath: The path to the PLY file.
        :type filepath: str

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        if not filepath.lower().endswith(".ply"):
            filepath += ".ply"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_ply_stream(file, scale_factor=scale_factor)

    def save_to_ply_stream(self, stream, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to a PLY stream.

        :param stream: A binary stream to write the PLY data.
        :type stream: BinaryFile

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        self.resize(scale_factor).to_trimesh().export(stream, "ply")

    def save_to_off_file(self, filepath: str, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to an OFF file.

        :param filepath: The path to the OFF file.
        :type filepath: str

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        if not filepath.lower().endswith(".off"):
            filepath += ".off"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_off_stream(file, scale_factor=scale_factor)

    def save_to_off_stream(self, stream, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to an OFF stream.

        :param stream: A binary stream to write the OFF data.
        :type stream: BinaryFile

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        self.resize(scale_factor).to_trimesh().export(stream, "off")

    def save_to_3mf_file(self, filepath: str, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to a 3MF file.

        :param filepath: The path to the 3MF file.
        :type filepath: str

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        if not filepath.lower().endswith(".3mf"):
            filepath += ".3mf"
            print(f"Changing name to {filepath}")

        with open(filepath, "wb") as file:
            self.save_to_3mf_stream(file, scale_factor=scale_factor)

    def save_to_3mf_stream(self, stream, scale_factor: float = 1000.0):
        """
        Save the 3D mesh to a 3MF stream.

        :param stream: A binary stream to write the 3MF data.
        :type stream: BinaryFile

        :param scale_factor: The scale factor to apply to the mesh (default is 1000.0).
        :type scale_factor: float, optional
        """
        self.resize(scale_factor).to_trimesh().export(stream, "3mf")
