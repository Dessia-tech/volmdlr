"""
Class for voxel representation of volmdlr models
"""
import math
import warnings
from typing import Dict, Iterable, List, Set, Tuple

import numpy as np
from dessia_common.core import PhysicalObject
from tqdm import tqdm

from volmdlr import Point2D, Point3D, Vector3D
from volmdlr.core import VolumeModel
from volmdlr.faces import PlaneFace3D, Triangle3D
from volmdlr.shells import ClosedShell3D, ClosedTriangleShell3D
from volmdlr.surfaces import PLANE3D_OXY, PLANE3D_OXZ, PLANE3D_OYZ, Surface2D
from volmdlr.wires import ClosedPolygon2D
from volmdlr.voxelization_compiled import triangle_intersects_voxel, aabb_intersecting_boxes

# Custom types
Point = Tuple[float, ...]
Triangle = Tuple[Point, ...]
Segment = Tuple[Point, ...]


class Voxelization(PhysicalObject):
    """
    Class for creation and manipulation of voxelization of volmdlr geometry.

    The voxelization is defined in an implicit grid of the 3D space. The grid is defined by the size of the voxels.
    The implicit grid is defined by discretizing the 3D space into axis aligned cubes with given size, with the
    origin of the global landmark (0, 0, 0) always being a corner point (and not inside a cube).

    For example, in 1D and with a voxel size of t, the set of voxels (defined by minimal and maximal point) can only
    be defined as {i ∈ N, t ∈ R | (i * t, (i+1) * t)}
    The corresponding set of voxel centers is defined by the following set: {i ∈ N, t ∈ R | (i + 0.5) * t}

    This approach allows to always define voxelization of same size in the same grid of the 3D space, which is very
    useful to perform very fast Boolean operations.
    """

    def __init__(self, voxels_centers: Set[Point], voxel_size: float, octree_root: "OctreeNode" = None, name: str = ""):
        """
        Initialize the Voxelization.

        :param voxels_centers: The set of points representing voxel centers.
        :type voxels_centers: set[tuple[float, float, float]]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param octree_root: The octree root used to create the voxelization, if so.
        :type octree_root: OctreeNode, optional
        :param name: The name of the Voxelization.
        :type name: str, optional
        """
        self.voxels_centers = voxels_centers
        self.voxel_size = voxel_size
        self.octree_root = octree_root

        PhysicalObject.__init__(self, name=name)

    def __eq__(self, other_voxelization: "Voxelization") -> bool:
        """
        Check if the current voxelization is equal to another voxelization.

        :param other_voxelization: The voxelization to compare.
        :type other_voxelization: Voxelization

        :return: True if the voxelizations are equal, False otherwise.
        :rtype: bool
        """
        return (
            self.voxels_centers == other_voxelization.voxels_centers
            and self.voxel_size == other_voxelization.voxel_size
        )

    def __add__(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Return the union of the current voxelization with another voxelization.

        :param other_voxelization: The voxelization to union with.
        :type other_voxelization: Voxelization

        :return: The union of the voxelizations.
        :rtype: Voxelization
        """
        return self.union(other_voxelization)

    def __sub__(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Return the difference between the current voxelization and another voxelization.

        :param other_voxelization: The voxelization to subtract.
        :type other_voxelization: Voxelization

        :return: The difference between the voxelizations.
        :rtype: Voxelization
        """
        return self.difference(other_voxelization)

    def __and__(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Return the intersection of the current voxelization with another voxelization.

        :param other_voxelization: The voxelization to intersect with.
        :type other_voxelization: Voxelization

        :return: The intersection of the voxelizations.
        :rtype: Voxelization
        """
        return self.intersection(other_voxelization)

    def __or__(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Return the union of the current voxelization with another voxelization.

        :param other_voxelization: The voxelization to union with.
        :type other_voxelization: Voxelization

        :return: The union of the voxelizations.
        :rtype: Voxelization
        """
        return self.union(other_voxelization)

    def __xor__(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Return the symmetric difference between the current voxelization and another voxelization.

        :param other_voxelization: The voxelization to calculate the symmetric difference with.
        :type other_voxelization: Voxelization
        :return: The symmetric difference between the voxelizations.
        :rtype: Voxelization
        """
        return self.symmetric_difference(other_voxelization)

    def __invert__(self) -> "Voxelization":
        """
        Return the inverse of the current voxelization.

        :return: The inverse Voxelization object.
        :rtype: Voxelization
        """
        return self.inverse()

    def __len__(self):
        """
        Return the number of voxels in the voxelization.

        :return: The number of voxels.
        :rtype: int
        """
        return len(self.voxels_centers)

    @classmethod
    def from_closed_triangle_shell(
        cls, closed_triangle_shell: ClosedTriangleShell3D, voxel_size: float, method: str = "iterative", name: str = ""
    ) -> "Voxelization":
        """
        Create a Voxelization from a ClosedTriangleShell3D.

        :param closed_triangle_shell: The ClosedTriangleShell3D to voxelize.
        :type closed_triangle_shell: ClosedTriangleShell3D
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param method: The method used to voxelize the geometry ("iterative" or "octree"). Default is "octree".
        :type method: str, optional
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: Voxelization
        """
        octree_root_node = None
        if method == "iterative":
            triangles = cls._closed_triangle_shell_to_triangles(closed_triangle_shell)
            voxels = cls._triangles_to_voxels(triangles, voxel_size)

        elif method == "octree":
            triangles = cls._closed_triangle_shell_to_triangles(closed_triangle_shell)
            octree_root_node = OctreeNode.octree_voxelization_size_based(triangles, voxel_size)
            voxels = set(octree_root_node.get_leaf_centers())

        else:
            raise ValueError("Invalid 'method' argument: must be 'iterative' or 'octree'")

        return cls(voxels_centers=voxels, voxel_size=voxel_size, octree_root=octree_root_node, name=name)

    @classmethod
    def from_closed_shell(
        cls, closed_shell: ClosedShell3D, voxel_size: float, method: str = "iterative", name: str = ""
    ) -> "Voxelization":
        """
        Create a Voxelization from a ClosedShell3D.

        :param closed_shell: The ClosedShell3D to voxelize.
        :type closed_shell: ClosedShell3D
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param method: The method used to voxelize the geometry ("iterative" or "octree"). Default is "octree".
        :type method: str, optional
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: Voxelization
        """
        octree_root_node = None
        if method == "iterative":
            triangles = cls._closed_shell_to_triangles(closed_shell)
            voxels = cls._triangles_to_voxels(triangles, voxel_size)

        elif method == "octree":
            triangles = cls._closed_shell_to_triangles(closed_shell)
            octree_root_node = OctreeNode.octree_voxelization_size_based(triangles, voxel_size)
            voxels = set(octree_root_node.get_leaf_centers())

        else:
            raise ValueError("Invalid 'method' argument: must be 'iterative' or 'octree'")

        return cls(voxels_centers=voxels, voxel_size=voxel_size, octree_root=octree_root_node, name=name)

    @classmethod
    def from_volume_model(
        cls, volume_model: VolumeModel, voxel_size: float, method: str = "iterative", name: str = ""
    ) -> "Voxelization":
        """
        Create a Voxelization from a VolumeModel.

        :param volume_model: The VolumeModel to voxelize.
        :type volume_model: VolumeModel
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param method: The method used to voxelize the geometry ("iterative" or "octree"). Default is "octree".
        :type method: str, optional
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: Voxelization
        """
        octree_root_node = None
        if method == "iterative":
            triangles = cls._volume_model_to_triangles(volume_model)
            voxels = cls._triangles_to_voxels(triangles, voxel_size)

        elif method == "octree":
            triangles = cls._volume_model_to_triangles(volume_model)
            octree_root_node = OctreeNode.octree_voxelization_size_based(triangles, voxel_size)
            voxels = set(octree_root_node.get_leaf_centers())

        else:
            raise ValueError("Invalid 'method' argument: must be 'iterative' or 'octree'")

        return cls(voxels_centers=voxels, voxel_size=voxel_size, octree_root=octree_root_node, name=name)

    @staticmethod
    def _closed_triangle_shell_to_triangles(closed_triangle_shell: ClosedTriangleShell3D) -> List[Triangle]:
        """
        Helper method to convert a ClosedTriangleShell3D to a list of triangles.

        :param closed_triangle_shell: The ClosedTriangleShell3D to convert to triangles.
        :type closed_triangle_shell: ClosedTriangleShell3D

        :return: The list of triangles extracted from the ClosedTriangleShell3D.
        :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        return list(
            tuple(tuple([point.x, point.y, point.z]) for point in [tr.point1, tr.point2, tr.point3])
            for tr in closed_triangle_shell.primitives
        )

    @staticmethod
    def _closed_shell_to_triangles(closed_shell: ClosedShell3D) -> List[Triangle]:
        """
        Helper method to convert a ClosedShell3D to a list of triangles.
        It uses the "triangulation" method to triangulate the ClosedShell3D.

        :param closed_shell: The ClosedShell3D to convert to triangles.
        :type closed_shell: ClosedShell3D

        :return: The list of triangles extracted from the triangulated ClosedShell3D.
        :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangulation = closed_shell.triangulation()
        return list(
            tuple(
                tuple([point.x, point.y, point.z])
                for point in [
                    triangulation.points[triangle[0]],
                    triangulation.points[triangle[1]],
                    triangulation.points[triangle[2]],
                ]
            )
            for triangle in triangulation.triangles
        )

    @staticmethod
    def _volume_model_to_triangles(volume_model: VolumeModel) -> List[Triangle]:
        """
        Helper method to convert a VolumeModel to a list of triangles.
        It uses the "triangulation" method to triangulate the primitives of the VolumeModel.

        :param volume_model: The VolumeModel to convert to triangles.
        :type volume_model: VolumeModel

        :return: The list of triangles extracted from the triangulated primitives of the VolumeModel.
        :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangles = []
        for primitive in volume_model.get_shells():
            triangulation = primitive.triangulation()
            triangles += list(
                tuple(
                    tuple([point.x, point.y, point.z])
                    for point in [
                        triangulation.points[triangle[0]],
                        triangulation.points[triangle[1]],
                        triangulation.points[triangle[2]],
                    ]
                )
                for triangle in triangulation.triangles
            )
        return triangles

    @staticmethod
    def _triangle_intersects_voxel(triangle: Triangle, voxel_center: Point, voxel_extents: List[float]) -> bool:
        """
        Helper method to compute if there is an intersection between a 3D triangle and a voxel.
        This method uses the "Separating Axis Theorem".

        This method has been implemented in Cython in voxelization_compiled.pyx file.

        :param triangle: The triangle to check if it intersects with the voxel.
        :type: triangle: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
        :param voxel_center: The center point of the voxel.
        :type voxel_center: tuple[float, float, float]
        :param voxel_extents: The extents of the voxel in each direction (half-size of the voxel size).
        :type voxel_extents: list[float, float, float]

        :return: True if there is an intersection, False otherwise.
        :rtype: bool
        """
        # Method ported from https://gist.github.com/zvonicek/fe73ba9903f49d57314cf7e8e0f05dcf
        # pylint: disable=invalid-name,too-many-locals,too-many-return-statements,too-many-statements,too-many-branches

        triangle = np.array(triangle)
        box_center = np.array(voxel_center)
        box_extents = np.array(voxel_extents)

        x, y, z = 0, 1, 2

        # Translate triangle as conceptually moving AABB to origin
        v0 = triangle[0] - box_center
        v1 = triangle[1] - box_center
        v2 = triangle[2] - box_center

        # Compute edge vectors for triangle
        f0 = triangle[1] - triangle[0]
        f1 = triangle[2] - triangle[1]
        f2 = triangle[0] - triangle[2]

        # REGION TEST AXES a00..a22 (CATEGORY 3)

        # Test axis a00
        a00 = np.array([0, -f0[z], f0[y]])
        p0 = np.dot(v0, a00)
        p1 = np.dot(v1, a00)
        p2 = np.dot(v2, a00)
        r = box_extents[y] * abs(f0[z]) + box_extents[z] * abs(f0[y])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a01
        a01 = np.array([0, -f1[z], f1[y]])
        p0 = np.dot(v0, a01)
        p1 = np.dot(v1, a01)
        p2 = np.dot(v2, a01)
        r = box_extents[y] * abs(f1[z]) + box_extents[z] * abs(f1[y])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a02
        a02 = np.array([0, -f2[z], f2[y]])
        p0 = np.dot(v0, a02)
        p1 = np.dot(v1, a02)
        p2 = np.dot(v2, a02)
        r = box_extents[y] * abs(f2[z]) + box_extents[z] * abs(f2[y])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a10
        a10 = np.array([f0[z], 0, -f0[x]])
        p0 = np.dot(v0, a10)
        p1 = np.dot(v1, a10)
        p2 = np.dot(v2, a10)
        r = box_extents[x] * abs(f0[z]) + box_extents[z] * abs(f0[x])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a11
        a11 = np.array([f1[z], 0, -f1[x]])
        p0 = np.dot(v0, a11)
        p1 = np.dot(v1, a11)
        p2 = np.dot(v2, a11)
        r = box_extents[x] * abs(f1[z]) + box_extents[z] * abs(f1[x])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a12
        a11 = np.array([f2[z], 0, -f2[x]])
        p0 = np.dot(v0, a11)
        p1 = np.dot(v1, a11)
        p2 = np.dot(v2, a11)
        r = box_extents[x] * abs(f2[z]) + box_extents[z] * abs(f2[x])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a20
        a20 = np.array([-f0[y], f0[x], 0])
        p0 = np.dot(v0, a20)
        p1 = np.dot(v1, a20)
        p2 = np.dot(v2, a20)
        r = box_extents[x] * abs(f0[y]) + box_extents[y] * abs(f0[x])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a21
        a21 = np.array([-f1[y], f1[x], 0])
        p0 = np.dot(v0, a21)
        p1 = np.dot(v1, a21)
        p2 = np.dot(v2, a21)
        r = box_extents[x] * abs(f1[y]) + box_extents[y] * abs(f1[x])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # Test axis a22
        a22 = np.array([-f2[y], f2[x], 0])
        p0 = np.dot(v0, a22)
        p1 = np.dot(v1, a22)
        p2 = np.dot(v2, a22)
        r = box_extents[x] * abs(f2[y]) + box_extents[y] * abs(f2[x])
        if (max(-max(p0, p1, p2), min(p0, p1, p2))) > r:
            return False

        # ENDREGION

        # REGION TEST THE THREE AXES CORRESPONDING TO THE FACE NORMALS OF AABB B (CATEGORY 1)

        # Exit if...
        # ... [-extents.X, extents.X] and [min(v0.X,v1.X,v2.X), max(v0.X,v1.X,v2.X)] do not overlap
        if max(v0[x], v1[x], v2[x]) < -box_extents[x] or min(v0[x], v1[x], v2[x]) > box_extents[x]:
            return False

        # ... [-extents.Y, extents.Y] and [min(v0.Y,v1.Y,v2.Y), max(v0.Y,v1.Y,v2.Y)] do not overlap
        if max(v0[y], v1[y], v2[y]) < -box_extents[y] or min(v0[y], v1[y], v2[y]) > box_extents[y]:
            return False

        # ... [-extents.Z, extents.Z] and [min(v0.Z,v1.Z,v2.Z), max(v0.Z,v1.Z,v2.Z)] do not overlap
        if max(v0[z], v1[z], v2[z]) < -box_extents[z] or min(v0[z], v1[z], v2[z]) > box_extents[z]:
            return False

        # ENDREGION

        # REGION TEST SEPARATING AXIS CORRESPONDING TO TRIANGLE FACE NORMAL (CATEGORY 2)

        plane_normal = np.cross(f0, f1)
        plane_distance = abs(np.dot(plane_normal, v0))

        # Compute the projection interval radius of b onto L(t) = b.c + t * p.n
        r = (
            box_extents[x] * abs(plane_normal[x])
            + box_extents[y] * abs(plane_normal[y])
            + box_extents[z] * abs(plane_normal[z])
        )

        # Intersection occurs when plane distance falls within [-r,+r] interval
        if plane_distance > r:
            return False

        # ENDREGION

        return True

    @staticmethod
    def _aabb_intersecting_boxes(min_point: Point, max_point: Point, voxel_size: float) -> List[Point]:
        """
        Helper method to compute the center of the voxels that intersect with a given axis aligned
        bounding box (defined by 2 points).

        This method has been implemented in Cython in voxelization_compiled.pyx file.

        :param min_point: The minimum point of the bounding box.
        :type min_point: tuple[float, float, float]
        :param max_point: The maximum point of the bounding box.
        :type max_point: tuple[float, float, float]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: A list of the centers of the intersecting voxels.
        :rtype: list[tuple[float, float, float]]
        """
        # Calculate the indices of the cubes that intersect with the bounding box
        x_indices = range(int(min_point[0] / voxel_size) - 1, int(max_point[0] / voxel_size) + 1)
        y_indices = range(int(min_point[1] / voxel_size) - 1, int(max_point[1] / voxel_size) + 1)
        z_indices = range(int(min_point[2] / voxel_size) - 1, int(max_point[2] / voxel_size) + 1)

        # Create a list of the centers of all the intersecting voxels
        centers = []
        for x in x_indices:
            for y in y_indices:
                for z in z_indices:
                    center = tuple(round((_ + 1 / 2) * voxel_size, 6) for _ in [x, y, z])
                    centers.append(center)

        return centers

    @staticmethod
    def _voxels_intersecting_voxels(voxel_centers_array: np.ndarray, voxel_size: float) -> Set[Point]:
        """
        Helper method to compute the center of the voxels that intersect with a given array of voxels.
        The returned voxels are part of an implicit 3D grid defined by the voxel size, whereas the given voxels centers
        are not.

        This method is used when translating or rotating the voxels get back to the implicit 3D grid and perform
        efficient Boolean operations (thanks to voxels defined in the same grid).

        :param voxel_centers_array: The array of voxel centers.
        :type voxel_centers_array: numpy.ndarray
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: A set of the centers of the intersecting voxels.
        :rtype: set[tuple[float, float, float]]
        """
        # Compute the voxel indices
        indices = np.floor(voxel_centers_array / voxel_size).astype(int)

        # Compute unique indices to avoid duplicates
        unique_indices = np.unique(indices, axis=0)

        # Convert back to voxel centers
        centers = np.around((unique_indices + 0.5) * voxel_size, 6)

        return set(tuple(center) for center in centers)

    @staticmethod
    def _triangles_to_voxels(triangles: List[Triangle], voxel_size: float) -> Set[Point]:
        """
        Helper method to compute all the voxels intersecting with a given list of triangles.

        :param triangles: The triangles to compute the intersecting voxels.
        :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: The centers of the voxels that intersect with the triangles.
        :rtype: set[tuple[float, float, float]]
        """
        voxel_centers = set()

        for triangle in tqdm(triangles):
            min_point = tuple(min(p[i] for p in triangle) for i in range(3))
            max_point = tuple(max(p[i] for p in triangle) for i in range(3))

            for bbox_center in aabb_intersecting_boxes(min_point, max_point, voxel_size):
                if bbox_center not in voxel_centers:
                    if triangle_intersects_voxel(
                        triangle,
                        bbox_center,
                        [0.5 * voxel_size for _ in range(3)],
                    ):
                        voxel_centers.add(bbox_center)

        return voxel_centers

    @staticmethod
    def _voxel_triangular_faces(voxel_center: Point, voxel_size: float) -> List[Triangle]:
        """
        Helper method to compute the 12 triangular faces that compose a voxel, for visualization.

        :param voxel_center: The voxel center point.
        :type voxel_center: tuple[float, float, float]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: The 12 triangles representing the 6 faces of the given voxel.
        """
        # pylint: disable=invalid-name

        x, y, z = voxel_center
        sx, sy, sz = voxel_size, voxel_size, voxel_size
        hx, hy, hz = sx / 2, sy / 2, sz / 2
        faces = [
            # Front face
            tuple(
                [
                    (x - hx, y - hy, z + hz),
                    (x - hx, y + hy, z + hz),
                    (x + hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z + hz),
                    (x + hx, y + hy, z + hz),
                    (x + hx, y - hy, z + hz),
                ]
            ),
            # Back face
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x - hx, y + hy, z - hz),
                    (x + hx, y + hy, z - hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x + hx, y + hy, z - hz),
                    (x + hx, y - hy, z - hz),
                ]
            ),
            # Left face
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x - hx, y + hy, z - hz),
                    (x - hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x - hx, y + hy, z + hz),
                    (x - hx, y - hy, z + hz),
                ]
            ),
            # Right face
            tuple(
                [
                    (x + hx, y - hy, z - hz),
                    (x + hx, y + hy, z - hz),
                    (x + hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x + hx, y - hy, z - hz),
                    (x + hx, y + hy, z + hz),
                    (x + hx, y - hy, z + hz),
                ]
            ),
            # Top face
            tuple(
                [
                    (x - hx, y + hy, z - hz),
                    (x + hx, y + hy, z - hz),
                    (x + hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y + hy, z - hz),
                    (x + hx, y + hy, z + hz),
                    (x - hx, y + hy, z + hz),
                ]
            ),
            # Bottom face
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x + hx, y - hy, z - hz),
                    (x + hx, y - hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x + hx, y - hy, z + hz),
                    (x - hx, y - hy, z + hz),
                ]
            ),
        ]

        triangular_faces = [tuple(tuple(round(_, 6) for _ in point) for point in face) for face in faces]
        return triangular_faces

    @staticmethod
    def _triangles_to_segments(
        triangles: Set[Triangle],
    ) -> Tuple[Dict[float, Set[Segment]], Dict[float, Set[Segment]], Dict[float, Set[Segment]]]:
        """
        Helper method to extract the segments of a given list of triangle representing a voxelization.
        The segments are sorted by plane: a plane is defined by its normal vector (X, Y or Z) and abscissa.

        Only the segments representing the relevant contours of a faces are returned.
        These relevant segments are the one that are only present in the list of triangles, because if a segment is
        present twice, this means that it's inside the face and is not interesting for defining the face contour.

        :param triangles: The triangles defining the voxelization.
        :type triangles: Set[Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]]

        :return: The extracted relevant segments, sorted by plane.
        :rtype: tuple[dict[float, set[Segment]], dict[float, set[Segment]], dict[float, set[Segment]]]
        """
        segments_x = {}
        segments_y = {}
        segments_z = {}

        for triangle in triangles:
            if triangle[0][0] == triangle[1][0] == triangle[2][0]:
                segments = segments_x.get(triangle[0][0], set())
                triangle_segments = [
                    (tuple(triangle[0][1:]), tuple(triangle[1][1:])),
                    (tuple(triangle[1][1:]), tuple(triangle[2][1:])),
                    (tuple(triangle[2][1:]), tuple(triangle[0][1:])),
                ]

                segments_x[triangle[0][0]] = Voxelization._add_triangle_segments_to_segments_set(
                    triangle_segments, segments
                )

            elif triangle[0][1] == triangle[1][1] == triangle[2][1]:
                segments = segments_y.get(triangle[0][1], set())
                triangle_segments = [
                    (
                        tuple([triangle[0][0], triangle[0][2]]),
                        tuple([triangle[1][0], triangle[1][2]]),
                    ),
                    (
                        tuple([triangle[1][0], triangle[1][2]]),
                        tuple([triangle[2][0], triangle[2][2]]),
                    ),
                    (
                        tuple([triangle[2][0], triangle[2][2]]),
                        tuple([triangle[0][0], triangle[0][2]]),
                    ),
                ]

                segments_y[triangle[0][1]] = Voxelization._add_triangle_segments_to_segments_set(
                    triangle_segments, segments
                )

            else:
                segments = segments_z.get(triangle[0][2], set())
                triangle_segments = [
                    (tuple(triangle[0][:2]), tuple(triangle[1][:2])),
                    (tuple(triangle[1][:2]), tuple(triangle[2][:2])),
                    (tuple(triangle[2][:2]), tuple(triangle[0][:2])),
                ]

                segments_z[triangle[0][2]] = Voxelization._add_triangle_segments_to_segments_set(
                    triangle_segments, segments
                )

        return segments_x, segments_y, segments_z

    @staticmethod
    def _add_triangle_segments_to_segments_set(
        triangle_segments: List[Segment], segments_set: Set[Segment]
    ) -> Set[Segment]:
        """
        Helper method to add the segments defining a triangle to a set of segments, avoiding non-relevant segments
        (i.e. segments that does not represent an inner or outer contour of the face, but that are inside the face).

        :param triangle_segments: The segments defining the triangle.
        :type triangle_segments: list[tuple[tuple[float, float], tuple[float, float]]]
        :param segments_set: The set in which the segments are added.
        :type segments_set: set[tuple[tuple[float, float], tuple[float, float]]]

        :return: The edited set of segments.
        :rtype: set[tuple[tuple[float, float], tuple[float, float]]]
        """
        for segment in triangle_segments:
            if segment not in segments_set and segment[::-1] not in segments_set:
                segments_set.add(segment)
            else:
                if segment in segments_set:
                    segments_set.remove(segment)
                else:
                    segments_set.remove(segment[::-1])

        return segments_set

    @staticmethod
    def _order_segments(segments: Iterable[Segment]) -> List[List[Segment]]:
        """
        Helper method to order a given set of segments defined in the same plane.
        The segments may define several contours, so this method define a list of lists of ordered segments,
        representing the multiple contours with ordered segments.

        :param segments: The segments to order.
        :type segments: Iterable[tuple[tuple[float, float], tuple[float, float]]]

        :return: The ordered segments lists.
        :rtype: list[list[tuple[tuple[float, float], tuple[float, float]]]]
        """
        segments = list(segments)
        ordered_segments_list = []

        while segments:
            # Start a new contour with the first segment
            ordered_segments = [segments.pop()]
            ordered_segments_list.append(ordered_segments)

            while True:
                last_segment = ordered_segments[-1]

                # Find the next segment
                for i, segment in enumerate(segments):
                    if segment[0] == last_segment[1]:
                        ordered_segments.append(segments.pop(i))
                        break

                    if segment[1] == last_segment[1]:
                        # If the segment is in the opposite direction, reverse it
                        ordered_segments.append((segment[1], segment[0]))
                        segments.pop(i)
                        break
                else:
                    # No more segments connect to this contour, start a new one
                    break

        return ordered_segments_list

    @staticmethod
    def _triangles_to_closed_polygons(triangles: Set[Triangle]) -> List[Dict[float, List[ClosedPolygon2D]]]:
        """
        Helper method to convert the triangles defined the faces of the voxelization to ClosedPolygon2D, sorted
        by plane (a plane is defined by its normal vector (X, Y or Z) and its abscissa).

        :param triangles: The triangles defining the voxelization.
        :type triangles: Set[Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]]

        :return: All the ClosedPolygon2D defining the faces of the voxelization, sorted by plane.
        :rtype: list[dict[float, list[ClosedPolygon2D]]]
        """
        segments_x, segments_y, segments_z = Voxelization._triangles_to_segments(triangles)
        polygons = [{}, {}, {}]

        for i, segments_i in enumerate([segments_x, segments_y, segments_z]):
            for abscissa, segments in segments_i.items():
                # Order the segments lists
                ordered_segments_list = Voxelization._order_segments(segments)

                # Split the self-intersecting polygons
                split_ordered_segments_list = []
                for ordered_segments in ordered_segments_list:
                    split_ordered_segments_list.extend(Voxelization._split_ordered_segments(ordered_segments))

                polygons[i][abscissa] = [
                    Voxelization._ordered_segments_to_closed_polygon_2d(split_ordered_segments)
                    for split_ordered_segments in split_ordered_segments_list
                ]

        return polygons

    @staticmethod
    def _split_ordered_segments(ordered_segments: List[Segment]) -> List[List[Segment]]:
        """
        Helper method to split the self-crossing polygons into multiple polygons to avoid triangulation errors.

        :param ordered_segments: The ordered segments defining the polygon.
        :type ordered_segments: list[tuple[tuple[float, float], tuple[float, float]]]

        :return: The split polygon into multiple polygons, defined by segment list.
        :rtype: list[list[tuple[tuple[float, float], tuple[float, float]]]]
        """
        # TODO: refactor this method to make it more understandable
        # pylint: disable=too-many-locals,too-many-branches

        points = [segment[0] for segment in ordered_segments]

        if len(points) == len(set(points)):
            # If no duplicate points, the polygon is not self-intersecting
            return [ordered_segments]

        graph = Graph(ordered_segments)
        ordered_segments_list = []

        for cycle in graph.find_cycles():
            ordered_segments_list.extend(Voxelization._order_segments(cycle))
            # print(Voxelization._order_segments(cycle))

        polygons = {}
        for _ordered_segments in ordered_segments_list:
            polygons[tuple(_ordered_segments)] = Voxelization._ordered_segments_to_closed_polygon_2d(_ordered_segments)

        # Check for polygon inclusion
        children = {}
        parents = {}
        for polygon in polygons.values():
            children[polygon] = []
            parents[polygon] = []

        for polygon_1 in polygons.values():
            for polygon_2 in polygons.values():
                if polygon_1 != polygon_2:
                    if polygon_1.is_inside(polygon_2):
                        children[polygon_1].append(polygon_2)
                        parents[polygon_2].append(polygon_1)
                    elif polygon_2.is_inside(polygon_1):
                        children[polygon_2].append(polygon_1)
                        parents[polygon_1].append(polygon_2)

        # Actual split polygons have no children in the cycle-defined polygons from the self-intersecting polygon
        candidate_polygons = {}
        for _ordered_segments, polygon in polygons.items():
            if len(children[polygon]) == 0:
                candidate_polygons[_ordered_segments] = polygon

        # If the polygon has all its segments part of the other polygons segments, it is an inner polygon
        ordered_segments_list = []
        for _ordered_segments in candidate_polygons:
            ordered_segments_list.append(list(_ordered_segments))

        ordered_segments_list_valid = []
        for i, _ordered_segments in enumerate(ordered_segments_list):
            all_current_polygon_segments = []
            all_other_polygons_segments = []

            for segments in _ordered_segments:
                all_current_polygon_segments.append(segments)
                all_current_polygon_segments.append(segments[::-1])

            for segments in ordered_segments_list[:i] + ordered_segments_list[i + 1 :]:
                for segment in segments:
                    all_other_polygons_segments.append(segment)
                    all_other_polygons_segments.append(segment[::-1])

            for segment in all_current_polygon_segments:
                if segment not in all_other_polygons_segments:
                    ordered_segments_list_valid.append(_ordered_segments)
                    break

        return ordered_segments_list_valid

    @staticmethod
    def _simplify_polygon_points(points: List[Point]) -> List[Point]:
        """
        Helper method to simplify the list of points defining a polygon, by removing the points that are not relevant
        for the polygon definition, i.e. the points that are not in a corner of the polygon.

        :param points: The points defining the polygon.
        :type points: list[tuple[float, float, float]]

        :return: The simplified points.
        :rtype: list[tuple[float, float, float]]
        """
        simplified_point = []

        for i in range(len(points) - 1):
            if not (
                points[i - 1][0] == points[i][0] == points[i + 1][0]
                or points[i - 1][1] == points[i][1] == points[i + 1][1]
            ):
                simplified_point.append(points[i])

        if not (points[-2][0] == points[-1][0] == points[0][0] or points[-2][1] == points[-1][1] == points[0][1]):
            simplified_point.append(points[-1])

        return simplified_point

    @staticmethod
    def _ordered_segments_to_closed_polygon_2d(ordered_segments: Iterable[Segment]) -> ClosedPolygon2D:
        """
        Helper method to convert an iterable of ordered segments to a ClosedPolygon2D.

        :param ordered_segments: The segments that compose the ClosedPolygon2D.
        :type ordered_segments: Iterable[tuple[tuple[float, float], tuple[float, float]]]

        :return: The created ClosedPolygon2D.
        :rtype: ClosedPolygon2D
        """
        return ClosedPolygon2D(
            points=[
                Point2D(*point)
                for point in Voxelization._simplify_polygon_points([segment[0] for segment in ordered_segments])
            ]
        )

    def to_triangles(self) -> Set[Triangle]:
        """
        Convert the voxelization to triangles for display purpose.

        Only the relevant faces are returned (i.e. the faces that are not at the interface of two different voxel,
        i.e. the faces that are only present once in the list of triangles representing the triangulated voxels).

        :return: The triangles representing the voxelization.
        :rtype: set[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangles = set()

        for voxel in self.voxels_centers:
            for triangle in self._voxel_triangular_faces(voxel, self.voxel_size):
                if triangle not in triangles:
                    triangles.add(triangle)
                else:
                    triangles.remove(triangle)

        return triangles

    def to_closed_shell(self) -> ClosedShell3D:
        """
        Convert the voxelization to a ClosedShell3D for display purpose.

        This ClosedShell3D is only composed of PlaneFace3D that represent the faces of the voxelized geometry, and
        are made to be the lightest geometrical representation of the voxelization (in terms of number of triangles
        in the babylonjs model).

        :return: The created ClosedShell3D with only PlaneFace3D.
        :rtype: ClosedShell3D
        """
        # TODO: refactor this method to make it more understandable
        # pylint: disable=too-many-nested-blocks

        polygons = self._triangles_to_closed_polygons(self.to_triangles())
        planes = [PLANE3D_OYZ, PLANE3D_OXZ, PLANE3D_OXY]
        faces = []

        for i, polygons_i in enumerate(polygons):
            for abscissa, polygons in polygons_i.items():
                offset = [0, 0, 0]
                offset[i] = abscissa
                offset = Vector3D(*offset)

                plane = planes[i].translation(offset)

                if len(polygons) == 1:
                    faces.append(PlaneFace3D(plane, Surface2D(polygons[0], [])))
                else:
                    # Check for polygon inclusion
                    children = {}
                    parents = {}
                    for polygon in polygons:
                        children[polygon] = []
                        parents[polygon] = []

                    for polygon_1 in polygons:
                        for polygon_2 in polygons:
                            if polygon_1 != polygon_2:
                                if polygon_1.is_inside(polygon_2):
                                    children[polygon_1].append(polygon_2)
                                    parents[polygon_2].append(polygon_1)
                                elif polygon_2.is_inside(polygon_1):
                                    children[polygon_2].append(polygon_1)
                                    parents[polygon_1].append(polygon_2)

                    for polygon in polygons:
                        if len(parents[polygon]) // 2 == 0:
                            inner_polygons = [
                                child
                                for child in children[polygon]
                                if set(children[polygon]).intersection(children[child]) == set()
                            ]
                            faces.append(PlaneFace3D(plane, Surface2D(polygon, inner_polygons)))

        return ClosedShell3D(faces, name=self.name)

    def to_closed_triangle_shell(self) -> ClosedTriangleShell3D:
        """
        Convert the voxelization to a ClosedTriangleShell3D for display purpose.

        To create the ClosedTriangleShell3D, this method triangulates each voxel and convert it to Triangle3D.
        The method is robust and fast but creates over-triangulated geometry.

        :return: The created ClosedShell3D with only Triangle3D.
        :rtype: ClosedTriangleShell3D
        """
        triangles3d = [
            Triangle3D(Point3D(*triangle[0]), Point3D(*triangle[1]), Point3D(*triangle[2]))
            for triangle in self.to_triangles()
        ]
        shell = ClosedTriangleShell3D(triangles3d, name=self.name)

        return shell

    def volmdlr_primitives(self, **kwargs):
        """
        Return a list of volmdlr primitives to build up volume model.
        It uses the simplified representation of the voxelization given by the "to_closed_shell" method.
        """
        if len(self.voxels_centers) == 0:
            warnings.warn("Empty voxelization.")
            return []

        return [self.to_closed_triangle_shell()]

    def intersection(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Create a voxelization that is the Boolean intersection of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean intersection with.
        :type other_voxelization: Voxelization

        :return: The created voxelization resulting from the Boolean intersection.
        :rtype: Voxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform intersection.")

        return Voxelization(self.voxels_centers.intersection(other_voxelization.voxels_centers), self.voxel_size)

    def is_intersecting(self, other_voxelization: "Voxelization") -> bool:
        """
        Check if two voxelizations are intersecting.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to check if there is an intersection with.
        :type other_voxelization: Voxelization

        :return: True if the voxelizations are intersecting, False otherwise.
        :rtype: bool
        """
        intersection = self.intersection(other_voxelization)

        return len(intersection) > 0

    def union(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Create a voxelization that is the Boolean union of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean union with.
        :type other_voxelization: Voxelization

        :return: The created voxelization resulting from the Boolean union.
        :rtype: Voxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform union.")

        return Voxelization(self.voxels_centers.union(other_voxelization.voxels_centers), self.voxel_size)

    def difference(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Create a voxelization that is the Boolean difference of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean difference with.
        :type other_voxelization: Voxelization

        :return: The created voxelization resulting from the Boolean difference.
        :rtype: Voxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform difference.")

        return Voxelization(self.voxels_centers.difference(other_voxelization.voxels_centers), self.voxel_size)

    def symmetric_difference(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Create a voxelization that is the Boolean symmetric difference (XOR) of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean symmetric difference with.
        :type other_voxelization: Voxelization

        :return: The created voxelization resulting from the Boolean symmetric difference.
        :rtype: Voxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform symmetric difference.")

        return Voxelization(
            self.voxels_centers.symmetric_difference(other_voxelization.voxels_centers), self.voxel_size
        )

    def interference(self, other_voxelization: "Voxelization") -> float:
        """
        Compute the percentage of interference between two voxelization.

        :param other_voxelization: The other voxelization to compute percentage of interference with.
        :type other_voxelization: Voxelization

        :return: The percentage of interference between the two voxelization.
        :rtype: float
        """
        return len(self.intersection(other_voxelization)) / len(self.union(other_voxelization))

    @staticmethod
    def _rotation_matrix(axis: Vector3D, angle: float) -> np.array:
        """
        Helper method that compute a rotation matrix from an axis and a radians angle.

        :param axis: The rotation axis.
        :type axis: Vector3D
        :param angle: The rotation angle.
        :type angle: float

        :return: The computed rotation matrix.
        :rtype: numpy.array
        """
        # pylint: disable=invalid-name,too-many-locals
        axis = np.array([axis.x, axis.y, axis.z])

        axis = axis / np.linalg.norm(axis)
        a = np.cos(angle / 2.0)

        b, c, d = -axis * np.sin(angle / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

        return np.array(
            [
                [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
            ]
        )

    def rotation(self, center: Point3D, axis: Vector3D, angle: float):
        """
        Rotate the voxelization around the specified center, axis, and angle.

        :param center: The center point of rotation.
        :type center: Point3D
        :param axis: The rotation axis.
        :type axis: Vector3D
        :param angle: The rotation angle in radians.
        :type angle: float

        :return: A new Voxelization object resulting from the rotation.
        :rtype: Voxelization
        """
        rotation_matrix = self._rotation_matrix(axis, angle)
        voxel_array = np.array(list(self.voxels_centers)) - np.array([center.x, center.y, center.z])
        rotated_voxels = np.dot(voxel_array, rotation_matrix.T)
        rotated_voxels += np.array([center.x, center.y, center.z])

        intersecting_voxels = self._voxels_intersecting_voxels(rotated_voxels, self.voxel_size)

        return Voxelization(intersecting_voxels, self.voxel_size)

    def translation(self, offset: Vector3D):
        """
        Translate the voxelization by the specified offset.

        :param offset: The translation offset.
        :type offset: Vector3D

        :return: A new Voxelization object resulting from the translation.
        :rtype: Voxelization
        """
        voxel_array = np.array(list(self.voxels_centers))
        translated_voxels = voxel_array + np.array([offset.x, offset.y, offset.z])

        intersecting_voxels = self._voxels_intersecting_voxels(translated_voxels, self.voxel_size)

        return Voxelization(intersecting_voxels, self.voxel_size)

    @staticmethod
    def voxel_center_in_implicit_grid(voxel_center: Point, voxel_size: float) -> bool:
        """
        Check if a given voxel center point is a voxel center of the implicit grid, defined by voxel_size.

        :param voxel_center: The voxel center point to check.
        :type voxel_center: tuple[float, float, float]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: True if the given voxel center point is a voxel center of the implicit grid, False otherwise.
        :rtype: bool
        """
        for coord in voxel_center:
            if not round((coord - 0.5 * voxel_size) / voxel_size, 6).is_integer():
                return False

        return True

    @classmethod
    def from_voxel_matrix(
        cls, voxel_matrix: List[List[List[bool]]], voxel_size: float, voxel_matrix_origin_center: Point
    ):
        """
        Create a Voxelization object from a voxel matrix.

        :param voxel_matrix: The voxel matrix representing the voxelization.
        :type voxel_matrix: list[list[list[bool]]]
        :param voxel_size: The size of the voxel edges.
        :type voxel_size: float
        :param voxel_matrix_origin_center: Voxel center of the origin of the voxel matrix, i.e 'voxel_matrix[0][0][0]'.
        :type voxel_matrix_origin_center: tuple[float, float, float]

        :return: A Voxelization object created from the voxel matrix.
        :rtype: Voxelization
        """
        if not cls.voxel_center_in_implicit_grid(voxel_matrix_origin_center, voxel_size):
            warnings.warn(
                """This voxel matrix is not defined in the implicit grid defined by the voxel_size. 
            Some methods like boolean operation or interference computing may not work as expected."""
            )

        voxels_centers = set()
        for i, row in enumerate(voxel_matrix):
            for j, col in enumerate(row):
                for k, voxel in enumerate(col):
                    if voxel:
                        center = (
                            round(voxel_matrix_origin_center[0] + i * voxel_size, 6),
                            round(voxel_matrix_origin_center[1] + j * voxel_size, 6),
                            round(voxel_matrix_origin_center[2] + k * voxel_size, 6),
                        )
                        voxels_centers.add(center)

        return cls(voxels_centers, voxel_size)

    def get_min_voxel_grid_center(self) -> Point:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.
        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, float, float]
        """
        min_x = min(point[0] for point in self.voxels_centers)
        min_y = min(point[1] for point in self.voxels_centers)
        min_z = min(point[2] for point in self.voxels_centers)

        return min_x, min_y, min_z

    def get_max_voxel_grid_center(self) -> Point:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.
        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, float, float]
        """
        max_x = max(point[0] for point in self.voxels_centers)
        max_y = max(point[1] for point in self.voxels_centers)
        max_z = max(point[2] for point in self.voxels_centers)

        return max_x, max_y, max_z

    def to_voxel_matrix(self) -> List[List[List[bool]]]:
        """
        Convert the voxelization to a voxel matrix.

        :return: The voxel matrix representing the voxelization.
        :rtype: list[list[list[bool]]]
        """
        min_center = self.get_min_voxel_grid_center()
        max_center = self.get_max_voxel_grid_center()

        dim_x = round((max_center[0] - min_center[0]) / self.voxel_size + 1)
        dim_y = round((max_center[1] - min_center[1]) / self.voxel_size + 1)
        dim_z = round((max_center[2] - min_center[2]) / self.voxel_size + 1)

        matrix = [[[False for _ in range(dim_z)] for _ in range(dim_y)] for _ in range(dim_x)]

        for voxel_center in self.voxels_centers:
            x = round((voxel_center[0] - min_center[0]) / self.voxel_size)
            y = round((voxel_center[1] - min_center[1]) / self.voxel_size)
            z = round((voxel_center[2] - min_center[2]) / self.voxel_size)

            matrix[x][y][z] = True

        return matrix

    def inverse(self) -> "Voxelization":
        """
        Create a new Voxelization object that is the inverse of the current voxelization.

        :return: The inverse Voxelization object.
        :rtype: Voxelization
        """
        voxel_matrix = self.to_voxel_matrix()
        inverted_matrix = np.logical_not(np.array(voxel_matrix)).tolist()
        min_voxel_center = self.get_min_voxel_grid_center()

        return Voxelization.from_voxel_matrix(inverted_matrix, self.voxel_size, min_voxel_center)


class OctreeNode:
    """Class representing an octree node for octree voxelization purpose."""

    def __init__(self, center: Point, size: float, depth: int, max_depth: int):
        self.children = []
        self.center = center
        self.size = size
        self.depth = depth
        self.max_depth = max_depth

    @classmethod
    def octree_voxelization_depth_based(cls, triangles: List[Triangle], max_depth: int) -> "OctreeNode":
        """Create a voxelization based on the depth of the tree."""
        # Compute the size of the bounding cube (root voxel)
        min_corner = np.min([np.min(triangle, axis=0) for triangle in triangles], axis=0)
        max_corner = np.max([np.max(triangle, axis=0) for triangle in triangles], axis=0)
        root_size = max(max_corner - min_corner)

        # Compute the center of the bounding cube (root voxel)
        corners = np.stack([min_corner, max_corner])
        center = corners.mean(axis=0)

        root = cls(center, root_size, 0, max_depth)
        root.subdivide(triangles)

        return root

    @classmethod
    def octree_voxelization_size_based(cls, triangles: List[Triangle], voxel_size: float) -> "OctreeNode":
        """Create a voxelization based on the size of the voxel."""
        min_corner = np.min([np.min(triangle, axis=0) for triangle in triangles], axis=0)
        max_corner = np.max([np.max(triangle, axis=0) for triangle in triangles], axis=0)

        # Compute the corners in the implicit grid defined by the voxel size
        min_corner = (min_corner // voxel_size - 2) * voxel_size
        max_corner = (max_corner // voxel_size + 2) * voxel_size

        root_size = max(max_corner - min_corner)

        # Compute the max depth corresponding the voxel_size
        max_depth = math.ceil(math.log2(root_size // voxel_size))

        # Compute the max corner to have voxel of given voxel size with the octree process
        max_corner = min_corner + ((2**max_depth) * voxel_size)
        root_size = max(max_corner - min_corner)

        corners = np.stack([min_corner, max_corner])
        center = np.round(corners.mean(axis=0), 6)

        root = cls(center, root_size, 0, max_depth)
        root.subdivide(triangles)

        return root

    def subdivide(self, triangles: List[Triangle]) -> None:
        """Recursive method to subdivide the voxelization in 8 until the wanted tree depth is reached."""
        children = []
        if self.depth < self.max_depth:
            half_size = round(self.size / 2, 6)
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        # calculate the center of the child node
                        child_center = (
                            np.round(self.center[0] + (i - 0.5) * half_size, 6),
                            np.round(self.center[1] + (j - 0.5) * half_size, 6),
                            np.round(self.center[2] + (k - 0.5) * half_size, 6),
                        )

                        # create a new OctreeNode for the child
                        child_node = OctreeNode(
                            child_center,
                            half_size,
                            self.depth + 1,
                            self.max_depth,
                        )

                        if self.depth != 0:
                            # add the child node to the list of children
                            # if the child node intersects with the triangles
                            self.children.extend(self.process_node(child_node, triangles))
                        else:
                            children.append(child_node)

            if self.depth == 0:
                for node in tqdm(self.process_node(node, triangles) for node in children):
                    self.children.extend(node)

        else:
            return  # reached max depth, do not subdivide further.

    def intersecting_triangles(self, triangles: List[Triangle]) -> List[Triangle]:
        """Given a list of triangle, return the triangles that are intersecting with the node."""
        intersecting_triangles = [
            triangle
            for triangle in triangles
            if triangle_intersects_voxel(triangle, self.center, [0.5 * self.size for _ in range(3)])
        ]
        return intersecting_triangles

    def get_leaf_centers(self) -> List[Point]:
        """Recursive method to extract all the leaf voxel center (voxels of minimal size)."""
        if self.depth == self.max_depth:  # if max depth reached, it is a leaf node
            return [self.center]

        centers = []
        for child in self.children:
            centers += child.get_leaf_centers()

        return centers

    @staticmethod
    def process_node(node, triangles):
        """
        Helper method to process a node .
        In other terms, adding it to the tree and recursively divide it if it intersects with the geometry.
        """
        intersecting_triangles = node.intersecting_triangles(triangles)
        if len(intersecting_triangles) > 0:
            # recursively subdivide the child node
            node.subdivide(intersecting_triangles)
            return [node]
        return []


class Graph:
    """Helper class for subdividing self-crossing polygons using graph theory algorithm."""

    def __init__(self, edges):
        """Initialize the graph from edges."""
        self.graph = self.build_graph(edges)
        self.visited = {node: False for node in self.graph}

    @staticmethod
    def build_graph(edges):
        """Builds adjacency list representation of graph from edge list."""
        graph = {}
        for edge in edges:
            if edge[0] not in graph:
                graph[edge[0]] = []
            if edge[1] not in graph:
                graph[edge[1]] = []
            graph[edge[0]].append(edge[1])
            graph[edge[1]].append(edge[0])  # assuming undirected graph
        return graph

    def find_cycles(self):
        """Finds all unique cycles in the graph."""
        cycles = []
        for node in self.graph:
            self._dfs(node, node, [], cycles)
        return self._convert_to_edge_cycles(cycles)

    def _dfs(self, node, start, path, cycles):
        """Performs depth-first search to find cycles."""
        self.visited[node] = True
        path.append(node)
        for neighbor in self.graph[node]:
            if not self.visited[neighbor]:
                self._dfs(neighbor, start, path, cycles)
            elif neighbor == start and len(path) >= 3:
                cycles.append(path[:])
        path.pop()
        self.visited[node] = False

    @staticmethod
    def _convert_to_edge_cycles(cycles):
        """Converts node-based cycles to edge-based cycles, removing duplicates."""
        edge_cycles_set = set()
        for cycle in cycles:
            edge_cycle = []
            for i, node in enumerate(cycle):
                edge_cycle.append(tuple(sorted([node, cycle[(i + 1) % len(cycle)]])))
            edge_cycle.sort()  # sort the edges in the cycle
            edge_cycles_set.add(tuple(edge_cycle))  # add it to the set

        return [list(cycle) for cycle in edge_cycles_set]
