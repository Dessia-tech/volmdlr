"""
Class for voxel representation of volmdlr models
"""
import warnings
from copy import copy
from typing import List, Set, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from dessia_common.core import PhysicalObject

from volmdlr import Point2D, Point3D, Vector3D
from volmdlr.core import BoundingBox, BoundingRectangle, VolumeModel
from volmdlr.faces import Triangle3D
from volmdlr.shells import ClosedTriangleShell3D, Shell3D
from volmdlr.voxelization_compiled import (
    flood_fill_matrix_2d,
    flood_fill_matrix_3d,
    line_segments_to_pixels,
    triangles_to_voxels,
    triangles_to_voxel_matrix,
)

# Custom types
Point = Tuple[float, float, float]
Triangle = Tuple[Point, Point, Point]
Segment = Tuple[Point, Point]


class PointVoxelization(PhysicalObject):
    """
    Class for creation and manipulation of voxelization of volmdlr geometry.

    This approach is used to create a Voxelization of the surfaces, and do not fill the volume.

    The voxelization is defined in an implicit grid of the 3D space. The grid is defined by the size of the voxels.
    The implicit grid is defined by discretizing the 3D space into axis aligned cubes with given size, with the
    origin of the global landmark (0, 0, 0) always being a corner point (and not inside a cube).

    For example, in 1D and with a voxel size of t, the set of voxels (defined by minimal and maximal point) can only
    be defined as {i ∈ N, t ∈ R | (i * t, (i+1) * t)}
    The corresponding set of voxel centers is defined by the following set: {i ∈ N, t ∈ R | (i + 0.5) * t}

    This approach allows to always define voxelization of same size in the same grid of the 3D space, which is very
    useful to perform very fast Boolean operations.
    """

    def __init__(self, voxel_centers: Set[Point], voxel_size: float, name: str = ""):
        """
        Initialize the Voxelization.

        :param voxel_centers: The set of points representing voxel centers.
        :type voxel_centers: set[tuple[float, float, float]]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the Voxelization.
        :type name: str, optional
        """
        self.voxel_centers = voxel_centers
        self.voxel_size = voxel_size

        PhysicalObject.__init__(self, name=name)

    def __eq__(self, other_voxelization: "PointVoxelization") -> bool:
        """
        Check if the current voxelization is equal to another voxelization.

        :param other_voxelization: The voxelization to compare.
        :type other_voxelization: PointVoxelization

        :return: True if the voxelizations are equal, False otherwise.
        :rtype: bool
        """
        return (
            self.voxel_centers == other_voxelization.voxel_centers and self.voxel_size == other_voxelization.voxel_size
        )

    def __add__(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Return the union of the current voxelization with another voxelization.

        :param other_voxelization: The voxelization to union with.
        :type other_voxelization: PointVoxelization

        :return: The union of the voxelizations.
        :rtype: PointVoxelization
        """
        return self.union(other_voxelization)

    def __sub__(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Return the difference between the current voxelization and another voxelization.

        :param other_voxelization: The voxelization to subtract.
        :type other_voxelization: PointVoxelization

        :return: The difference between the voxelizations.
        :rtype: PointVoxelization
        """
        return self.difference(other_voxelization)

    def __and__(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Return the intersection of the current voxelization with another voxelization.

        :param other_voxelization: The voxelization to intersect with.
        :type other_voxelization: PointVoxelization

        :return: The intersection of the voxelizations.
        :rtype: PointVoxelization
        """
        return self.intersection(other_voxelization)

    def __or__(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Return the union of the current voxelization with another voxelization.

        :param other_voxelization: The voxelization to union with.
        :type other_voxelization: PointVoxelization

        :return: The union of the voxelizations.
        :rtype: PointVoxelization
        """
        return self.union(other_voxelization)

    def __xor__(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Return the symmetric difference between the current voxelization and another voxelization.

        :param other_voxelization: The voxelization to calculate the symmetric difference with.
        :type other_voxelization: PointVoxelization
        :return: The symmetric difference between the voxelizations.
        :rtype: PointVoxelization
        """
        return self.symmetric_difference(other_voxelization)

    def __invert__(self) -> "PointVoxelization":
        """
        Return the inverse of the current voxelization.

        :return: The inverse Voxelization object.
        :rtype: PointVoxelization
        """
        return self.inverse()

    def __len__(self):
        """
        Return the number of voxels in the voxelization.

        :return: The number of voxels.
        :rtype: int
        """
        return len(self.voxel_centers)

    @classmethod
    def from_shell(
        cls, shell: Shell3D, voxel_size: float, name: str = ""
    ) -> "PointVoxelization":
        """
        Create a Voxelization from a ClosedTriangleShell3D.

        :param shell: The Shell3D to voxelize.
        :type shell: Shell3D
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: PointVoxelization
        """
        triangles = _shell_to_triangles(shell)
        voxels = cls._triangles_to_voxels(triangles, voxel_size)

        return cls(voxel_centers=voxels, voxel_size=voxel_size, name=name)

    @classmethod
    def from_volume_model(
        cls, volume_model: VolumeModel, voxel_size: float, name: str = ""
    ) -> "PointVoxelization":
        """
        Create a Voxelization from a VolumeModel.

        :param volume_model: The VolumeModel to voxelize.
        :type volume_model: VolumeModel
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: PointVoxelization
        """
        triangles = _volume_model_to_triangles(volume_model)
        voxels = cls._triangles_to_voxels(triangles, voxel_size)

        return cls(voxel_centers=voxels, voxel_size=voxel_size, name=name)

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
        return triangles_to_voxels(triangles, voxel_size).union(
            PointVoxelization._triangles_to_voxels_at_interface(triangles, voxel_size)
        )

    @staticmethod
    def _triangles_to_voxels_at_interface(triangles: List[Triangle], voxel_size: float) -> Set[Point]:
        voxel_centers = set()
        for triangle in triangles:
            voxel_centers = voxel_centers.union(PointVoxelization._triangle_interface_voxels(triangle, voxel_size))

        return voxel_centers

    @staticmethod
    def _triangle_interface_voxels(triangle: Triangle, voxel_size: float) -> Set[Point]:
        """
        Calculate the set of voxel centers that intersect the interface of a triangle within the voxelization.

        :param triangle: The triangle to calculate the voxel centers for.
        :type triangle: Triangle
        :param voxel_size: The size of the voxels in the voxelization.
        :type voxel_size: float
        :return: The set of voxel centers that intersect the triangle interface.
        :rtype: Set[Point]
        """
        voxel_centers = set()

        for i in range(3):
            # Check if the triangle is defined in the XY, XZ or YZ plane
            if round(triangle[0][i], 6) == round(triangle[1][i], 6) == round(triangle[2][i], 6):
                abscissa = round(triangle[0][i], 6)

                # Check if this plane is defined is at the interface between voxels
                if round(abscissa / voxel_size, 6).is_integer():
                    # Define the 3D triangle in 2D

                    v0 = triangle[0][:i] + triangle[0][i + 1 :]
                    v1 = triangle[1][:i] + triangle[1][i + 1 :]
                    v2 = triangle[2][:i] + triangle[2][i + 1 :]

                    triangle_2d = np.array([v0, v1, v2])

                    pixelization = Pixelization.from_polygon(triangle_2d, voxel_size).fill_enclosed_pixels()

                    for center in pixelization.pixel_centers:
                        center_left = list(copy(center))
                        center_left.insert(i, round(abscissa - voxel_size / 2, 6))
                        voxel_centers.add(tuple(center_left))

                        center_right = list(copy(center))
                        center_right.insert(i, round(abscissa + voxel_size / 2, 6))
                        voxel_centers.add(tuple(center_right))

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

    def to_triangles(self) -> Set[Triangle]:
        """
        Convert the voxelization to triangles for display purpose.

        Only the relevant faces are returned (i.e. the faces that are not at the interface of two different voxel,
        i.e. the faces that are only present once in the list of triangles representing the triangulated voxels).

        :return: The triangles representing the voxelization.
        :rtype: set[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangles = set()

        for voxel in self.voxel_centers:
            for triangle in self._voxel_triangular_faces(voxel, self.voxel_size):
                if triangle not in triangles:
                    triangles.add(triangle)
                else:
                    triangles.remove(triangle)

        return triangles

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
        if len(self.voxel_centers) == 0:
            raise ValueError("Empty voxelization.")

        return [self.to_closed_triangle_shell()]

    def intersection(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Create a voxelization that is the Boolean intersection of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean intersection with.
        :type other_voxelization: PointVoxelization

        :return: The created voxelization resulting from the Boolean intersection.
        :rtype: PointVoxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform intersection.")

        return PointVoxelization(self.voxel_centers.intersection(other_voxelization.voxel_centers), self.voxel_size)

    def is_intersecting(self, other_voxelization: "PointVoxelization") -> bool:
        """
        Check if two voxelizations are intersecting.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to check if there is an intersection with.
        :type other_voxelization: PointVoxelization

        :return: True if the voxelizations are intersecting, False otherwise.
        :rtype: bool
        """
        intersection = self.intersection(other_voxelization)

        return len(intersection) > 0

    def union(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Create a voxelization that is the Boolean union of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean union with.
        :type other_voxelization: PointVoxelization

        :return: The created voxelization resulting from the Boolean union.
        :rtype: PointVoxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform union.")

        return PointVoxelization(self.voxel_centers.union(other_voxelization.voxel_centers), self.voxel_size)

    def difference(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Create a voxelization that is the Boolean difference of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean difference with.
        :type other_voxelization: PointVoxelization

        :return: The created voxelization resulting from the Boolean difference.
        :rtype: PointVoxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform difference.")

        return PointVoxelization(self.voxel_centers.difference(other_voxelization.voxel_centers), self.voxel_size)

    def symmetric_difference(self, other_voxelization: "PointVoxelization") -> "PointVoxelization":
        """
        Create a voxelization that is the Boolean symmetric difference (XOR) of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean symmetric difference with.
        :type other_voxelization: PointVoxelization

        :return: The created voxelization resulting from the Boolean symmetric difference.
        :rtype: PointVoxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform symmetric difference.")

        return PointVoxelization(self.voxel_centers.symmetric_difference(other_voxelization.voxel_centers), self.voxel_size)

    def interference(self, other_voxelization: "PointVoxelization") -> float:
        """
        Compute the percentage of interference between two voxelization.

        :param other_voxelization: The other voxelization to compute percentage of interference with.
        :type other_voxelization: PointVoxelization

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
        :rtype: PointVoxelization
        """
        rotation_matrix = self._rotation_matrix(axis, angle)
        voxel_array = np.array(list(self.voxel_centers)) - np.array([center.x, center.y, center.z])
        rotated_voxels = np.dot(voxel_array, rotation_matrix.T)
        rotated_voxels += np.array([center.x, center.y, center.z])

        intersecting_voxels = self._voxels_intersecting_voxels(rotated_voxels, self.voxel_size)

        return PointVoxelization(intersecting_voxels, self.voxel_size)

    def translation(self, offset: Vector3D):
        """
        Translate the voxelization by the specified offset.

        :param offset: The translation offset.
        :type offset: Vector3D

        :return: A new Voxelization object resulting from the translation.
        :rtype: PointVoxelization
        """
        voxel_array = np.array(list(self.voxel_centers))
        translated_voxels = voxel_array + np.array([offset.x, offset.y, offset.z])

        intersecting_voxels = self._voxels_intersecting_voxels(translated_voxels, self.voxel_size)

        return PointVoxelization(intersecting_voxels, self.voxel_size)

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
    def from_voxel_matrix(cls, voxel_matrix: "VoxelMatrix"):
        """
        Create a Voxelization object from a voxel matrix.

        :param voxel_matrix: The voxel matrix object representing the voxelization.
        :type voxel_matrix: VoxelMatrix

        :return: A Voxelization object created from the voxel matrix.
        :rtype: PointVoxelization
        """
        if not cls.voxel_center_in_implicit_grid(voxel_matrix.matrix_origin_center, voxel_matrix.voxel_size):
            warnings.warn(
                """This voxel matrix is not defined in the implicit grid defined by the voxel_size. 
            Some methods like boolean operation or interference computing may not work as expected."""
            )

        indices = np.argwhere(voxel_matrix.matrix)
        voxel_centers = voxel_matrix.matrix_origin_center + indices * voxel_matrix.voxel_size

        return cls(set(map(tuple, np.round(voxel_centers, 6))), voxel_matrix.voxel_size)

    def _get_min_voxel_grid_center(self) -> Point:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.
        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, float, float]
        """
        min_x = min_y = min_z = float("inf")

        for point in self.voxel_centers:
            min_x = min(min_x, point[0])
            min_y = min(min_y, point[1])
            min_z = min(min_z, point[2])

        return min_x, min_y, min_z

    min_voxel_grid_center = property(_get_min_voxel_grid_center)

    def _get_max_voxel_grid_center(self) -> Point:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.
        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, float, float]
        """
        max_x = max_y = max_z = -float("inf")

        for point in self.voxel_centers:
            max_x = max(max_x, point[0])
            max_y = max(max_y, point[1])
            max_z = max(max_z, point[2])

        return max_x, max_y, max_z

    max_voxel_grid_center = property(_get_max_voxel_grid_center)

    def to_voxel_matrix(self) -> "VoxelMatrix":
        """
        Convert the voxelization to a voxel matrix object.

        :return: The voxel matrix representing the voxelization.
        :rtype: VoxelMatrix
        """
        min_center = self.min_voxel_grid_center
        max_center = self.max_voxel_grid_center

        dim_x = round((max_center[0] - min_center[0]) / self.voxel_size + 1)
        dim_y = round((max_center[1] - min_center[1]) / self.voxel_size + 1)
        dim_z = round((max_center[2] - min_center[2]) / self.voxel_size + 1)

        indices = np.round((np.array(list(self.voxel_centers)) - min_center) / self.voxel_size).astype(int)

        matrix = np.zeros((dim_x, dim_y, dim_z), dtype=np.bool_)
        matrix[indices[:, 0], indices[:, 1], indices[:, 2]] = True

        return VoxelMatrix(matrix, min_center, self.voxel_size)

    def inverse(self) -> "PointVoxelization":
        """
        Create a new Voxelization object that is the inverse of the current voxelization.

        :return: The inverse Voxelization object.
        :rtype: PointVoxelization
        """
        inverted_voxel_matrix = self.to_voxel_matrix().inverse()

        return PointVoxelization.from_voxel_matrix(inverted_voxel_matrix)

    def _get_bounding_box(self):
        """
        Get the bounding box of the voxelization.

        :return: The bounding box of the voxelization.
        :rtype: BoundingBox
        """
        min_point = (
            np.array([self.min_voxel_grid_center]) - np.array([self.voxel_size, self.voxel_size, self.voxel_size])
        )[0]
        max_point = (
            np.array([self.max_voxel_grid_center]) + np.array([self.voxel_size, self.voxel_size, self.voxel_size])
        )[0]

        return BoundingBox(min_point[0], max_point[0], min_point[1], max_point[1], min_point[2], max_point[2])

    bounding_box = property(_get_bounding_box)

    def _point_to_local_grid_index(self, point: Point) -> Tuple[int, ...]:
        """
        Convert a point to the local grid index within the voxelization.

        :param point: The point to convert.
        :type point: Point

        :return: The local grid index of the point.
        :rtype: Tuple[int, ...]

        :raises ValueError: If the point is not within the voxelization's bounding box.
        """
        if not self.bounding_box.point_belongs(Point3D(*point)):
            raise ValueError("Point not in local voxel grid.")

        x_index = int((point[0] - self.bounding_box.xmin) // self.voxel_size)
        y_index = int((point[1] - self.bounding_box.ymin) // self.voxel_size)
        z_index = int((point[2] - self.bounding_box.zmin) // self.voxel_size)

        return x_index, y_index, z_index

    def flood_fill(self, start_point: Point, fill_with: bool) -> "PointVoxelization":
        """
        Perform a flood fill operation within the voxelization starting from the given point.

        :param start_point: The starting point of the flood fill operation.
        :type start_point: Point
        :param fill_with: The value to fill the voxels during the flood fill operation.
        :type fill_with: bool

        :return: A new Voxelization object with the filled voxels.
        :rtype: PointVoxelization
        """
        start = self._point_to_local_grid_index(start_point)
        voxel_matrix = self.to_voxel_matrix()
        filled_voxel_matrix = voxel_matrix.flood_fill(start, fill_with)

        return self.from_voxel_matrix(filled_voxel_matrix)

    def fill_outer_voxels(self) -> "PointVoxelization":
        return self.from_voxel_matrix(self.to_voxel_matrix().fill_outer_voxels())

    def fill_enclosed_voxels(self) -> "PointVoxelization":
        return self.from_voxel_matrix(self.to_voxel_matrix().fill_enclosed_voxels())


class VoxelMatrix:
    """Class to manipulate voxel matrix."""

    def __init__(
        self,
        voxel_matrix: np.ndarray[np.bool_, np.ndim == 3],
        voxel_matrix_origin_center: Point,
        voxel_size: float,
    ):
        """
        :param voxel_matrix: The voxel numpy matrix object representing the voxelization.
        :type voxel_matrix: np.ndarray[np.bool_, np.ndim == 3]
        :param voxel_size: The size of the voxel edges.
        :type voxel_size: float
        :param voxel_matrix_origin_center: Voxel center of the origin of the voxel matrix, i.e 'matrix[0][0][0]'.
        :type voxel_matrix_origin_center: tuple[float, float, float]
        """
        self.matrix = voxel_matrix
        self.matrix_origin_center = voxel_matrix_origin_center
        self.voxel_size = voxel_size

    def __eq__(self, other_voxel_matrix: "VoxelMatrix") -> bool:
        return (
            self.voxel_size == other_voxel_matrix.voxel_size
            and self.matrix_origin_center == other_voxel_matrix.matrix_origin_center
            and np.array_equal(self.matrix, other_voxel_matrix.matrix)
        )

    def __len__(self) -> int:
        """
        Return the number of True value in the 3D voxel matrix.
        
        :return: The number of True value in the 3D voxel matrix.
        :rtype: int
        """
        return len(np.argwhere(self.matrix))

    def union(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        return self._logical_operation(other_voxel_matrix, np.logical_or)

    def __add__(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        """
        Return the union of the current voxel matrix with another voxel matrix.

        :param other_voxel_matrix: The voxel matrix to union with.
        :type other_voxel_matrix: VoxelMatrix

        :return: The union of the voxel matrices.
        :rtype: VoxelMatrix
        """
        return self.union(other_voxel_matrix)

    def difference(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        return self._logical_operation(other_voxel_matrix, lambda a, b: np.logical_and(a, np.logical_not(b)))

    def __sub__(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        """
        Return the difference between the current voxel matrix and another voxel matrix.

        :param other_voxel_matrix: The voxel matrix to subtract.
        :type other_voxel_matrix: VoxelMatrix

        :return: The difference between the voxel matrices.
        :rtype: VoxelMatrix
        """
        return self.difference(other_voxel_matrix)

    def intersection(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        return self._logical_operation(other_voxel_matrix, np.logical_and)

    def __and__(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        """
        Return the intersection of the current voxel matrix with another voxel matrix.

        :param other_voxel_matrix: The voxel matrix to intersect with.
        :type other_voxel_matrix: VoxelMatrix

        :return: The intersection of the voxel matrices.
        :rtype: VoxelMatrix
        """
        return self.intersection(other_voxel_matrix)

    def symmetric_difference(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        return self._logical_operation(other_voxel_matrix, np.logical_xor)

    def __xor__(self, other_voxel_matrix: "VoxelMatrix") -> "VoxelMatrix":
        """
        Return the symmetric difference between the current voxel matrix and another voxel matrix.

        :param other_voxel_matrix: The voxel matrix to calculate the symmetric difference with.
        :type other_voxel_matrix: VoxelMatrix
        :return: The symmetric difference between the voxel matrices.
        :rtype: VoxelMatrix
        """
        return self.symmetric_difference(other_voxel_matrix)

    def inverse(self) -> "VoxelMatrix":
        inverted_matrix = np.logical_not(self.matrix)
        return VoxelMatrix(inverted_matrix, self.matrix_origin_center, self.voxel_size)

    def __invert__(self) -> "VoxelMatrix":
        """
        Return the inverse of the current voxel matrix.

        :return: The inverse VoxelMatrix object.
        :rtype: VoxelMatrix
        """
        return self.inverse()

    def _get_extents(self):
        extents_min = np.round(np.array(self.matrix_origin_center), 6)
        extents_max = np.round(
            self.matrix_origin_center + self.matrix.shape * np.array([self.voxel_size for _ in range(3)]), 6
        )
        return extents_min, extents_max

    def _logical_operation(self, other: "VoxelMatrix", logical_operation):
        if self.voxel_size != other.voxel_size:
            raise ValueError("Voxel sizes must be the same to perform boolean operations.")

        self_min, self_max = self._get_extents()
        other_min, other_max = other._get_extents()

        global_min = np.min([self_min, other_min], axis=0)
        global_max = np.max([self_max, other_max], axis=0)

        new_shape = np.round((global_max - global_min) / self.voxel_size, 6).astype(int)

        new_self = np.zeros(new_shape, dtype=bool)
        new_other = np.zeros(new_shape, dtype=bool)

        self_start = np.round((self_min - global_min) / self.voxel_size, 6).astype(int)
        other_start = np.round((other_min - global_min) / self.voxel_size, 6).astype(int)

        new_self[
            self_start[0] : self_start[0] + self.matrix.shape[0],
            self_start[1] : self_start[1] + self.matrix.shape[1],
            self_start[2] : self_start[2] + self.matrix.shape[2],
        ] = self.matrix

        new_other[
            other_start[0] : other_start[0] + other.matrix.shape[0],
            other_start[1] : other_start[1] + other.matrix.shape[1],
            other_start[2] : other_start[2] + other.matrix.shape[2],
        ] = other.matrix

        result_matrix = logical_operation(new_self, new_other)

        return VoxelMatrix(result_matrix, tuple(global_min), self.voxel_size)._crop_matrix()

    def _crop_matrix(self) -> "VoxelMatrix":
        """
        Crop the VoxelMatrix to the smallest possible size.
        """
        # Find the indices of the True voxels
        true_voxels = np.argwhere(self.matrix)

        if len(true_voxels) == 0:
            # Can't crop a matrix of False values
            return self

        # Find the minimum and maximum indices along each axis
        min_voxel_coords, max_voxel_coords = np.round(np.min(true_voxels, axis=0), 6), np.round(
            np.max(true_voxels, axis=0), 6
        )

        # Crop the matrix to the smallest possible size
        cropped_matrix = self.matrix[
            min_voxel_coords[0] : max_voxel_coords[0] + 1,
            min_voxel_coords[1] : max_voxel_coords[1] + 1,
            min_voxel_coords[2] : max_voxel_coords[2] + 1,
        ]

        # Calculate new matrix_origin_center
        new_origin_center = np.round(self.matrix_origin_center + min_voxel_coords * self.voxel_size, 6)

        return self.__class__(cropped_matrix, tuple(new_origin_center), self.voxel_size)

    def flood_fill(self, start, fill_with) -> "VoxelMatrix":
        return VoxelMatrix(
            flood_fill_matrix_3d(self.matrix, start, fill_with), self.matrix_origin_center, self.voxel_size
        )

    def fill_outer_voxels(self) -> "VoxelMatrix":
        expanded_voxel_matrix = self._expand()
        outer_filled_expanded_voxel_matrix = expanded_voxel_matrix.flood_fill((0, 0, 0), True)
        outer_filled_voxel_matrix = outer_filled_expanded_voxel_matrix._reduce()

        return outer_filled_voxel_matrix

    def fill_enclosed_voxels(self) -> "VoxelMatrix":
        outer_filled_voxel_matrix = self.fill_outer_voxels()
        inner_filled_voxel_matrix = self + outer_filled_voxel_matrix.inverse()

        # if inner_filled_voxel_matrix == self:
        #     warnings.warn("This voxelization doesn't have any enclosed voxels.")

        return inner_filled_voxel_matrix

    def _expand(self) -> "VoxelMatrix":
        current_shape = self.matrix.shape
        new_shape = tuple(dim + 2 for dim in current_shape)
        expanded_matrix = np.zeros(new_shape, dtype="bool")
        slices = tuple(slice(1, -1) for _ in current_shape)
        expanded_matrix[slices] = self.matrix.copy()

        return VoxelMatrix(
            expanded_matrix,
            tuple(round(coord - self.voxel_size, 6) for coord in self.matrix_origin_center),
            self.voxel_size,
        )

    def _reduce(self) -> "VoxelMatrix":
        current_shape = self.matrix.shape
        slices = tuple(slice(1, -1) for _ in current_shape)
        reduced_matrix = self.matrix.copy()[slices]

        return VoxelMatrix(
            reduced_matrix,
            tuple(round(coord + self.voxel_size, 6) for coord in self.matrix_origin_center),
            self.voxel_size,
        )

    @classmethod
    def from_voxelization(cls, voxelization: "PointVoxelization") -> "VoxelMatrix":
        return voxelization.to_voxel_matrix()

    @classmethod
    def from_shell(cls, shell: Shell3D, voxel_size: float, name: str = ""):
        triangles = _shell_to_triangles(shell)
        matrix, matrix_origin_center = triangles_to_voxel_matrix(triangles, voxel_size)

        return cls(matrix, matrix_origin_center, voxel_size)._crop_matrix()

    @classmethod
    def from_volume_model(cls, volume_model: VolumeModel, voxel_size: float, name: str = ""):
        triangles = _volume_model_to_triangles(volume_model)
        matrix, matrix_origin_center = triangles_to_voxel_matrix(triangles, voxel_size)

        return cls(matrix, matrix_origin_center, voxel_size)._crop_matrix()


    def to_voxelization(self) -> "PointVoxelization":
        return PointVoxelization.from_voxel_matrix(self)

    def to_triangles(self) -> Set[Triangle]:
        return self.to_voxelization().to_triangles()


class Pixelization:
    """Voxelization but in 2D"""

    def __init__(self, pixel_centers, pixel_size):
        self.pixel_centers = pixel_centers
        self.pixel_size = pixel_size

    def __eq__(self, other_pixelization: "Pixelization") -> bool:
        """
        Check if the current pixelization is equal to another pixelization.

        :param other_pixelization: The pixelization to compare.
        :type other_pixelization: Pixelization

        :return: True if the pixelizations are equal, False otherwise.
        :rtype: bool
        """
        return (
            self.pixel_centers == other_pixelization.pixel_centers and self.pixel_size == other_pixelization.pixel_size
        )

    def __add__(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Return the union of the current pixelization with another pixelization.

        :param other_pixelization: The pixelization to union with.
        :type other_pixelization: Pixelization

        :return: The union of the pixelizations.
        :rtype: Pixelization
        """
        return self.union(other_pixelization)

    def __sub__(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Return the difference between the current pixelization and another pixelization.

        :param other_pixelization: The pixelization to subtract.
        :type other_pixelization: Pixelization

        :return: The difference between the pixelizations.
        :rtype: Pixelization
        """
        return self.difference(other_pixelization)

    def __and__(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Return the intersection of the current pixelization with another pixelization.

        :param other_pixelization: The pixelization to intersect with.
        :type other_pixelization: Pixelization

        :return: The intersection of the pixelizations.
        :rtype: Pixelization
        """
        return self.intersection(other_pixelization)

    def __or__(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Return the union of the current pixelization with another pixelization.

        :param other_pixelization: The pixelization to union with.
        :type other_pixelization: Pixelization

        :return: The union of the pixelizations.
        :rtype: Pixelization
        """
        return self.union(other_pixelization)

    def __xor__(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Return the symmetric difference between the current pixelization and another pixelization.

        :param other_pixelization: The pixelization to calculate the symmetric difference with.
        :type other_pixelization: Pixelization
        :return: The symmetric difference between the pixelizations.
        :rtype: Pixelization
        """
        return self.symmetric_difference(other_pixelization)

    def __invert__(self) -> "Pixelization":
        """
        Return the inverse of the current pixelization.

        :return: The inverse Pixelization object.
        :rtype: Pixelization
        """
        return self.inverse()

    def __len__(self):
        """
        Return the number of pixels in the pixelization.

        :return: The number of pixels.
        :rtype: int
        """
        return len(self.pixel_centers)

    def intersection(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Create a pixelization that is the Boolean intersection of two pixelizations.
        Both pixelizations must have the same pixel size.

        :param other_pixelization: The other pixelization to compute the Boolean intersection with.
        :type other_pixelization: Pixelization

        :return: The created pixelization resulting from the Boolean intersection.
        :rtype: Pixelization
        """
        if self.pixel_size != other_pixelization.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform intersection.")

        return Pixelization(self.pixel_centers.intersection(other_pixelization.pixel_centers), self.pixel_size)

    def is_intersecting(self, other_pixelization: "Pixelization") -> bool:
        """
        Check if two pixelizations are intersecting.
        Both pixelizations must have the same pixel size.

        :param other_pixelization: The other pixelization to check if there is an intersection with.
        :type other_pixelization: Pixelization

        :return: True if the pixelizations are intersecting, False otherwise.
        :rtype: bool
        """
        intersection = self.intersection(other_pixelization)

        return len(intersection) > 0

    def union(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Create a pixelization that is the Boolean union of two pixelizations.
        Both pixelizations must have the same pixel size.

        :param other_pixelization: The other pixelization to compute the Boolean union with.
        :type other_pixelization: Pixelization

        :return: The created pixelization resulting from the Boolean union.
        :rtype: Pixelization
        """
        if self.pixel_size != other_pixelization.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform union.")

        return Pixelization(self.pixel_centers.union(other_pixelization.pixel_centers), self.pixel_size)

    def difference(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Create a pixelization that is the Boolean difference of two pixelizations.
        Both pixelizations must have the same pixel size.

        :param other_pixelization: The other pixelization to compute the Boolean difference with.
        :type other_pixelization: Pixelization

        :return: The created pixelization resulting from the Boolean difference.
        :rtype: Pixelization
        """
        if self.pixel_size != other_pixelization.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform difference.")

        return Pixelization(self.pixel_centers.difference(other_pixelization.pixel_centers), self.pixel_size)

    def symmetric_difference(self, other_pixelization: "Pixelization") -> "Pixelization":
        """
        Create a pixelization that is the Boolean symmetric difference (XOR) of two pixelizations.
        Both pixelizations must have the same pixel size.

        :param other_pixelization: The other pixelization to compute the Boolean symmetric difference with.
        :type other_pixelization: Pixelization

        :return: The created pixelization resulting from the Boolean symmetric difference.
        :rtype: Pixelization
        """
        if self.pixel_size != other_pixelization.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform symmetric difference.")

        return Pixelization(self.pixel_centers.symmetric_difference(other_pixelization.pixel_centers), self.pixel_size)

    def interference(self, other_pixelization: "Pixelization") -> float:
        """
        Compute the percentage of interference between two pixelizations.

        :param other_pixelization: The other pixelization to compute the percentage of interference with.
        :type other_pixelization: Pixelization

        :return: The percentage of interference between the two pixelizations.
        :rtype: float
        """
        return len(self.intersection(other_pixelization)) / len(self.union(other_pixelization))

    def plot(self, ax=None):
        """Plots the pixels on a 2D plane"""
        if ax is None:
            _, ax = plt.subplots()

        for center in self.pixel_centers:
            x, y = center
            ax.add_patch(
                patches.Rectangle(
                    (x - self.pixel_size / 2, y - self.pixel_size / 2),  # Bottom left corner
                    self.pixel_size,  # Width
                    self.pixel_size,  # Height
                    color="black",
                )
            )

        # Setting the x and y limits to contain all the pixels
        ax.set_xlim(
            min(x[0] for x in self.pixel_centers) - self.pixel_size,
            max(x[0] for x in self.pixel_centers) + self.pixel_size,
        )
        ax.set_ylim(
            min(x[1] for x in self.pixel_centers) - self.pixel_size,
            max(x[1] for x in self.pixel_centers) + self.pixel_size,
        )

        ax.set_aspect("equal")  # Ensuring equal scaling for both axes
        plt.show()

        return ax

    @classmethod
    def from_line_segment(cls, line_segment, pixel_size):
        return cls(cls._line_segments_to_pixels([line_segment], pixel_size), pixel_size)

    @classmethod
    def from_polygon(cls, polygon, pixel_size):
        line_segments = [(polygon[i - 1], polygon[i]) for i in range(len(polygon))]
        return cls(cls._line_segments_to_pixels(line_segments, pixel_size), pixel_size)

    @staticmethod
    def _line_segments_to_pixels(line_segments, pixel_size):
        return line_segments_to_pixels(line_segments, pixel_size)

    @classmethod
    def from_pixel_matrix(
        cls, pixel_matrix: "PixelMatrix", pixel_size: float, pixel_matrix_origin_center: Tuple[float, float]
    ):
        indices = np.argwhere(pixel_matrix.matrix)
        pixel_centers = pixel_matrix_origin_center + indices * pixel_size
        return cls(set(map(tuple, np.round(pixel_centers, 6))), pixel_size)

    def _get_min_pixel_grid_center(self) -> Tuple[float, float]:
        min_x = min_y = float("inf")
        for point in self.pixel_centers:
            min_x = min(min_x, point[0])
            min_y = min(min_y, point[1])
        return min_x, min_y

    min_pixel_grid_center = property(_get_min_pixel_grid_center)

    def _get_max_pixel_grid_center(self) -> Tuple[float, float]:
        max_x = max_y = -float("inf")
        for point in self.pixel_centers:
            max_x = max(max_x, point[0])
            max_y = max(max_y, point[1])
        return max_x, max_y

    max_pixel_grid_center = property(_get_max_pixel_grid_center)

    def to_pixel_matrix(self) -> "PixelMatrix":
        min_center = self.min_pixel_grid_center
        max_center = self.max_pixel_grid_center

        dim_x = round((max_center[0] - min_center[0]) / self.pixel_size + 1)
        dim_y = round((max_center[1] - min_center[1]) / self.pixel_size + 1)

        indices = np.round((np.array(list(self.pixel_centers)) - min_center) / self.pixel_size).astype(int)

        matrix = np.zeros((dim_x, dim_y), dtype=np.bool_)
        matrix[indices[:, 0], indices[:, 1]] = True

        return PixelMatrix(matrix)

    def inverse(self) -> "Pixelization":
        """
        Create a new Pixelization object that is the inverse of the current pixelization.

        :return: The inverse Pixelization object.
        :rtype: Pixelization
        """
        inverted_pixel_matrix = self.to_pixel_matrix().inverse()
        min_pixel_center = self.min_pixel_grid_center

        return Pixelization.from_pixel_matrix(inverted_pixel_matrix, self.pixel_size, min_pixel_center)

    def _get_bounding_rectangle(self):
        min_point = np.array([self.min_pixel_grid_center]) - np.array([self.pixel_size, self.pixel_size])
        max_point = np.array([self.max_pixel_grid_center]) + np.array([self.pixel_size, self.pixel_size])

        return BoundingRectangle(min_point[0], max_point[0], min_point[1], max_point[1])

    bounding_rectangle = property(_get_bounding_rectangle)

    def _point_to_local_grid_index(self, point: Tuple[float, float]) -> Tuple[int, int]:
        if not self.bounding_rectangle.point_belongs(Point2D(*point)):
            raise ValueError("Point not in local pixel grid.")

        x_index = int((point[0] - self.bounding_rectangle.xmin) // self.pixel_size)
        y_index = int((point[1] - self.bounding_rectangle.ymin) // self.pixel_size)

        return x_index, y_index

    def flood_fill(self, start_point: Tuple[float, float], fill_with: bool) -> "Pixelization":
        start = self._point_to_local_grid_index(start_point)
        pixel_matrix = self.to_pixel_matrix()
        filled_pixel_matrix = pixel_matrix.flood_fill(start, fill_with)

        return self.from_pixel_matrix(filled_pixel_matrix, self.pixel_size, self.min_pixel_grid_center)

    def fill_outer_pixels(self) -> "Pixelization":
        return self.from_pixel_matrix(
            self.to_pixel_matrix().fill_outer_pixels(), self.pixel_size, self.min_pixel_grid_center
        )

    def fill_enclosed_pixels(self) -> "Pixelization":
        return self.from_pixel_matrix(
            self.to_pixel_matrix().fill_enclosed_pixels(), self.pixel_size, self.min_pixel_grid_center
        )


class PixelMatrix:
    """Class to manipulate pixel matrix."""

    def __init__(self, numpy_pixel_matrix: np.ndarray):  # np.ndarray[np.bool_, np.ndim == 2]
        self.matrix = numpy_pixel_matrix

    def __eq__(self, other_pixel_matrix: "PixelMatrix") -> bool:
        return np.array_equal(self.matrix, other_pixel_matrix.matrix)

    def __add__(self, other_pixel_matrix: "PixelMatrix") -> "PixelMatrix":
        return PixelMatrix(self.matrix + other_pixel_matrix.matrix)

    def inverse(self) -> "PixelMatrix":
        inverted_matrix = np.logical_not(self.matrix)
        return PixelMatrix(inverted_matrix)

    def flood_fill(self, start: Tuple[int, int], fill_with: bool) -> "PixelMatrix":
        return PixelMatrix(flood_fill_matrix_2d(self.matrix, start, fill_with))

    def _expand(self) -> "PixelMatrix":
        current_shape = self.matrix.shape
        new_shape = tuple(dim + 2 for dim in current_shape)
        expanded_matrix = np.zeros(new_shape, dtype="bool")
        slices = tuple(slice(1, -1) for _ in current_shape)
        expanded_matrix[slices] = self.matrix.copy()

        return PixelMatrix(expanded_matrix)

    def _reduce(self) -> "PixelMatrix":
        current_shape = self.matrix.shape
        slices = tuple(slice(1, -1) for _ in current_shape)
        reduced_matrix = self.matrix.copy()[slices]

        return PixelMatrix(reduced_matrix)

    def fill_outer_pixels(self) -> "PixelMatrix":
        expanded_pixel_matrix = self._expand()
        outer_filled_expanded_pixel_matrix = expanded_pixel_matrix.flood_fill((0, 0), True)
        outer_filled_pixel_matrix = outer_filled_expanded_pixel_matrix._reduce()

        return outer_filled_pixel_matrix

    def fill_enclosed_pixels(self) -> "PixelMatrix":
        outer_filled_pixel_matrix = self.fill_outer_pixels()
        inner_filled_pixel_matrix = self + outer_filled_pixel_matrix.inverse()

        # if inner_filled_pixel_matrix == self:
        #     warnings.warn("This pixel matrix doesn't have any enclosed pixels.")

        return inner_filled_pixel_matrix


def _shell_to_triangles(shell: Shell3D) -> List[Triangle]:
    """
    Helper function to convert a Shell3D to a list of triangles.
    It uses the "triangulation" method to triangulate the Shell3D.

    :param shell: The Shell3D to convert to triangles.
    :type shell: Shell3D

    :return: The list of triangles extracted from the triangulated Shell3D.
    :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
    """
    triangulation = shell.triangulation()
    return [
        tuple(
            tuple([float(point.x), float(point.y), float(point.z)])
            for point in [
                triangulation.points[triangle[0]],
                triangulation.points[triangle[1]],
                triangulation.points[triangle[2]],
            ]
        )
        for triangle in triangulation.triangles
    ]


def _volume_model_to_triangles(volume_model: VolumeModel) -> List[Triangle]:
    """
    Helper function to convert a VolumeModel to a list of triangles.
    It uses the "triangulation" method to triangulate the shells of the VolumeModel.

    :param volume_model: The VolumeModel to convert to triangles.
    :type volume_model: VolumeModel

    :return: The list of triangles extracted from the triangulated primitives of the VolumeModel.
    :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
    """
    triangles = []
    for shell in volume_model.get_shells():
        triangles.extend(_shell_to_triangles(shell))

    return triangles
