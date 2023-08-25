"""
Class for voxel representation of volmdlr models
"""
# pylint: disable=no-name-in-module

import warnings
from abc import ABC, abstractmethod
from typing import List, Set, Tuple, TypeVar

from matplotlib import patches
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from dessia_common.core import PhysicalObject

from volmdlr import Point2D, Point3D, Vector3D
from volmdlr.core import BoundingBox, BoundingRectangle, VolumeModel
from volmdlr.faces import Triangle3D
from volmdlr.shells import ClosedTriangleShell3D, Shell3D
from volmdlr.voxelization_compiled import (
    flood_fill_matrix_2d,
    flood_fill_matrix_3d,
    line_segments_to_pixels,
    triangles_to_voxel_matrix,
    voxel_triangular_faces,
)

# CUSTOM TYPES
Point = Tuple[float, float, float]
Triangle = Tuple[Point, Point, Point]
Segment = Tuple[Point, Point]

# GLOBAL VARIABLE
DECIMALS = 9  # Used to round numbers and avoid floating point arithmetic imprecision


class Voxelization(ABC, PhysicalObject):
    """
    Abstract base class for creating and manipulating voxelizations of volmdlr geometries.

    This approach is used to create a voxelization of the surfaces, without filling the volume.

    The voxelization is defined on an implicit 3D grid, where each voxel is defined by its size.
    The implicit grid consists of axis-aligned cubes with a given size, ensuring that the global
    origin (0, 0, 0) is always a corner point of a voxel.

    For example, in 1D with a voxel size of `t`, the set of voxels (defined by minimum and maximum points)
    is: {i ∈ N, t ∈ R | (i * t, (i+1) * t)}
    The corresponding set of voxel centers is: {i ∈ N, t ∈ R | (i + 0.5) * t}

    This approach enables consistent voxelization across the 3D space, facilitating fast Boolean operations.
    """

    VoxelizationType = TypeVar("VoxelizationType", bound="Voxelization")

    def __init__(self, voxel_size: float, name: str):
        self.voxel_size = voxel_size

        PhysicalObject.__init__(self, name=name)

    @property
    @abstractmethod
    def voxel_centers(self) -> Set[Point]:
        """
        Get the center point of each voxel.

        :return: The center point of each voxel.
        :rtype: set[tuple[float, float, float]]
        """

    def __str__(self):
        """
        Return a custom string representation of the voxelization.

        :return: A string representation of the voxelization.
        :rtype: str
        """
        return f"Voxelization: voxel size={self.voxel_size}, number of voxels={self.__len__()}, name={self.name}"

    @abstractmethod
    def __eq__(self, other_voxelization: VoxelizationType) -> bool:
        """
        Check if two voxelizations are equal.

        :param other_voxelization: Another voxelization to compare with.
        :type other_voxelization: VoxelizationType

        :return: True if the voxelizations are equal, False otherwise.
        :rtype: bool
        """

    @abstractmethod
    def __len__(self) -> int:
        """
        Get the number of voxels in the voxelization.

        :return: The number of voxels in the voxelization.
        :rtype: int
        """

    @property
    def volume(self) -> float:
        """
        Calculate the volume of the voxelization.

        :return: The volume of the voxelization.
        :rtype: float
        """
        return len(self) * self.voxel_size**3

    @property
    @abstractmethod
    def min_voxel_grid_center(self) -> Point:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, float, float]
        """

    @property
    @abstractmethod
    def max_voxel_grid_center(self) -> Point:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, float, float]
        """

    @property
    def bounding_box(self):
        """
        Get the bounding box of the voxelization.

        :return: The bounding box of the voxelization.
        :rtype: BoundingBox
        """
        min_point = np.round((np.array([self.min_voxel_grid_center]) - self.voxel_size)[0], DECIMALS)
        max_point = np.round((np.array([self.max_voxel_grid_center]) + self.voxel_size)[0], DECIMALS)

        return BoundingBox(min_point[0], max_point[0], min_point[1], max_point[1], min_point[2], max_point[2])

    # CLASS METHODS
    @classmethod
    @abstractmethod
    def from_shell(cls, shell: Shell3D, voxel_size: float, name: str = "") -> VoxelizationType:
        """
        Create a voxelization from a Shell3D.

        :param shell: The Shell3D to create the voxelization from.
        :type shell: Shell3D
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the Shell3D.
        :rtype: VoxelizationType
        """

    @classmethod
    @abstractmethod
    def from_volume_model(cls, volume_model: VolumeModel, voxel_size: float, name: str = "") -> VoxelizationType:
        """
        Create a voxelization from a VolumeModel.

        :param volume_model: The VolumeModel to create the voxelization from.
        :type volume_model: VolumeModel
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the VolumeModel.
        :rtype: VoxelizationType
        """

    # BOOLEAN OPERATIONS
    def __add__(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Overloaded '+' operator for performing a union operation.

        :param other_voxelization: The voxelization to perform the union with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the union operation.
        :rtype: VoxelizationType
        """
        return self.union(other_voxelization)

    def __sub__(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Overloaded '-' operator for performing a difference operation.

        :param other_voxelization: The voxelization to perform the difference with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the difference operation.
        :rtype: VoxelizationType
        """
        return self.difference(other_voxelization)

    def __and__(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Overloaded '&' operator for performing an intersection operation.

        :param other_voxelization: The voxelization to perform the intersection with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the intersection operation.
        :rtype: VoxelizationType
        """
        return self.intersection(other_voxelization)

    def __xor__(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Overloaded '^' operator for performing a symmetric difference operation.

        :param other_voxelization: The voxelization to perform the symmetric difference with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the symmetric difference operation.
        :rtype: VoxelizationType
        """
        return self.symmetric_difference(other_voxelization)

    def __invert__(self) -> VoxelizationType:
        """
        Overloaded '~' operator for computing the inverse.

        :return: A new voxelization representing the inverse.
        :rtype: Voxelization
        """
        return self.inverse()

    @abstractmethod
    def union(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Perform a union operation with another voxelization.

        :param other_voxelization: The voxelization to perform the union with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the union operation.
        :rtype: VoxelizationType
        """

    @abstractmethod
    def difference(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Perform a difference operation with another voxelization.

        :param other_voxelization: The voxelization to perform the difference with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the difference operation.
        :rtype: VoxelizationType
        """

    @abstractmethod
    def intersection(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Perform an intersection operation with another voxelization.

        :param other_voxelization: The voxelization to perform the intersection with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the intersection operation.
        :rtype: VoxelizationType
        """

    def symmetric_difference(self, other_voxelization: VoxelizationType) -> VoxelizationType:
        """
        Perform a symmetric difference operation with another voxelization.

        :param other_voxelization: The voxelization to perform the symmetric difference with.
        :type other_voxelization: VoxelizationType

        :return: A new voxelization resulting from the symmetric difference operation.
        :rtype: VoxelizationType
        """

    @abstractmethod
    def inverse(self) -> VoxelizationType:
        """
        Compute the inverse of the voxelization.

        :return: A new voxelization representing the inverse.
        :rtype: VoxelizationType
        """

    def interference(self, other_voxelization: VoxelizationType) -> float:
        """
        Compute the percentage of interference between two voxelizations.

        :param other_voxelization: The other voxelization to compute interference with.
        :type other_voxelization: Voxelization

        :return: The percentage of interference between the two voxelizations.
        :rtype: float
        """
        return len(self.intersection(other_voxelization)) / len(self.union(other_voxelization))

    def is_intersecting(self, other_voxelization: VoxelizationType) -> bool:
        """
        Check if two voxelizations are intersecting.

        :param other_voxelization: The other voxelization to check if there is an intersection with.
        :type other_voxelization: VoxelizationType

        :return: True if the voxelizations are intersecting, False otherwise.
        :rtype: bool
        """
        intersection = self.intersection(other_voxelization)

        return len(intersection) > 0

    # FILLING METHODS
    @abstractmethod
    def flood_fill(self, start, fill_with: bool) -> VoxelizationType:
        """
        Perform a flood fill operation on the voxelization.

        :param start: The starting point for the flood fill.
        :param fill_with: The value to fill the voxels with during the operation.
        :type fill_with: bool

        :return: A new voxelization resulting from the flood fill operation.
        :rtype: VoxelizationType
        """

    @abstractmethod
    def fill_outer_voxels(self) -> VoxelizationType:
        """
        Fill the outer voxels of the voxelization.

        :return: A new voxelization with outer voxels filled.
        :rtype: VoxelizationType
        """

    @abstractmethod
    def fill_enclosed_voxels(self) -> VoxelizationType:
        """
        Fill the enclosed voxels of the voxelization.

        :return: A new voxelization with enclosed voxels filled.
        :rtype: VoxelizationType
        """

    # DISPLAY METHODS
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
            for triangle in voxel_triangular_faces(voxel, self.voxel_size):
                if triangle not in triangles:
                    triangles.add(triangle)
                else:
                    triangles.remove(triangle)

        return triangles

    def to_closed_triangle_shell(self) -> ClosedTriangleShell3D:
        """
        Generate a closed triangle shell representing the voxelization.

        :return: A closed triangle shell representation of the voxelization.
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
        Generate volmdlr primitives.

        :param kwargs: Additional keyword arguments.

        :return: A list of volmdlr primitives.
        :rtype: List[ClosedTriangleShell3D]
        """
        return [self.to_closed_triangle_shell()]

    # HELPER METHODS
    @staticmethod
    def check_voxel_center_is_in_implicit_grid(voxel_center: Point, voxel_size: float) -> bool:
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
            if not round((coord - 0.5 * voxel_size) / voxel_size, DECIMALS).is_integer():
                return False

        return True

    def is_hollow(self) -> bool:
        """
        Check if the voxelization is hollow.

        A hollow voxelization is one that has enclosed voxels that are not filled.

        :return: True if the voxelization is hollow, False otherwise.
        :rtype: bool
        """
        return self.fill_enclosed_voxels() == self

    @staticmethod
    def _check_voxel_size_number_of_decimals(voxel_size: float):
        """
        Check the number of decimal places in the voxel size.

        If the voxel size has more decimal places than a specified threshold (DECIMALS),
        a warning is issued, as some functions may not work as intended.

        This is due to the use of rounding functions to avoid floating point arithmetic imprecision.

        :param voxel_size: The size of the voxels.
        :type voxel_size: float

        :raises ValueError: If voxel_size is not a float.
        """
        if isinstance(voxel_size, float):
            decimal_part = abs(voxel_size - int(voxel_size))
            if decimal_part == 0:
                pass
            else:
                decimals = len(str(decimal_part).split(".")[1])
                if decimals >= DECIMALS:
                    warnings.warn(
                        f"""voxel_size has too many decimals: some functions may not work as intended.
                        Consider using a voxel_size with less than {DECIMALS}."""
                    )
        else:
            raise ValueError("voxel_size is not a float")

    def _check_other_voxelization_type(self, other_voxelization):
        """
        Check if the provided 'other_voxelization' is an instance of the same Voxelization subclass.

        :param other_voxelization: Another voxelization to be checked.
        :type other_voxelization: Voxelization

        :raises ValueError: If 'other_voxelization' is not an instance of the same Voxelization subclass.
        """
        if not isinstance(other_voxelization, self.__class__):
            raise ValueError(f"'other_voxelization' must be an instance of '{self.__class__.__name__}'")

    def _check_other_voxelization_voxel_size(self, other_voxelization: VoxelizationType):
        """
        Check if the provided 'other_voxelization' has the same voxel size.

        :param other_voxelization: Another voxelization to be checked.
        :type other_voxelization: Voxelization

        :raises ValueError: If 'other_voxelization' has not the same voxel size.
        """
        if not self.voxel_size == other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform this operation.")

    @staticmethod
    def _shell_to_triangles(shell: Shell3D) -> List[Triangle]:
        """
        Helper method to convert a Shell3D to a list of triangles.

        It uses the "triangulation" method to triangulate the Shell3D.

        :param shell: The Shell3D to convert to triangles.
        :type shell: Shell3D

        :return: The list of triangles extracted from the triangulated Shell3D.
        :rtype: List[Triangle]
        """
        triangulation = shell.triangulation()
        return [
            (
                (
                    float(triangulation.points[triangle[0]].x),
                    float(triangulation.points[triangle[0]].y),
                    float(triangulation.points[triangle[0]].z),
                ),
                (
                    float(triangulation.points[triangle[1]].x),
                    float(triangulation.points[triangle[1]].y),
                    float(triangulation.points[triangle[1]].z),
                ),
                (
                    float(triangulation.points[triangle[2]].x),
                    float(triangulation.points[triangle[2]].y),
                    float(triangulation.points[triangle[2]].z),
                ),
            )
            for triangle in triangulation.triangles
        ]

    @staticmethod
    def _volume_model_to_triangles(volume_model: VolumeModel) -> List[Triangle]:
        """
        Helper method to convert a VolumeModel to a list of triangles.

        It uses the "triangulation" method to triangulate the shells of the VolumeModel.

        :param volume_model: The VolumeModel to convert to triangles.
        :type volume_model: VolumeModel

        :return: The list of triangles extracted from the triangulated primitives of the VolumeModel.
        :rtype: List[Triangle]
        """
        triangles = []
        for shell in volume_model.get_shells():
            triangles.extend(Voxelization._shell_to_triangles(shell))

        return triangles


class PointBasedVoxelization(Voxelization):
    """Voxelization implemented as a set of points, representing each voxel center."""

    def __init__(self, voxel_centers: Set[Point], voxel_size: float, name: str = ""):
        """
        Initialize the PointBasedVoxelization.

        :param voxel_centers: The set of points representing voxel centers.
        :type voxel_centers: set[tuple[float, float, float]]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the Voxelization.
        :type name: str, optional
        """
        self._check_voxel_size_number_of_decimals(voxel_size)

        self._voxel_centers = voxel_centers

        Voxelization.__init__(self, voxel_size=voxel_size, name=name)

    @property
    def voxel_centers(self) -> Set[Point]:
        """
        Get the center point of each voxel.

        :return: The center point of each voxel.
        :rtype: set[tuple[float, float, float]]
        """
        return self._voxel_centers

    def __eq__(self, other_voxelization: "PointBasedVoxelization") -> bool:
        """
        Check if two voxelizations are equal.

        :param other_voxelization: Another voxelization to compare with.
        :type other_voxelization: PointBasedVoxelization

        :return: True if the voxelizations are equal, False otherwise.
        :rtype: bool
        """
        self._check_other_voxelization_type(other_voxelization)

        return (
            self.voxel_centers == other_voxelization.voxel_centers and self.voxel_size == other_voxelization.voxel_size
        )

    def __len__(self) -> int:
        """
        Get the number of voxels in the voxelization.

        :return: The number of voxels in the voxelization (i.e. the number of voxel centers).
        :rtype: int
        """
        return len(self.voxel_centers)

    @property
    def min_voxel_grid_center(self) -> Point:
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

    @property
    def max_voxel_grid_center(self) -> Point:
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

    # CLASS METHODS
    @classmethod
    def from_shell(cls, shell: Shell3D, voxel_size: float, name: str = "") -> "PointBasedVoxelization":
        """
        Create a PointBasedVoxelization from a Shell3D.

        :param shell: The Shell3D to create the voxelization from.
        :type shell: Shell3D
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the Shell3D.
        :rtype: PointBasedVoxelization
        """
        voxels = MatrixBasedVoxelization.from_shell(shell, voxel_size).voxel_centers

        return cls(voxel_centers=voxels, voxel_size=voxel_size, name=name)

    @classmethod
    def from_volume_model(
        cls, volume_model: VolumeModel, voxel_size: float, name: str = ""
    ) -> "PointBasedVoxelization":
        """
        Create a PointBasedVoxelization from a VolumeModel.

        :param volume_model: The VolumeModel to create the voxelization from.
        :type volume_model: VolumeModel
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the VolumeModel.
        :rtype: PointBasedVoxelization
        """
        voxels = MatrixBasedVoxelization.from_volume_model(volume_model, voxel_size).voxel_centers

        return cls(voxel_centers=voxels, voxel_size=voxel_size, name=name)

    @classmethod
    def from_matrix_based_voxelization(
        cls, matrix_based_voxelization: "MatrixBasedVoxelization"
    ) -> "PointBasedVoxelization":
        """
        Create a PointBasedVoxelization object from a MatrixBasedVoxelization.

        :param matrix_based_voxelization: The MatrixBasedVoxelization object representing the voxelization.
        :type matrix_based_voxelization: MatrixBasedVoxelization

        :return: A PointBasedVoxelization object created from the MatrixBasedVoxelization.
        :rtype: PointBasedVoxelization
        """
        if not cls.check_voxel_center_is_in_implicit_grid(
            matrix_based_voxelization.min_voxel_grid_center, matrix_based_voxelization.voxel_size
        ):
            warnings.warn(
                """This matrix based voxelization is not defined in the implicit grid defined by the voxel_size. 
            Some methods like boolean operation or interference computing may not work as expected."""
            )

        return matrix_based_voxelization.to_point_voxelization()

    # BOOLEAN OPERATIONS
    def union(self, other_voxelization: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Perform a union operation with another PointBasedVoxelization.

        :param other_voxelization: The PointBasedVoxelization to perform the union with.
        :type other_voxelization: PointBasedVoxelization

        :return: A new PointBasedVoxelization resulting from the union operation.
        :rtype: PointBasedVoxelization
        """
        self._check_other_voxelization_type(other_voxelization)
        self._check_other_voxelization_voxel_size(other_voxelization)

        return self.__class__(self.voxel_centers.union(other_voxelization.voxel_centers), self.voxel_size)

    def difference(self, other_voxelization: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Perform an intersection operation with another PointBasedVoxelization.

        :param other_voxelization: The PointBasedVoxelization to perform the intersection with.
        :type other_voxelization: PointBasedVoxelization

        :return: A new PointBasedVoxelization resulting from the intersection operation.
        :rtype: PointBasedVoxelization
        """
        self._check_other_voxelization_type(other_voxelization)
        self._check_other_voxelization_voxel_size(other_voxelization)

        return self.__class__(self.voxel_centers.difference(other_voxelization.voxel_centers), self.voxel_size)

    def intersection(self, other_voxelization: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Create a voxelization that is the Boolean intersection of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the Boolean intersection with.
        :type other_voxelization: PointBasedVoxelization

        :return: The created voxelization resulting from the Boolean intersection.
        :rtype: PointBasedVoxelization
        """
        self._check_other_voxelization_type(other_voxelization)
        self._check_other_voxelization_voxel_size(other_voxelization)

        return self.__class__(self.voxel_centers.intersection(other_voxelization.voxel_centers), self.voxel_size)

    def symmetric_difference(self, other_voxelization: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Perform a symmetric difference operation with another PointBasedVoxelization.

        :param other_voxelization: The PointBasedVoxelization to perform the symmetric difference with.
        :type other_voxelization: PointBasedVoxelization

        :return: A new PointBasedVoxelization resulting from the symmetric difference operation.
        :rtype: PointBasedVoxelization
        """
        self._check_other_voxelization_type(other_voxelization)
        self._check_other_voxelization_voxel_size(other_voxelization)

        return self.__class__(
            self.voxel_centers.symmetric_difference(other_voxelization.voxel_centers), self.voxel_size
        )

    def inverse(self) -> "PointBasedVoxelization":
        """
        Compute the inverse of the voxelization.

        :return: A new voxelization representing the inverse.
        :rtype: PointBasedVoxelization
        """
        inverted_voxel_matrix = self.to_matrix_based_voxelization().inverse()

        return self.__class__.from_matrix_based_voxelization(inverted_voxel_matrix)

    # FILLING METHODS
    def flood_fill(self, start: Point, fill_with: bool) -> "PointBasedVoxelization":
        """
        Perform a flood fill operation on the voxelization.

        :param start: The coordinates of the starting point for the flood fill.
        :type start: tuple[float, float, float]
        :param fill_with: The value to fill the voxels with during the operation.
        :type fill_with: bool

        :return: A new voxelization resulting from the flood fill operation.
        :rtype: PointBasedVoxelization
        """
        start = self._point_to_local_grid_index(start)
        voxel_matrix = self.to_matrix_based_voxelization()
        filled_voxel_matrix = voxel_matrix.flood_fill(start, fill_with)

        return self.from_matrix_based_voxelization(filled_voxel_matrix)

    def fill_outer_voxels(self) -> "PointBasedVoxelization":
        """
        Fill the outer voxels of the voxelization.

        :return: A new voxelization with outer voxels filled.
        :rtype: PointBasedVoxelization
        """
        return self.from_matrix_based_voxelization(self.to_matrix_based_voxelization().fill_outer_voxels())

    def fill_enclosed_voxels(self) -> "PointBasedVoxelization":
        """
        Fill the enclosed voxels of the voxelization.

        :return: A new voxelization with enclosed voxels filled.
        :rtype: VoxelizationType
        """
        return self.from_matrix_based_voxelization(self.to_matrix_based_voxelization().fill_enclosed_voxels())

    # MOVING METHODS
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
        :rtype: PointBasedVoxelization
        """
        rotation_matrix = self._rotation_matrix(axis, angle)
        voxel_array = np.array(list(self.voxel_centers)) - np.array([center.x, center.y, center.z])
        rotated_voxels = np.dot(voxel_array, rotation_matrix.T)
        rotated_voxels += np.array([center.x, center.y, center.z])

        intersecting_voxels = self._voxels_intersecting_voxels(rotated_voxels, self.voxel_size)

        return self.__class__(intersecting_voxels, self.voxel_size)

    def translation(self, offset: Vector3D):
        """
        Translate the voxelization by the specified offset.

        :param offset: The translation offset.
        :type offset: Vector3D

        :return: A new Voxelization object resulting from the translation.
        :rtype: PointBasedVoxelization
        """
        voxel_array = np.array(list(self.voxel_centers))
        translated_voxels = voxel_array + np.array([offset.x, offset.y, offset.z])

        intersecting_voxels = self._voxels_intersecting_voxels(translated_voxels, self.voxel_size)

        return self.__class__(intersecting_voxels, self.voxel_size)

    # HELPER METHODS
    def to_matrix_based_voxelization(self) -> "MatrixBasedVoxelization":
        """
        Convert the point based voxelization to a matrix based voxelization.

        :return: The matrix based voxelization.
        :rtype: MatrixBasedVoxelization
        """
        min_center = self.min_voxel_grid_center
        max_center = self.max_voxel_grid_center

        dim_x = round((max_center[0] - min_center[0]) / self.voxel_size + 1)
        dim_y = round((max_center[1] - min_center[1]) / self.voxel_size + 1)
        dim_z = round((max_center[2] - min_center[2]) / self.voxel_size + 1)

        indices = np.round((np.array(list(self.voxel_centers)) - min_center) / self.voxel_size).astype(int)

        matrix = np.zeros((dim_x, dim_y, dim_z), dtype=np.bool_)
        matrix[indices[:, 0], indices[:, 1], indices[:, 2]] = True

        return MatrixBasedVoxelization(matrix, min_center, self.voxel_size)

    def _point_to_local_grid_index(self, point: Point) -> Tuple[int, int, int]:
        """
        Convert a point to the local grid index within the voxelization.

        :param point: The point to convert.
        :type point: tuple[float, float, float]

        :return: The local grid index of the point.
        :rtype: Tuple[int, int, int]

        :raises ValueError: If the point is not within the voxelization's bounding box.
        """
        if not self.bounding_box.point_belongs(Point3D(*point)):
            raise ValueError("Point not in local voxel grid.")

        x_index = int((point[0] - self.bounding_box.xmin) // self.voxel_size)
        y_index = int((point[1] - self.bounding_box.ymin) // self.voxel_size)
        z_index = int((point[2] - self.bounding_box.zmin) // self.voxel_size)

        return x_index, y_index, z_index

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
        centers = np.around((unique_indices + 0.5) * voxel_size, DECIMALS)

        return set(tuple(center) for center in centers)


class MatrixBasedVoxelization(Voxelization):
    """Voxelization implemented as a 3D matrix."""

    def __init__(
        self,
        voxel_matrix: NDArray[np.bool_],
        min_voxel_grid_center: Point,
        voxel_size: float,
        name: str = "",
    ):
        """
        :param voxel_matrix: The voxel numpy matrix object representing the voxelization.
        :type voxel_matrix: np.ndarray[np.bool_, np.ndim == 3]
        :param voxel_size: The size of the voxel edges.
        :type voxel_size: float
        :param min_voxel_grid_center: Minimum voxel center point of the voxel grid matrix, i.e 'matrix[0][0][0]'.
        This point may not be a voxel of the voxelization, because it's the minimum center in each direction (X, Y, Z).
        :type min_voxel_grid_center: tuple[float, float, float]
        """
        self._check_voxel_size_number_of_decimals(voxel_size)

        self.matrix = voxel_matrix
        self._min_voxel_grid_center = min_voxel_grid_center

        Voxelization.__init__(self, voxel_size=voxel_size, name=name)

    @property
    def voxel_centers(self) -> Set[Point]:
        """
        Get the center point of each voxel.

        :return: The center point of each voxel.
        :rtype: set[tuple[float, float, float]]
        """
        indices = np.argwhere(self.matrix)
        voxel_centers = self._min_voxel_grid_center + indices * self.voxel_size

        return set(map(tuple, np.round(voxel_centers, DECIMALS)))

    def __eq__(self, other_voxelization: "MatrixBasedVoxelization") -> bool:
        """
        Check if two MatrixBasedVoxelization are equal.

        :param other_voxelization: Another MatrixBasedVoxelization to compare with.
        :type other_voxelization: MatrixBasedVoxelization

        :return: True if the MatrixBasedVoxelization are equal, False otherwise.
        :rtype: bool
        """
        return (
            self.voxel_size == other_voxelization.voxel_size
            and self._min_voxel_grid_center == other_voxelization._min_voxel_grid_center
            and np.array_equal(self.matrix, other_voxelization.matrix)
        )

    def __len__(self) -> int:
        """
        Get the number of voxels in the voxelization (i.e. the number of True value in the 3D voxel matrix).

        :return: The number of voxels in the voxelization.
        :rtype: int
        """
        return len(np.argwhere(self.matrix))

    @property
    def min_voxel_grid_center(self) -> Point:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, float, float]
        """
        return self._min_voxel_grid_center

    @property
    def max_voxel_grid_center(self) -> Point:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, float, float]
        """
        return tuple(
            np.round(
                np.array(self.min_voxel_grid_center) + (np.array(self.matrix.shape) - 1) * self.voxel_size,
                6,
            )
        )

    # CLASS METHODS
    @classmethod
    def from_shell(cls, shell: Shell3D, voxel_size: float, name: str = "") -> "MatrixBasedVoxelization":
        """
        Create a voxelization from a Shell3D.

        :param shell: The Shell3D to create the voxelization from.
        :type shell: Shell3D
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the Shell3D.
        :rtype: MatrixBasedVoxelization
        """
        triangles = Voxelization._shell_to_triangles(shell)
        matrix, matrix_origin_center = triangles_to_voxel_matrix(triangles, voxel_size)

        return cls(matrix, matrix_origin_center, voxel_size, name).crop_matrix()

    @classmethod
    def from_volume_model(
        cls, volume_model: VolumeModel, voxel_size: float, name: str = ""
    ) -> "MatrixBasedVoxelization":
        """
        Create a voxelization from a VolumeModel.

        :param volume_model: The VolumeModel to create the voxelization from.
        :type volume_model: VolumeModel
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the VolumeModel.
        :rtype: MatrixBasedVoxelization
        """
        triangles = Voxelization._volume_model_to_triangles(volume_model)
        matrix, matrix_origin_center = triangles_to_voxel_matrix(triangles, voxel_size)

        return cls(matrix, matrix_origin_center, voxel_size, name).crop_matrix()

    @classmethod
    def from_point_voxelization(cls, voxelization: "PointBasedVoxelization") -> "MatrixBasedVoxelization":
        return voxelization.to_matrix_based_voxelization()

    # BOOLEAN OPERATIONS
    def union(self, other_voxelization: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform a union operation with another MatrixBasedVoxelization.

        :param other_voxelization: The MatrixBasedVoxelization to perform the union with.
        :type other_voxelization: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the union operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other_voxelization, np.logical_or)

    def difference(self, other_voxelization: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform a difference operation with another voxelization.

        :param other_voxelization: The MatrixBasedVoxelization to perform the difference with.
        :type other_voxelization: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the difference operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other_voxelization, lambda a, b: np.logical_and(a, np.logical_not(b)))

    def intersection(self, other_voxelization: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform an intersection operation with another MatrixBasedVoxelization.

        :param other_voxelization: The MatrixBasedVoxelization to perform the intersection with.
        :type other_voxelization: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the intersection operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other_voxelization, np.logical_and)

    def symmetric_difference(self, other_voxelization: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform a symmetric difference operation with another MatrixBasedVoxelization.

        :param other_voxelization: The MatrixBasedVoxelization to perform the symmetric difference with.
        :type other_voxelization: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the symmetric difference operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other_voxelization, np.logical_xor)

    def inverse(self) -> "MatrixBasedVoxelization":
        """
        Compute the inverse of the voxelization.

        :return: A new voxelization representing the inverse.
        :rtype: MatrixBasedVoxelization
        """
        inverted_matrix = np.logical_not(self.matrix)
        return self.__class__(inverted_matrix, self._min_voxel_grid_center, self.voxel_size)

    def flood_fill(self, start: Tuple[int, int, int], fill_with: bool) -> "MatrixBasedVoxelization":
        """
        Perform a flood fill operation on the voxelization.

        :param start: The indexes of the starting voxel in the 3D matrix for the flood fill.
        :type start: tuple[int, int, int]
        :param fill_with: The value to fill the voxels with during the operation.
        :type fill_with: bool

        :return: A new voxelization resulting from the flood fill operation.
        :rtype: MatrixBasedVoxelization
        """
        return self.__class__(
            flood_fill_matrix_3d(self.matrix, start, fill_with), self._min_voxel_grid_center, self.voxel_size
        )

    def fill_outer_voxels(self) -> "MatrixBasedVoxelization":
        """
        Fill the outer voxels of the voxelization.

        :return: A new voxelization with outer voxels filled.
        :rtype: MatrixBasedVoxelization
        """
        # pylint: disable=protected-access

        expanded_voxel_matrix = self._expand()
        outer_filled_expanded_voxel_matrix = expanded_voxel_matrix.flood_fill((0, 0, 0), True)
        outer_filled_voxel_matrix = outer_filled_expanded_voxel_matrix._reduce()

        return outer_filled_voxel_matrix

    def fill_enclosed_voxels(self) -> "MatrixBasedVoxelization":
        """
        Fill the enclosed voxels of the voxelization.

        :return: A new voxelization with enclosed voxels filled.
        :rtype: MatrixBasedVoxelization
        """
        outer_filled_voxel_matrix = self.fill_outer_voxels()
        inner_filled_voxel_matrix = self + outer_filled_voxel_matrix.inverse()

        return inner_filled_voxel_matrix

    # HELPER METHODS
    def to_point_voxelization(self) -> "PointBasedVoxelization":
        """
        Convert the MatrixBasedVoxelization to a PointBasedVoxelization.

        :return: A PointBasedVoxelization representation of the current voxelization.
        :rtype: PointBasedVoxelization
        """
        return PointBasedVoxelization(self.voxel_centers, self.voxel_size)

    def _expand(self) -> "MatrixBasedVoxelization":
        """
        Expand the voxelization matrix by adding a single layer of False voxels around the existing matrix.

        :return: A new MatrixBasedVoxelization with an expanded voxelization matrix.
        :rtype: MatrixBasedVoxelization
        """
        current_shape = self.matrix.shape
        new_shape = tuple(dim + 2 for dim in current_shape)
        expanded_matrix = np.zeros(new_shape, dtype="bool")
        slices = tuple(slice(1, -1) for _ in current_shape)
        expanded_matrix[slices] = self.matrix.copy()

        return self.__class__(
            expanded_matrix,
            tuple(np.round(np.array(self.min_voxel_grid_center) - self.voxel_size, DECIMALS)),
            self.voxel_size,
        )

    def _reduce(self) -> "MatrixBasedVoxelization":
        """
        Reduce the size of the voxelization matrix by removing a single layer of voxels around the edges.

        :return: A new MatrixBasedVoxelization with a reduced voxelization matrix.
        :rtype: MatrixBasedVoxelization
        """
        current_shape = self.matrix.shape
        slices = tuple(slice(1, -1) for _ in current_shape)
        reduced_matrix = self.matrix.copy()[slices]

        return self.__class__(
            reduced_matrix,
            tuple(np.round(np.array(self.min_voxel_grid_center) + self.voxel_size, DECIMALS)),
            self.voxel_size,
        )

    def _logical_operation(self, other_voxelization: "MatrixBasedVoxelization", logical_operation):
        """
        Perform a logical operation (e.g., union, intersection, etc.) between two voxelizations.

        :param other_voxelization: The other MatrixBasedVoxelization to perform the operation with.
        :type other_voxelization: MatrixBasedVoxelization
        :param logical_operation: The logical operation function to apply (e.g., np.logical_or).

        :return: A new MatrixBasedVoxelization resulting from the logical operation.
        :rtype: MatrixBasedVoxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Voxel sizes must be the same to perform boolean operations.")

        self_min, self_max = np.array(self.min_voxel_grid_center), np.array(self.min_voxel_grid_center) + 1
        other_min, other_max = (
            np.array(other_voxelization.min_voxel_grid_center),
            np.array(other_voxelization.min_voxel_grid_center) + 1,
        )

        global_min = np.min([self_min, other_min], axis=0)
        global_max = np.max([self_max, other_max], axis=0)

        new_shape = np.round((global_max - global_min) / self.voxel_size, DECIMALS).astype(int)

        new_self = np.zeros(new_shape, dtype=bool)
        new_other = np.zeros(new_shape, dtype=bool)

        self_start = np.round((self_min - global_min) / self.voxel_size, DECIMALS).astype(int)
        other_start = np.round((other_min - global_min) / self.voxel_size, DECIMALS).astype(int)

        new_self[
            self_start[0] : self_start[0] + self.matrix.shape[0],
            self_start[1] : self_start[1] + self.matrix.shape[1],
            self_start[2] : self_start[2] + self.matrix.shape[2],
        ] = self.matrix

        new_other[
            other_start[0] : other_start[0] + other_voxelization.matrix.shape[0],
            other_start[1] : other_start[1] + other_voxelization.matrix.shape[1],
            other_start[2] : other_start[2] + other_voxelization.matrix.shape[2],
        ] = other_voxelization.matrix

        result_matrix = logical_operation(new_self, new_other)

        return self.__class__(result_matrix, tuple(global_min), self.voxel_size).crop_matrix()

    def crop_matrix(self) -> "MatrixBasedVoxelization":
        """
        Crop the voxel matrix to the smallest possible size.

        :return: The MatrixBasedVoxelization with cropped voxel matrix.
        :rtype: MatrixBasedVoxelization
        """
        # Find the indices of the True voxels
        true_voxels = np.argwhere(self.matrix)

        if len(true_voxels) == 0:
            # Can't crop a matrix of False values
            return self

        # Find the minimum and maximum indices along each axis
        min_voxel_coords, max_voxel_coords = np.round(np.min(true_voxels, axis=0), DECIMALS), np.round(
            np.max(true_voxels, axis=0), 6
        )

        # Crop the matrix to the smallest possible size
        cropped_matrix = self.matrix[
            min_voxel_coords[0] : max_voxel_coords[0] + 1,
            min_voxel_coords[1] : max_voxel_coords[1] + 1,
            min_voxel_coords[2] : max_voxel_coords[2] + 1,
        ]

        # Calculate new matrix_origin_center
        new_origin_center = np.round(self._min_voxel_grid_center + min_voxel_coords * self.voxel_size, DECIMALS)

        return self.__class__(cropped_matrix, tuple(new_origin_center), self.voxel_size)


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
        return cls(set(map(tuple, np.round(pixel_centers, DECIMALS))), pixel_size)

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

    def __init__(self, numpy_pixel_matrix: NDArray[np.bool_]):
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
        # pylint: disable=protected-access

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
