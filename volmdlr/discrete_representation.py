"""
Class for discrete representations of volmdlr models (voxelization for 3D geometries, pixelization for 2D geometries).
"""
# pylint: disable=no-name-in-module

import warnings
from abc import ABC, abstractmethod
from typing import List, Set, Tuple, TypeVar, Dict, Any

import matplotlib.pyplot as plt
import numpy as np
from dessia_common.core import DessiaObject, PhysicalObject
from dessia_common.serialization import JsonSerializable
from matplotlib import patches
from numpy.typing import NDArray

from volmdlr import Point2D, Point3D, Vector3D
from volmdlr.core import BoundingBox, BoundingRectangle, VolumeModel
from volmdlr.discrete_representation_compiled import (
    flood_fill_matrix_2d,
    flood_fill_matrix_3d,
    line_segments_to_pixels,
    triangles_to_voxel_matrix,
    voxel_triangular_faces,
)
from volmdlr.edges import LineSegment2D
from volmdlr.faces import Triangle3D
from volmdlr.shells import ClosedTriangleShell3D, Shell3D
from volmdlr.wires import ClosedPolygon2D

# CUSTOM TYPES
_Point3D = Tuple[float, float, float]
_Triangle3D = Tuple[_Point3D, _Point3D, _Point3D]

_Point2D = Tuple[float, float]
_Segment2D = Tuple[_Point2D, _Point2D]


# GLOBAL VARIABLE
DECIMALS = 9  # Used to round numbers and avoid floating point arithmetic imprecision


class DiscreteRepresentation(ABC):
    """
    Abstract base class for discrete representation in any dimension.

    It is used for:
        - Voxelization in 3D
        - Pixelization in 2D

    Any discrete representation follows the same approach:
    The representation is defined on an implicit grid, where each element is defined by its size.
    The implicit grid consists of axis-aligned elements with a given size, ensuring that the global
    origin (0, ..., 0) is always a corner point of an element.

    For example, in 1D with an element size of `t`, the set of elements (defined by minimum and maximum points)
    is: {i ∈ N, t ∈ R | (i * t, (i+1) * t)}
    The corresponding set of element centers is: {i ∈ N, t ∈ R | (i + 0.5) * t}

    This approach enables consistent representation across the space, facilitating fast Boolean operations.
    """

    DiscreteRepresentationType = TypeVar("DiscreteRepresentationType", bound="DiscreteRepresentation")

    def __init__(self, element_size: float):
        """
        Initialize the discrete representation.

        :param element_size: The element size.
        :type element_size: float
        """
        self.element_size = element_size

    def __str__(self):
        """
        Return a custom string representation of the discrete representation.

        :return: A string representation of the discrete representation.
        :rtype: str
        """
        return f"{self.__class__}: element size={self.element_size}, number of elements={len(self)}"

    @abstractmethod
    def __eq__(self, other: DiscreteRepresentationType) -> bool:
        """
        Check if two discrete representations are equal.

        :param other: Another discrete representation to compare with.
        :type other: DiscreteRepresentationType

        :return: True if the discrete representations are equal, False otherwise.
        :rtype: bool
        """

    @abstractmethod
    def __len__(self) -> int:
        """
        Get the number of elements in the discrete representation.

        :return: The number of elements in the discrete representation.
        :rtype: int
        """

    @abstractmethod
    def _get_element_centers(self) -> Set[Tuple]:
        """
        Get the center point of each element.

        :return: The center point of each element.
        :rtype: set[tuple[float, ...]]
        """

    @property
    @abstractmethod
    def min_grid_center(self) -> Tuple:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, ...]
        """

    @property
    @abstractmethod
    def max_grid_center(self) -> Tuple:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, ...]
        """

    # BOOLEAN OPERATIONS
    def __add__(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Overloaded '+' operator for performing a union operation.

        :param other: The discrete representation to perform the union with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the union operation.
        :rtype: DiscreteRepresentationType
        """
        return self.union(other)

    def __sub__(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Overloaded '-' operator for performing a difference operation.

        :param other: The discrete representation to perform the difference with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the difference operation.
        :rtype: DiscreteRepresentationType
        """
        return self.difference(other)

    def __and__(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Overloaded '&' operator for performing an intersection operation.

        :param other: The discrete representation to perform the intersection with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the intersection operation.
        :rtype: DiscreteRepresentationType
        """
        return self.intersection(other)

    def __xor__(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Overloaded '^' operator for performing a symmetric difference operation.

        :param other: The discrete representation to perform the symmetric difference with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the symmetric difference operation.
        :rtype: DiscreteRepresentationType
        """
        return self.symmetric_difference(other)

    def __invert__(self) -> DiscreteRepresentationType:
        """
        Overloaded '~' operator for computing the inverse.

        :return: A new discrete representation representing the inverse.
        :rtype: DiscreteRepresentation
        """
        return self.inverse()

    @abstractmethod
    def union(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform a union operation with another discrete representation.

        :param other: The discrete representation to perform the union with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the union operation.
        :rtype: DiscreteRepresentationType
        """

    @abstractmethod
    def difference(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform a difference operation with another discrete representation.

        :param other: The discrete representation to perform the difference with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the difference operation.
        :rtype: DiscreteRepresentationType
        """

    @abstractmethod
    def intersection(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform an intersection operation with another discrete representation.

        :param other: The discrete representation to perform the intersection with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the intersection operation.
        :rtype: DiscreteRepresentationType
        """

    def symmetric_difference(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform a symmetric difference operation with another discrete representation.

        :param other: The discrete representation to perform the symmetric difference with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the symmetric difference operation.
        :rtype: DiscreteRepresentationType
        """

    @abstractmethod
    def inverse(self) -> DiscreteRepresentationType:
        """
        Compute the inverse of the discrete representation.

        :return: A new discrete representation representing the inverse.
        :rtype: DiscreteRepresentationType
        """

    def interference(self, other: DiscreteRepresentationType) -> float:
        """
        Compute the percentage of interference between two discrete representations.

        :param other: The other discrete representation to compute interference with.
        :type other: DiscreteRepresentation

        :return: The percentage of interference between the two discrete representations.
        :rtype: float
        """
        return len(self.intersection(other)) / len(self.union(other))

    def is_intersecting(self, other: DiscreteRepresentationType) -> bool:
        """
        Check if two discrete representations are intersecting.

        :param other: The other discrete representation to check if there is an intersection with.
        :type other: DiscreteRepresentationType

        :return: True if the discrete representations are intersecting, False otherwise.
        :rtype: bool
        """
        intersection = self.intersection(other)

        return len(intersection) > 0

    # FILLING METHODS
    @abstractmethod
    def flood_fill(self, start, fill_with: bool) -> DiscreteRepresentationType:
        """
        Perform a flood fill operation on the discrete representation.

        :param start: The starting point for the flood fill.
        :param fill_with: The value to fill the elements with during the operation.
        :type fill_with: bool

        :return: A new discrete representation resulting from the flood fill operation.
        :rtype: DiscreteRepresentationType
        """

    @abstractmethod
    def _fill_outer_elements(self) -> DiscreteRepresentationType:
        """
        Fill the outer elements of the discrete representation.

        :return: A new discrete representation with outer elements filled.
        :rtype: DiscreteRepresentationType
        """

    @abstractmethod
    def _fill_enclosed_elements(self) -> DiscreteRepresentationType:
        """
        Fill the enclosed elements of the discrete representation.

        :return: A new discrete representation with enclosed elements filled.
        :rtype: DiscreteRepresentationType
        """

    # HELPER METHODS
    @staticmethod
    def check_center_is_in_implicit_grid(element_center: Tuple[float, ...], element_size: float) -> bool:
        """
        Check if a given element center point is an element center of the implicit grid, defined by element_size.

        :param element_center: The element center point to check.
        :type element_center: tuple[float, ...]
        :param element_size: The element edges size.
        :type element_size: float

        :return: True if the given element center point is an element center of the implicit grid, False otherwise.
        :rtype: bool
        """
        for coord in element_center:
            if not round((coord - 0.5 * element_size) / element_size, DECIMALS).is_integer():
                return False

        return True

    def is_hollow(self) -> bool:
        """
        Check if the discrete representation is hollow.

        A hollow discrete representation is one that has enclosed elements that are not filled.

        :return: True if the discrete representation is hollow, False otherwise.
        :rtype: bool
        """
        return self._fill_enclosed_elements() == self

    @staticmethod
    def _check_element_size_number_of_decimals(element_size: float):
        """
        Check the number of decimal places in the element size.

        If the element size has more decimal places than a specified threshold (DECIMALS),
        a warning is issued, as some functions may not work as intended.

        This is due to the use of rounding functions to avoid floating point arithmetic imprecision.

        :param element_size: The size of the elements.
        :type element_size: float

        :raises ValueError: If element_size is not a float.
        """
        if isinstance(element_size, float):
            decimal_part = abs(element_size - int(element_size))
            if decimal_part == 0:
                pass
            else:
                decimals = len(str(decimal_part).split(".")[1])
                if decimals >= DECIMALS:
                    warnings.warn(
                        f"""Element size has too many decimals: some functions may not work as intended.
                        Consider using an element size with less than {DECIMALS}."""
                    )
        else:
            raise ValueError("Element size is not a float")

    def _check_other_type(self, other):
        """
        Check if the provided 'other' is an instance of the same DiscreteRepresentation subclass.

        :param other: Another discrete representation to be checked.
        :type other: DiscreteRepresentation

        :raises ValueError: If 'other' is not an instance of the same DiscreteRepresentation subclass.
        """
        if not isinstance(other, self.__class__):
            raise ValueError(f"'other' must be an instance of '{self.__class__.__name__}'")

    def _check_other_element_size(self, other: DiscreteRepresentationType):
        """
        Check if the provided 'other' has the same element size.

        :param other: Another discrete representation to be checked.
        :type other: DiscreteRepresentation

        :raises ValueError: If 'other' has not the same element size.
        """
        if not self.element_size == other.element_size:
            raise ValueError(f"Both {self.__class__} must have same element size to perform this operation.")


class Voxelization(DiscreteRepresentation, PhysicalObject):
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
        """
        Initialize the voxelization.

        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the voxelization.
        :type name: str, optional
        """
        DiscreteRepresentation.__init__(self, element_size=voxel_size)
        PhysicalObject.__init__(self, name=name)

    @property
    def voxel_size(self) -> float:
        """
        Get the voxel size.

        :return: The voxel size.
        :rtype: float
        """
        return self.element_size

    def get_voxel_centers(self) -> Set[_Point3D]:
        """
        Get the center point of each voxel.

        :return: The center point of each voxel.
        :rtype: set[tuple[float, float, float]]
        """
        return self._get_element_centers()

    @property
    def volume(self) -> float:
        """
        Calculate the volume of the voxelization.

        :return: The volume of the voxelization.
        :rtype: float
        """
        return len(self) * self.voxel_size**3

    @property
    def bounding_box(self) -> BoundingBox:
        """
        Get the bounding box of the voxelization.

        :return: The bounding box of the voxelization.
        :rtype: BoundingBox
        """
        min_point = np.round((np.array([self.min_grid_center]) - self.voxel_size)[0], DECIMALS)
        max_point = np.round((np.array([self.max_grid_center]) + self.voxel_size)[0], DECIMALS)

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

    # FILLING METHODS
    def fill_outer_voxels(self) -> VoxelizationType:
        """
        Fill the outer voxels of the voxelization.

        :return: A new voxelization with outer voxels filled.
        :rtype: VoxelizationType
        """
        return self._fill_outer_elements()

    def fill_enclosed_voxels(self) -> VoxelizationType:
        """
        Fill the enclosed voxels of the voxelization.

        :return: A new voxelization with enclosed voxels filled.
        :rtype: VoxelizationType
        """
        return self._fill_enclosed_elements()

    # DISPLAY METHODS
    def to_triangles(self) -> Set[_Triangle3D]:
        """
        Convert the voxelization to triangles for display purpose.

        Only the relevant faces are returned (i.e. the faces that are not at the interface of two different voxel,
        i.e. the faces that are only present once in the list of triangles representing the triangulated voxels).

        :return: The triangles representing the voxelization.
        :rtype: set[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangles = set()

        for voxel in self.get_voxel_centers():
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
    def _shell_to_triangles(shell: Shell3D) -> List[_Triangle3D]:
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
    def _volume_model_to_triangles(volume_model: VolumeModel) -> List[_Triangle3D]:
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

    def __init__(self, voxel_centers: Set[_Point3D], voxel_size: float, name: str = ""):
        """
        Initialize the PointBasedVoxelization.

        :param voxel_centers: The set of points representing voxel centers.
        :type voxel_centers: set[tuple[float, float, float]]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the voxelization.
        :type name: str, optional
        """
        self._check_element_size_number_of_decimals(voxel_size)

        self.voxel_centers = voxel_centers

        Voxelization.__init__(self, voxel_size=voxel_size, name=name)

    def _get_element_centers(self) -> Set[_Point3D]:
        """
        Get the center point of each voxel.

        :return: The center point of each voxel.
        :rtype: set[tuple[float, float, float]]
        """
        return self.voxel_centers

    def __eq__(self, other: "PointBasedVoxelization") -> bool:
        """
        Check if two voxelizations are equal.

        :param other: Another voxelization to compare with.
        :type other: PointBasedVoxelization

        :return: True if the voxelizations are equal, False otherwise.
        :rtype: bool
        """
        self._check_other_type(other)

        return self.voxel_centers == other.voxel_centers and self.voxel_size == other.voxel_size

    def __len__(self) -> int:
        """
        Get the number of voxels in the voxelization.

        :return: The number of voxels in the voxelization (i.e. the number of voxel centers).
        :rtype: int
        """
        return len(self.voxel_centers)

    @property
    def min_grid_center(self) -> _Point3D:
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
    def max_grid_center(self) -> _Point3D:
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
        voxels = MatrixBasedVoxelization.from_shell(shell, voxel_size).get_voxel_centers()

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
        voxels = MatrixBasedVoxelization.from_volume_model(volume_model, voxel_size).get_voxel_centers()

        return cls(voxel_centers=voxels, voxel_size=voxel_size, name=name)

    @classmethod
    def from_matrix_based_voxelization(
        cls, matrix_based_voxelization: "MatrixBasedVoxelization",
    ) -> "PointBasedVoxelization":
        """
        Create a PointBasedVoxelization object from a MatrixBasedVoxelization.

        :param matrix_based_voxelization: The MatrixBasedVoxelization object representing the voxelization.
        :type matrix_based_voxelization: MatrixBasedVoxelization
        :param name: object's name.

        :return: A PointBasedVoxelization object created from the MatrixBasedVoxelization.
        :rtype: PointBasedVoxelization
        """
        if not cls.check_center_is_in_implicit_grid(
            matrix_based_voxelization.min_grid_center, matrix_based_voxelization.voxel_size
        ):
            warnings.warn(
                """This matrix based voxelization is not defined in the implicit grid defined by the voxel_size.
            Some methods like boolean operation or interference computing may not work as expected."""
            )

        return matrix_based_voxelization.to_point_based_voxelization()

    # BOOLEAN OPERATIONS
    def union(self, other: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Perform a union operation with another PointBasedVoxelization.

        :param other: The PointBasedVoxelization to perform the union with.
        :type other: PointBasedVoxelization

        :return: A new PointBasedVoxelization resulting from the union operation.
        :rtype: PointBasedVoxelization
        """
        self._check_other_type(other)
        self._check_other_element_size(other)

        return self.__class__(self.voxel_centers.union(other.voxel_centers), self.voxel_size)

    def difference(self, other: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Perform a difference operation with another PointBasedVoxelization.

        :param other: The PointBasedVoxelization to perform the difference with.
        :type other: PointBasedVoxelization

        :return: A new PointBasedVoxelization resulting from the difference operation.
        :rtype: PointBasedVoxelization
        """
        self._check_other_type(other)
        self._check_other_element_size(other)

        return self.__class__(self.voxel_centers.difference(other.voxel_centers), self.voxel_size)

    def intersection(self, other: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Perform an intersection operation with another PointBasedVoxelization.

        :param other: The PointBasedVoxelization to perform the intersection with.
        :type other: PointBasedVoxelization

        :return: A new PointBasedVoxelization resulting from the intersection operation.
        :rtype: PointBasedVoxelization
        """
        self._check_other_type(other)
        self._check_other_element_size(other)

        return self.__class__(self.voxel_centers.intersection(other.voxel_centers), self.voxel_size)

    def symmetric_difference(self, other: "PointBasedVoxelization") -> "PointBasedVoxelization":
        """
        Perform a symmetric difference operation with another PointBasedVoxelization.

        :param other: The PointBasedVoxelization to perform the symmetric difference with.
        :type other: PointBasedVoxelization

        :return: A new PointBasedVoxelization resulting from the symmetric difference operation.
        :rtype: PointBasedVoxelization
        """
        self._check_other_type(other)
        self._check_other_element_size(other)

        return self.__class__(self.voxel_centers.symmetric_difference(other.voxel_centers), self.voxel_size)

    def inverse(self) -> "PointBasedVoxelization":
        """
        Compute the inverse of the voxelization.

        :return: A new voxelization representing the inverse.
        :rtype: PointBasedVoxelization
        """
        inverted_voxel_matrix = self.to_matrix_based_voxelization().inverse()

        return self.__class__.from_matrix_based_voxelization(inverted_voxel_matrix)

    # FILLING METHODS
    def flood_fill(self, start: _Point3D, fill_with: bool) -> "PointBasedVoxelization":
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

    def _fill_outer_elements(self) -> "PointBasedVoxelization":
        """
        Fill the outer voxels of the voxelization.

        :return: A new voxelization with outer voxels filled.
        :rtype: PointBasedVoxelization
        """
        return self.from_matrix_based_voxelization(self.to_matrix_based_voxelization().fill_outer_voxels())

    def _fill_enclosed_elements(self) -> "PointBasedVoxelization":
        """
        Fill the enclosed voxels of the voxelization.

        :return: A new voxelization with enclosed voxels filled.
        :rtype: PointBasedVoxelization
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

    # SERIALIZATION
    def to_dict(
        self, use_pointers: bool = True, memo=None, path: str = "#", id_method=True, id_memo=None
    ) -> JsonSerializable:
        """Specific 'to_dict' method to allow serialization of a set."""
        dict_ = self.base_dict()

        dict_["voxel_centers"] = list(self.voxel_centers)
        dict_["voxel_size"] = self.voxel_size
        dict_["name"] = self.name

        return dict_

    @classmethod
    def dict_to_object(
        cls,
        dict_: JsonSerializable,
        force_generic: bool = False,
        global_dict=None,
        pointers_memo: Dict[str, Any] = None,
        path: str = "#",
    ) -> "PointBasedVoxelization":
        """Specific 'dict_to_object' method to allow deserialization of a set."""

        voxel_centers = set(tuple(voxel_center) for voxel_center in dict_["voxel_centers"])
        voxel_size = dict_["voxel_size"]
        name = dict_["name"]

        return cls(voxel_centers, voxel_size, name)

    # HELPER METHODS
    def to_matrix_based_voxelization(self) -> "MatrixBasedVoxelization":
        """
        Convert the point based voxelization to a matrix based voxelization.

        :return: The matrix based voxelization.
        :rtype: MatrixBasedVoxelization
        """
        min_center = self.min_grid_center
        max_center = self.max_grid_center

        dim_x = round((max_center[0] - min_center[0]) / self.voxel_size + 1)
        dim_y = round((max_center[1] - min_center[1]) / self.voxel_size + 1)
        dim_z = round((max_center[2] - min_center[2]) / self.voxel_size + 1)

        indices = np.round((np.array(list(self.voxel_centers)) - min_center) / self.voxel_size).astype(int)

        matrix = np.zeros((dim_x, dim_y, dim_z), dtype=np.bool_)
        matrix[indices[:, 0], indices[:, 1], indices[:, 2]] = True

        return MatrixBasedVoxelization(matrix, min_center, self.voxel_size, self.name)

    def _point_to_local_grid_index(self, point: _Point3D) -> Tuple[int, int, int]:
        """
        Convert a point to the local grid index within the voxelization.

        :param point: The point to convert.
        :type point: tuple[float, float, float]

        :return: The local grid index of the point.
        :rtype: Tuple[int, int, int]

        :raises ValueError: If the point is not within the bounding box of the voxelization.
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
    def _voxels_intersecting_voxels(voxel_centers_array: np.ndarray, voxel_size: float) -> Set[_Point3D]:
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
        min_grid_center: _Point3D,
        voxel_size: float,
        name: str = "",
    ):
        """
        Initialize a MatrixBasedVoxelization.

        :param voxel_matrix: The voxel numpy matrix object representing the voxelization.
        :type voxel_matrix: np.ndarray[np.bool_, np.ndim == 3]
        :param voxel_size: The size of the voxel edges.
        :type voxel_size: float
        :param min_grid_center: Minimum voxel center point of the voxel grid matrix, i.e 'matrix[0][0][0]'.
        This point may not be a voxel of the voxelization, because it's the minimum center in each direction (X, Y, Z).
        :type min_grid_center: tuple[float, float, float]
        """
        self._check_element_size_number_of_decimals(voxel_size)

        self.matrix = voxel_matrix
        self._min_grid_center = min_grid_center

        Voxelization.__init__(self, voxel_size=voxel_size, name=name)

    def _get_element_centers(self) -> Set[_Point3D]:
        """
        Get the center point of each voxel.

        :return: The center point of each voxel.
        :rtype: set[tuple[float, float, float]]
        """
        indices = np.argwhere(self.matrix)
        voxel_centers = self.min_grid_center + indices * self.voxel_size

        return set(map(tuple, np.round(voxel_centers, DECIMALS)))

    def __eq__(self, other: "MatrixBasedVoxelization") -> bool:
        """
        Check if two MatrixBasedVoxelization are equal.

        :param other: Another MatrixBasedVoxelization to compare with.
        :type other: MatrixBasedVoxelization

        :return: True if the MatrixBasedVoxelization are equal, False otherwise.
        :rtype: bool
        """
        return (
            self.voxel_size == other.voxel_size
            and self.min_grid_center == other.min_grid_center
            and np.array_equal(self.matrix, other.matrix)
        )

    def __len__(self) -> int:
        """
        Get the number of voxels in the voxelization (i.e. the number of True value in the 3D voxel matrix).

        :return: The number of voxels in the voxelization.
        :rtype: int
        """
        return len(np.argwhere(self.matrix))

    @property
    def min_grid_center(self) -> _Point3D:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, float, float]
        """
        return self._min_grid_center

    @property
    def max_grid_center(self) -> _Point3D:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, float, float]
        """
        return tuple(
            np.round(
                np.array(self.min_grid_center) + (np.array(self.matrix.shape) - 1) * self.voxel_size,
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
    def from_point_based_voxelization(
        cls, point_based_voxelization: "PointBasedVoxelization"
    ) -> "MatrixBasedVoxelization":
        """
        Create a MatrixBasedVoxelization from a PointBasedVoxelization.

        :param point_based_voxelization: The PointBasedVoxelization to create the MatrixBasedVoxelization from.
        :type point_based_voxelization: PointBasedVoxelization

        :return: A MatrixBasedVoxelization created from the PointBasedVoxelization.
        :rtype: MatrixBasedVoxelization
        """
        return point_based_voxelization.to_matrix_based_voxelization()

    # BOOLEAN OPERATIONS
    def union(self, other: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform a union operation with another MatrixBasedVoxelization.

        :param other: The MatrixBasedVoxelization to perform the union with.
        :type other: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the union operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other, np.logical_or)

    def difference(self, other: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform a difference operation with another voxelization.

        :param other: The MatrixBasedVoxelization to perform the difference with.
        :type other: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the difference operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other, lambda a, b: np.logical_and(a, np.logical_not(b)))

    def intersection(self, other: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform an intersection operation with another MatrixBasedVoxelization.

        :param other: The MatrixBasedVoxelization to perform the intersection with.
        :type other: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the intersection operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other, np.logical_and)

    def symmetric_difference(self, other: "MatrixBasedVoxelization") -> "MatrixBasedVoxelization":
        """
        Perform a symmetric difference operation with another MatrixBasedVoxelization.

        :param other: The MatrixBasedVoxelization to perform the symmetric difference with.
        :type other: MatrixBasedVoxelization

        :return: A new MatrixBasedVoxelization resulting from the symmetric difference operation.
        :rtype: MatrixBasedVoxelization
        """
        return self._logical_operation(other, np.logical_xor)

    def inverse(self) -> "MatrixBasedVoxelization":
        """
        Compute the inverse of the voxelization.

        :return: A new voxelization representing the inverse.
        :rtype: MatrixBasedVoxelization
        """
        inverted_matrix = np.logical_not(self.matrix)
        return self.__class__(inverted_matrix, self.min_grid_center, self.voxel_size)

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
            flood_fill_matrix_3d(self.matrix, start, fill_with), self.min_grid_center, self.voxel_size
        )

    def _fill_outer_elements(self) -> "MatrixBasedVoxelization":
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

    def _fill_enclosed_elements(self) -> "MatrixBasedVoxelization":
        """
        Fill the enclosed voxels of the voxelization.

        :return: A new voxelization with enclosed voxels filled.
        :rtype: MatrixBasedVoxelization
        """
        outer_filled_voxel_matrix = self.fill_outer_voxels()
        inner_filled_voxel_matrix = self + outer_filled_voxel_matrix.inverse()

        return inner_filled_voxel_matrix

    # SERIALIZATION
    def to_dict(
        self, use_pointers: bool = True, memo=None, path: str = "#", id_method=True, id_memo=None
    ) -> JsonSerializable:
        """Specific 'to_dict' method to allow serialization of a numpy array."""
        dict_ = self.base_dict()

        dict_["matrix"] = self.matrix.tolist()
        dict_["min_grid_center"] = self.min_grid_center
        dict_["voxel_size"] = self.voxel_size
        dict_["name"] = self.name

        return dict_

    @classmethod
    def dict_to_object(
        cls,
        dict_: JsonSerializable,
        force_generic: bool = False,
        global_dict=None,
        pointers_memo: Dict[str, Any] = None,
        path: str = "#",
    ) -> "MatrixBasedVoxelization":
        """Specific 'dict_to_object' method to allow deserialization of a numpy array."""

        matrix = np.array(dict_["matrix"])
        min_grid_center = tuple(dict_["min_grid_center"])
        voxel_size = dict_["voxel_size"]
        name = dict_["name"]

        return cls(matrix, min_grid_center, voxel_size, name)

    # HELPER METHODS
    def to_point_based_voxelization(self) -> "PointBasedVoxelization":
        """
        Convert the MatrixBasedVoxelization to a PointBasedVoxelization.

        :return: A PointBasedVoxelization representation of the current voxelization.
        :rtype: PointBasedVoxelization
        """
        return PointBasedVoxelization(self.get_voxel_centers(), self.voxel_size, self.name)

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
            tuple(np.round(np.array(self.min_grid_center) - self.voxel_size, DECIMALS)),
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
            tuple(np.round(np.array(self.min_grid_center) + self.voxel_size, DECIMALS)),
            self.voxel_size,
        )

    def _logical_operation(self, other: "MatrixBasedVoxelization", logical_operation) -> "MatrixBasedVoxelization":
        """
        Perform a logical operation (e.g., union, intersection, etc.) between two voxelizations.

        :param other: The other MatrixBasedVoxelization to perform the operation with.
        :type other: MatrixBasedVoxelization
        :param logical_operation: The logical operation function to apply (e.g., np.logical_or).

        :return: A new MatrixBasedVoxelization resulting from the logical operation.
        :rtype: MatrixBasedVoxelization
        """
        if self.voxel_size != other.voxel_size:
            raise ValueError("Voxel sizes must be the same to perform boolean operations.")

        self_min, self_max = np.array(self.min_grid_center), np.array(self.min_grid_center) + 1
        other_min, other_max = (
            np.array(other.min_grid_center),
            np.array(other.min_grid_center) + 1,
        )

        global_min = np.min([self_min, other_min], axis=0)
        global_max = np.max([self_max, other_max], axis=0)

        new_shape = np.round((global_max - global_min) / self.voxel_size, DECIMALS).astype(int) + 2

        new_self = np.zeros(new_shape, dtype=np.bool_)
        new_other = np.zeros(new_shape, dtype=np.bool_)

        self_start = np.round((self_min - global_min) / self.voxel_size, DECIMALS).astype(int)
        other_start = np.round((other_min - global_min) / self.voxel_size, DECIMALS).astype(int)

        new_self[
            self_start[0]: self_start[0] + self.matrix.shape[0],
            self_start[1]: self_start[1] + self.matrix.shape[1],
            self_start[2]: self_start[2] + self.matrix.shape[2],
        ] = self.matrix

        new_other[
            other_start[0]: other_start[0] + other.matrix.shape[0],
            other_start[1]: other_start[1] + other.matrix.shape[1],
            other_start[2]: other_start[2] + other.matrix.shape[2],
        ] = other.matrix

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
            min_voxel_coords[0]: max_voxel_coords[0] + 1,
            min_voxel_coords[1]: max_voxel_coords[1] + 1,
            min_voxel_coords[2]: max_voxel_coords[2] + 1,
        ]

        # Calculate new matrix_origin_center
        new_origin_center = np.round(self.min_grid_center + min_voxel_coords * self.voxel_size, DECIMALS)

        return self.__class__(cropped_matrix, tuple(new_origin_center), self.voxel_size, self.name)


class Pixelization(DiscreteRepresentation, DessiaObject):
    """
    Abstract base class for creating and manipulating pixelizations of volmdlr geometries.

    This approach is used to create a pixelization of the contour, without filling the surface.

    The pixelization is defined on an implicit 2D grid, where each pixel is defined by its size.
    The implicit grid consists of axis-aligned squares with a given size, ensuring that the global
    origin (0, 0) is always a corner point of a pixel.

    For example, in 1D with a pixel size of `t`, the set of pixels (defined by minimum and maximum points)
    is: {i ∈ N, t ∈ R | (i * t, (i+1) * t)}
    The corresponding set of pixel centers is: {i ∈ N, t ∈ R | (i + 0.5) * t}

    This approach enables consistent pixelization across the 3D space, facilitating fast Boolean operations.
    """

    PixelizationType = TypeVar("PixelizationType", bound="Pixelization")

    def __init__(self, pixel_size: float, name: str):
        """
        Initialize the pixelization.

        :param pixel_size: The pixel edges size.
        :type pixel_size: float
        :param name: The name of the pixelization.
        :type name: str, optional
        """
        DiscreteRepresentation.__init__(self, element_size=pixel_size)
        DessiaObject.__init__(self, name=name)

    @property
    def pixel_size(self) -> float:
        """
        Get the pixel size.

        :return: The pixel size.
        :rtype: float
        """
        return self.element_size

    def get_pixel_centers(self) -> Set[_Point2D]:
        """
        Get the center point of each pixel.

        :return: The center point of each pixel.
        :rtype: set[tuple[float, float]]
        """
        return self._get_element_centers()

    @property
    def area(self) -> float:
        """
        Calculate the area of the pixelization.

        :return: The volume of the pixelization.
        :rtype: float
        """
        return len(self) * self.pixel_size**2

    @property
    def bounding_rectangle(self) -> BoundingRectangle:
        """
        Get the bounding rectangle of the pixelization.

        :return: The bounding rectangle of the pixelization.
        :rtype: BoundingRectangle
        """
        min_point = np.round((np.array([self.min_grid_center]) - self.pixel_size)[0], DECIMALS)
        max_point = np.round((np.array([self.max_grid_center]) + self.pixel_size)[0], DECIMALS)

        return BoundingRectangle(min_point[0], max_point[0], min_point[1], max_point[1])

    # CLASS METHODS
    @classmethod
    @abstractmethod
    def from_line_segment(cls, line_segment: LineSegment2D, pixel_size: float, name: str = "") -> PixelizationType:
        """
        Create a pixelization from a LineSegment2D.

        :param line_segment: The LineSegment2D to create the pixelization from.
        :type line_segment: LineSegment2D
        :param pixel_size: The size of each pixel.
        :type pixel_size: float
        :param name: Optional name for the pixelization.
        :type name: str

        :return: A pixelization created from the LineSegment2D.
        :rtype: PixelizationType
        """

    @classmethod
    @abstractmethod
    def from_closed_polygon(
        cls, closed_polygon: ClosedPolygon2D, pixel_size: float, name: str = ""
    ) -> PixelizationType:
        """
        Create a pixelization from a ClosedPolygon2D.

        :param closed_polygon: The ClosedPolygon2D to create the pixelization from.
        :type closed_polygon: ClosedPolygon2D
        :param pixel_size: The size of each pixel.
        :type pixel_size: float
        :param name: Optional name for the pixelization.
        :type name: str

        :return: A pixelization created from the ClosedPolygon2D.
        :rtype: PixelizationType
        """

    # FILLING METHODS
    def fill_outer_pixels(self) -> PixelizationType:
        """
        Fill the outer pixels of the pixelization.

        :return: A new pixelization with outer pixels filled.
        :rtype: PixelizationType
        """
        return self._fill_outer_elements()

    def fill_enclosed_pixels(self) -> PixelizationType:
        """
        Fill the enclosed pixels of the pixelization.

        :return: A new pixelization with enclosed pixels filled.
        :rtype: PixelizationType
        """
        return self._fill_enclosed_elements()

    # DISPLAY METHOD
    def plot(self, ax=None, color: str = "b", **kwargs):
        """
        Plots the pixels on a 2D plane.

        :param ax: An existing figure to plot on.
        :param color: The color of the pixels.
        :param kwargs: Additional keyword arguments.

        :return: The plotted MatPlotLib figure.
        """
        # pylint: disable=arguments-renamed

        if ax is None:
            _, ax = plt.subplots()

        for center in self.get_pixel_centers():
            x, y = center
            ax.add_patch(
                patches.Rectangle(
                    (x - self.pixel_size / 2, y - self.pixel_size / 2),  # Bottom left corner
                    self.pixel_size,  # Width
                    self.pixel_size,  # Height
                    color=color,
                )
            )

        # Setting the x and y limits to contain all the pixels
        ax.set_xlim(
            min(x[0] for x in self.get_pixel_centers()) - self.pixel_size,
            max(x[0] for x in self.get_pixel_centers()) + self.pixel_size,
        )
        ax.set_ylim(
            min(x[1] for x in self.get_pixel_centers()) - self.pixel_size,
            max(x[1] for x in self.get_pixel_centers()) + self.pixel_size,
        )

        ax.set_aspect("equal")  # Ensuring equal scaling for both axes
        plt.show()

        return ax

    # HELPER METHODS
    @staticmethod
    def _extract_segment_from_line_segment(
        line_segment: LineSegment2D,
    ) -> _Segment2D:
        """
        Extract the segment coordinates from the given LineSegment2D.

        :param line_segment: The LineSegment2D to extract from.
        :type line_segment: LineSegment2D

        :return: The extracted segment coordinates.
        :rtype: tuple[tuple[float, float], tuple[float, float]]
        """
        return (line_segment.start.x, line_segment.start.y), (line_segment.end.x, line_segment.end.y)

    @staticmethod
    def _extract_segments_from_closed_polygon(
        closed_polygon: ClosedPolygon2D,
    ) -> List[_Segment2D]:
        """
        Extract the segments coordinates from the given ClosedPolygon2D.

        :param closed_polygon: The ClosedPolygon2D to extract from.
        :type closed_polygon: ClosedPolygon2D

        :return: The extracted segments coordinates.
        :rtype: list[tuple[tuple[float, float], tuple[float, float]]]
        """
        return [
            (
                (closed_polygon.points[i - 1].x, closed_polygon.points[i - 1].y),
                (closed_polygon.points[i].x, closed_polygon.points[i].y),
            )
            for i in range(len(closed_polygon.points))
        ]


class PointBasedPixelization(Pixelization):
    """Pixelization implemented as a set of points, representing each pixel center."""

    def __init__(self, pixel_centers: Set[_Point2D], pixel_size: float, name: str = ""):
        """
        Initialize the PointBasedPixelization.

        :param pixel_centers: The set of points representing pixel centers.
        :type pixel_centers: set[tuple[float, float]]
        :param pixel_size: The pixel edges size.
        :type pixel_size: float
        :param name: The name of the pixelization.
        :type name: str, optional
        """
        self._check_element_size_number_of_decimals(pixel_size)

        self.pixel_centers = pixel_centers

        Pixelization.__init__(self, pixel_size=pixel_size, name=name)

    def _get_element_centers(self) -> Set[_Point2D]:
        """
        Get the center point of each pixel.

        :return: The center point of each pixel.
        :rtype: set[tuple[float, float]]
        """
        return self.pixel_centers

    def __eq__(self, other: "PointBasedPixelization") -> bool:
        """
        Check if two pixelizations are equal.

        :param other: Another pixelization to compare with.
        :type other: PointBasedPixelization

        :return: True if the pixelizations are equal, False otherwise.
        :rtype: bool
        """
        return self.pixel_size == other.pixel_size and self.pixel_centers == other.pixel_centers

    def __len__(self):
        """
        Get the number of pixel in the pixelization.

        :return: The number of pixels in the pixelization (i.e. the number of pixel centers).
        :rtype: int
        """
        return len(self.pixel_centers)

    @property
    def min_grid_center(self) -> _Point2D:
        """
        Get the minimum center point from the set of pixel centers, in the pixel 3D grid.

        This point may not be a pixel of the pixelization, because it is the minimum center in each direction (X, Y).

        :return: The minimum center point.
        :rtype: tuple[float, float]
        """
        min_x = min_y = float("inf")

        for point in self.pixel_centers:
            min_x = min(min_x, point[0])
            min_y = min(min_y, point[1])

        return min_x, min_y

    @property
    def max_grid_center(self) -> _Point2D:
        """
        Get the maximum center point from the set of pixel centers, in the pixel 3D grid.

        This point may not be a pixel of the pixelization, because it is the maximum center in each direction (X, Y).

        :return: The maximum center point.
        :rtype: tuple[float, float]
        """
        max_x = max_y = -float("inf")

        for point in self.pixel_centers:
            max_x = max(max_x, point[0])
            max_y = max(max_y, point[1])

        return max_x, max_y

    # CLASS METHODS
    @classmethod
    def from_line_segment(
        cls, line_segment: LineSegment2D, pixel_size: float, name: str = ""
    ) -> "PointBasedPixelization":
        """
        Create a pixelization from a LineSegment2D.

        :param line_segment: The LineSegment2D to create the pixelization from.
        :type line_segment: LineSegment2D
        :param pixel_size: The size of each pixel.
        :type pixel_size: float
        :param name: Optional name for the pixelization.
        :type name: str

        :return: A pixelization created from the LineSegment2D.
        :rtype: PointBasedPixelization
        """

        return cls(
            line_segments_to_pixels([cls._extract_segment_from_line_segment(line_segment)], pixel_size),
            pixel_size,
            name,
        )

    @classmethod
    def from_closed_polygon(
        cls, closed_polygon: ClosedPolygon2D, pixel_size: float, name: str = ""
    ) -> "PointBasedPixelization":
        """
        Create a pixelization from a ClosedPolygon2D.

        :param closed_polygon: The ClosedPolygon2D to create the pixelization from.
        :type closed_polygon: ClosedPolygon2D
        :param pixel_size: The size of each pixel.
        :type pixel_size: float
        :param name: Optional name for the pixelization.
        :type name: str

        :return: A pixelization created from the ClosedPolygon2D.
        :rtype: PointBasedPixelization
        """

        line_segments = cls._extract_segments_from_closed_polygon(closed_polygon)
        return cls(line_segments_to_pixels(line_segments, pixel_size), pixel_size, name)

    @classmethod
    def from_matrix_based_pixelization(
        cls, matrix_based_pixelization: "MatrixBasedPixelization"
    ) -> "PointBasedPixelization":
        """
        Create a PointBasedPixelization object from a MatrixBasedPixelization.

        :param matrix_based_pixelization: The MatrixBasedPixelization object representing the pixelization.
        :type matrix_based_pixelization: MatrixBasedPixelization

        :return: A PointBasedPixelization object created from the MatrixBasedPixelization.
        :rtype: PointBasedPixelization
        """
        if not cls.check_center_is_in_implicit_grid(
            matrix_based_pixelization.min_grid_center, matrix_based_pixelization.pixel_size
        ):
            warnings.warn(
                """This matrix based pixelization is not defined in the implicit grid defined by the pixel_size.
            Some methods like boolean operation or interference computing may not work as expected."""
            )

        return matrix_based_pixelization.to_point_based_pixelization()

    # BOOLEAN OPERATIONS
    def union(self, other: "PointBasedPixelization") -> "PointBasedPixelization":
        """
        Perform a union operation with another PointBasedPixelization.

        :param other: The PointBasedPixelization to perform the union with.
        :type other: PointBasedPixelization

        :return: A new PointBasedPixelization resulting from the union operation.
        :rtype: PointBasedPixelization
        """
        if self.pixel_size != other.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform union.")

        return PointBasedPixelization(self.pixel_centers.union(other.pixel_centers), self.pixel_size)

    def difference(self, other: "PointBasedPixelization") -> "PointBasedPixelization":
        """
        Perform an intersection operation with another PointBasedPixelization.

        :param other: The PointBasedPixelization to perform the intersection with.
        :type other: PointBasedPixelization

        :return: A new PointBasedPixelization resulting from the intersection operation.
        :rtype: PointBasedPixelization
        """
        if self.pixel_size != other.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform difference.")

        return PointBasedPixelization(self.pixel_centers.difference(other.pixel_centers), self.pixel_size)

    def intersection(self, other: "PointBasedPixelization") -> "PointBasedPixelization":
        """
        Create a pixelization that is the Boolean intersection of two pixelization.

        Both pixelization must have same pixel size.

        :param other: The other pixelization to compute the Boolean intersection with.
        :type other: PointBasedPixelization

        :return: The created pixelization resulting from the Boolean intersection.
        :rtype: PointBasedPixelization
        """
        if self.pixel_size != other.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform intersection.")

        return PointBasedPixelization(self.pixel_centers.intersection(other.pixel_centers), self.pixel_size)

    def symmetric_difference(self, other: "PointBasedPixelization") -> "PointBasedPixelization":
        """
        Perform a symmetric difference operation with another PointBasedPixelization.

        :param other: The PointBasedPixelization to perform the symmetric difference with.
        :type other: PointBasedPixelization

        :return: A new PointBasedPixelization resulting from the symmetric difference operation.
        :rtype: PointBasedPixelization
        """
        if self.pixel_size != other.pixel_size:
            raise ValueError("Both pixelizations must have the same pixel_size to perform symmetric difference.")

        return PointBasedPixelization(self.pixel_centers.symmetric_difference(other.pixel_centers), self.pixel_size)

    def inverse(self) -> "PointBasedPixelization":
        """
        Compute the inverse of the pixelization.

        :return: A new pixelization representing the inverse.
        :rtype: PointBasedPixelization
        """
        inverted_pixel_matrix = self.to_matrix_based_pixelization().inverse()

        return PointBasedPixelization.from_matrix_based_pixelization(inverted_pixel_matrix)

    # FILLING METHODS
    def flood_fill(self, start: _Point2D, fill_with: bool) -> "PointBasedPixelization":
        """
        Perform a flood fill operation on the pixelization.

        :param start: The coordinates of the starting point for the flood fill.
        :type start: tuple[float, float]
        :param fill_with: The value to fill the pixels with during the operation.
        :type fill_with: bool

        :return: A new pixelization resulting from the flood fill operation.
        :rtype: PointBasedPixelization
        """
        start = self._point_to_local_grid_index(start)
        pixel_matrix = self.to_matrix_based_pixelization()
        filled_pixel_matrix = pixel_matrix.flood_fill(start, fill_with)

        return self.from_matrix_based_pixelization(filled_pixel_matrix)

    def _fill_outer_elements(self) -> "PointBasedPixelization":
        """
        Fill the outer pixels of the pixelization.

        :return: A new pixelization with outer pixels filled.
        :rtype: PointBasedPixelization
        """
        return self.from_matrix_based_pixelization(self.to_matrix_based_pixelization().fill_outer_pixels())

    def _fill_enclosed_elements(self) -> "PointBasedPixelization":
        """
        Fill the enclosed pixels of the pixelization.

        :return: A new pixelization with enclosed pixels filled.
        :rtype: PointBasedPixelization
        """
        return self.from_matrix_based_pixelization(self.to_matrix_based_pixelization().fill_enclosed_pixels())

    # SERIALIZATION
    def to_dict(
        self, use_pointers: bool = True, memo=None, path: str = "#", id_method=True, id_memo=None
    ) -> JsonSerializable:
        """Specific 'to_dict' method to allow serialization of a set."""
        dict_ = self.base_dict()

        dict_["pixel_centers"] = list(self.pixel_centers)
        dict_["pixel_size"] = self.pixel_size
        dict_["name"] = self.name

        return dict_

    @classmethod
    def dict_to_object(
        cls,
        dict_: JsonSerializable,
        force_generic: bool = False,
        global_dict=None,
        pointers_memo: Dict[str, Any] = None,
        path: str = "#",
    ) -> "PointBasedPixelization":
        """Specific 'dict_to_object' method to allow deserialization of a set."""

        pixel_centers = set(tuple(pixel_center) for pixel_center in dict_["pixel_centers"])
        pixel_size = dict_["pixel_size"]
        name = dict_["name"]

        return cls(pixel_centers, pixel_size, name)

    # HELPER METHODS
    def to_matrix_based_pixelization(self) -> "MatrixBasedPixelization":
        """
        Convert the point based pixelization to a matrix based pixelization.

        :return: The matrix based pixelization.
        :rtype: MatrixBasedPixelization
        """
        min_center = self.min_grid_center
        max_center = self.max_grid_center

        dim_x = round((max_center[0] - min_center[0]) / self.pixel_size + 1)
        dim_y = round((max_center[1] - min_center[1]) / self.pixel_size + 1)

        indices = np.round((np.array(list(self.pixel_centers)) - min_center) / self.pixel_size).astype(int)

        matrix = np.zeros((dim_x, dim_y), dtype=np.bool_)
        matrix[indices[:, 0], indices[:, 1]] = True

        return MatrixBasedPixelization(matrix, min_center, self.pixel_size, self.name)

    def _point_to_local_grid_index(self, point: _Point2D) -> Tuple[int, int]:
        """
        Convert a point to the local grid index within the pixelization.

        :param point: The point to convert.
        :type point: tuple[float, float]

        :return: The local grid index of the point.
        :rtype: Tuple[int, int]

        :raises ValueError: If the point is not within the bounding rectangle of the pixelization.
        """
        if not self.bounding_rectangle.point_belongs(Point2D(*point)):
            raise ValueError("Point not in local pixel grid.")

        x_index = int((point[0] - self.bounding_rectangle.xmin) // self.pixel_size)
        y_index = int((point[1] - self.bounding_rectangle.ymin) // self.pixel_size)

        return x_index, y_index


class MatrixBasedPixelization(Pixelization):
    """Pixelization implemented as a 2D matrix."""

    def __init__(
        self,
        pixel_matrix: NDArray[np.bool_],
        min_grid_center: _Point2D,
        pixel_size: float,
        name: str = "",
    ):
        """
        Initialize a MatrixBasedPixelization.

        :param pixel_matrix: The pixel numpy matrix object representing the pixelization.
        :type pixel_matrix: np.ndarray[np.bool_, np.ndim == 3]
        :param pixel_size: The size of the pixel edges.
        :type pixel_size: float
        :param min_grid_center: Minimum pixel center point of the pixel grid matrix, i.e 'matrix[0][0]'.
        This point may not be a pixel of the pixelization, because it's the minimum center in each direction (X, Y).
        :type min_grid_center: tuple[float, float]
        """
        self._check_element_size_number_of_decimals(pixel_size)

        self.matrix = pixel_matrix
        self._min_grid_center = min_grid_center

        Pixelization.__init__(self, pixel_size=pixel_size, name=name)

    def _get_element_centers(self) -> Set[_Point2D]:
        """
        Get the center point of each pixel.

        :return: The center point of each pixel.
        :rtype: set[tuple[float, float]]
        """
        indices = np.argwhere(self.matrix)
        pixel_centers = self.min_grid_center + indices * self.pixel_size

        return set(map(tuple, np.round(pixel_centers, DECIMALS)))

    def __eq__(self, other: "MatrixBasedPixelization") -> bool:
        """
        Check if two MatrixBasedPixelization are equal.

        :param other: Another MatrixBasedPixelization to compare with.
        :type other: MatrixBasedPixelization

        :return: True if the MatrixBasedPixelization are equal, False otherwise.
        :rtype: bool
        """
        return (
            self.pixel_size == other.pixel_size
            and self.min_grid_center == other.min_grid_center
            and np.array_equal(self.matrix, other.matrix)
        )

    def __len__(self) -> int:
        """
        Get the number of pixels in the pixelization (i.e. the number of True value in the 2D pixel matrix).

        :return: The number of pixels in the pixelization.
        :rtype: int
        """
        return len(np.argwhere(self.matrix))

    @property
    def min_grid_center(self) -> _Point2D:
        """
        Get the minimum center point from the set of pixel centers, in the pixel 3D grid.

        This point may not be a pixel of the pixelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, float]
        """
        return self._min_grid_center

    @property
    def max_grid_center(self) -> _Point2D:
        """
        Get the maximum center point from the set of pixel centers, in the pixel 3D grid.

        This point may not be a pixel of the pixelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, float]
        """
        return tuple(
            np.round(
                np.array(self.min_grid_center) + (np.array(self.matrix.shape) - 1) * self.pixel_size,
                6,
            )
        )

    # CLASS METHODS
    @classmethod
    def from_line_segment(
        cls, line_segment: LineSegment2D, pixel_size: float, name: str = ""
    ) -> "MatrixBasedPixelization":
        """
        Create a pixelization from a LineSegment2D.

        :param line_segment: The LineSegment2D to create the pixelization from.
        :type line_segment: LineSegment2D
        :param pixel_size: The size of each pixel.
        :type pixel_size: float
        :param name: Optional name for the pixelization.
        :type name: str

        :return: A pixelization created from the LineSegment2D.
        :rtype: MatrixBasedPixelization
        """

        return PointBasedPixelization.from_line_segment(line_segment, pixel_size, name).to_matrix_based_pixelization()

    @classmethod
    def from_closed_polygon(
        cls, closed_polygon: ClosedPolygon2D, pixel_size: float, name: str = ""
    ) -> "MatrixBasedPixelization":
        """
        Create a pixelization from a ClosedPolygon2D.

        :param closed_polygon: The ClosedPolygon2D to create the pixelization from.
        :type closed_polygon: ClosedPolygon2D
        :param pixel_size: The size of each pixel.
        :type pixel_size: float
        :param name: Optional name for the pixelization.
        :type name: str

        :return: A pixelization created from the ClosedPolygon2D.
        :rtype: MatrixBasedPixelization
        """

        return PointBasedPixelization.from_closed_polygon(
            closed_polygon, pixel_size, name
        ).to_matrix_based_pixelization()

    @classmethod
    def from_point_based_pixelization(
        cls, point_based_pixelization: "PointBasedPixelization"
    ) -> "MatrixBasedPixelization":
        """
        Create a MatrixBasedPixelization from a PointBasedPixelization.

        :param point_based_pixelization: The PointBasedPixelization to create the MatrixBasedPixelization from.
        :type point_based_pixelization: PointBasedPixelization

        :return: A MatrixBasedPixelization created from the PointBasedPixelization.
        :rtype: MatrixBasedPixelization
        """
        return point_based_pixelization.to_matrix_based_pixelization()

    # BOOLEAN OPERATIONS
    def union(self, other: "MatrixBasedPixelization") -> "MatrixBasedPixelization":
        """
        Perform a union operation with another MatrixBasedPixelization.

        :param other: The MatrixBasedPixelization to perform the union with.
        :type other: MatrixBasedPixelization

        :return: A new MatrixBasedPixelization resulting from the union operation.
        :rtype: MatrixBasedPixelization
        """
        return self._logical_operation(other, np.logical_or)

    def difference(self, other: "MatrixBasedPixelization") -> "MatrixBasedPixelization":
        """
        Perform a difference operation with another pixelization.

        :param other: The MatrixBasedPixelization to perform the difference with.
        :type other: MatrixBasedPixelization

        :return: A new MatrixBasedPixelization resulting from the difference operation.
        :rtype: MatrixBasedPixelization
        """
        return self._logical_operation(other, lambda a, b: np.logical_and(a, np.logical_not(b)))

    def intersection(self, other: "MatrixBasedPixelization") -> "MatrixBasedPixelization":
        """
        Perform an intersection operation with another MatrixBasedPixelization.

        :param other: The MatrixBasedPixelization to perform the intersection with.
        :type other: MatrixBasedPixelization

        :return: A new MatrixBasedPixelization resulting from the intersection operation.
        :rtype: MatrixBasedPixelization
        """
        return self._logical_operation(other, np.logical_and)

    def symmetric_difference(self, other: "MatrixBasedPixelization") -> "MatrixBasedPixelization":
        """
        Perform a symmetric difference operation with another MatrixBasedPixelization.

        :param other: The MatrixBasedPixelization to perform the symmetric difference with.
        :type other: MatrixBasedPixelization

        :return: A new MatrixBasedPixelization resulting from the symmetric difference operation.
        :rtype: MatrixBasedPixelization
        """
        return self._logical_operation(other, np.logical_xor)

    def inverse(self) -> "MatrixBasedPixelization":
        """
        Compute the inverse of the pixelization.

        :return: A new pixelization representing the inverse.
        :rtype: MatrixBasedPixelization
        """
        inverted_matrix = np.logical_not(self.matrix)
        return self.__class__(inverted_matrix, self.min_grid_center, self.pixel_size)

    def flood_fill(self, start: Tuple[int, int], fill_with: bool) -> "MatrixBasedPixelization":
        """
        Perform a flood fill operation on the pixelization.

        :param start: The indexes of the starting pixel in the 3D matrix for the flood fill.
        :type start: tuple[int, int]
        :param fill_with: The value to fill the pixels with during the operation.
        :type fill_with: bool

        :return: A new pixelization resulting from the flood fill operation.
        :rtype: MatrixBasedPixelization
        """
        return self.__class__(
            flood_fill_matrix_2d(self.matrix, start, fill_with), self.min_grid_center, self.pixel_size
        )

    def _fill_outer_elements(self) -> "MatrixBasedPixelization":
        """
        Fill the outer pixels of the pixelization.

        :return: A new pixelization with outer pixels filled.
        :rtype: MatrixBasedPixelization
        """
        # pylint: disable=protected-access

        expanded_pixel_matrix = self._expand()
        outer_filled_expanded_pixel_matrix = expanded_pixel_matrix.flood_fill((0, 0), True)
        outer_filled_pixel_matrix = outer_filled_expanded_pixel_matrix._reduce()

        return outer_filled_pixel_matrix

    def _fill_enclosed_elements(self) -> "MatrixBasedPixelization":
        """
        Fill the enclosed pixels of the pixelization.

        :return: A new pixelization with enclosed pixels filled.
        :rtype: MatrixBasedPixelization
        """
        outer_filled_pixel_matrix = self.fill_outer_pixels()
        inner_filled_pixel_matrix = self + outer_filled_pixel_matrix.inverse()

        return inner_filled_pixel_matrix

    # SERIALIZATION
    def to_dict(
        self, use_pointers: bool = True, memo=None, path: str = "#", id_method=True, id_memo=None
    ) -> JsonSerializable:
        """Specific 'to_dict' method to allow serialization of a numpy array."""
        dict_ = self.base_dict()

        dict_["matrix"] = self.matrix.tolist()
        dict_["min_grid_center"] = self.min_grid_center
        dict_["pixel_size"] = self.pixel_size
        dict_["name"] = self.name

        return dict_

    @classmethod
    def dict_to_object(
        cls,
        dict_: JsonSerializable,
        force_generic: bool = False,
        global_dict=None,
        pointers_memo: Dict[str, Any] = None,
        path: str = "#",
    ) -> "MatrixBasedPixelization":
        """Specific 'dict_to_object' method to allow deserialization of a numpy array."""

        matrix = np.array(dict_["matrix"])
        min_grid_center = tuple(dict_["min_grid_center"])
        pixel_size = dict_["pixel_size"]
        name = dict_["name"]

        return cls(matrix, min_grid_center, pixel_size, name)

    # HELPER METHODS
    def to_point_based_pixelization(self) -> "PointBasedPixelization":
        """
        Convert the MatrixBasedPixelization to a PointBasedPixelization.

        :return: A PointBasedPixelization representation of the current pixelization.
        :rtype: PointBasedPixelization
        """
        return PointBasedPixelization(self.get_pixel_centers(), self.pixel_size, self.name)

    def _expand(self) -> "MatrixBasedPixelization":
        """
        Expand the pixelization matrix by adding a single layer of False pixels around the existing matrix.

        :return: A new MatrixBasedPixelization with an expanded pixelization matrix.
        :rtype: MatrixBasedPixelization
        """
        current_shape = self.matrix.shape
        new_shape = tuple(dim + 2 for dim in current_shape)
        expanded_matrix = np.zeros(new_shape, dtype="bool")
        slices = tuple(slice(1, -1) for _ in current_shape)
        expanded_matrix[slices] = self.matrix.copy()

        return self.__class__(
            expanded_matrix,
            tuple(np.round(np.array(self.min_grid_center) - self.pixel_size, DECIMALS)),
            self.pixel_size,
        )

    def _reduce(self) -> "MatrixBasedPixelization":
        """
        Reduce the size of the pixelization matrix by removing a single layer of pixels around the edges.

        :return: A new MatrixBasedPixelization with a reduced pixelization matrix.
        :rtype: MatrixBasedPixelization
        """
        current_shape = self.matrix.shape
        slices = tuple(slice(1, -1) for _ in current_shape)
        reduced_matrix = self.matrix.copy()[slices]

        return self.__class__(
            reduced_matrix,
            tuple(np.round(np.array(self.min_grid_center) + self.pixel_size, DECIMALS)),
            self.pixel_size,
        )

    def _logical_operation(self, other: "MatrixBasedPixelization", logical_operation) -> "MatrixBasedPixelization":
        """
        Perform a logical operation (e.g., union, intersection, etc.) between two pixelizations.

        :param other: The other MatrixBasedPixelization to perform the operation with.
        :type other: MatrixBasedPixelization
        :param logical_operation: The logical operation function to apply (e.g., np.logical_or).

        :return: A new MatrixBasedPixelization resulting from the logical operation.
        :rtype: MatrixBasedPixelization
        """
        if self.pixel_size != other.pixel_size:
            raise ValueError("Pixel sizes must be the same to perform boolean operations.")

        self_min, self_max = np.array(self.min_grid_center), np.array(self.max_grid_center) + 1
        other_min, other_max = (
            np.array(other.min_grid_center),
            np.array(other.max_grid_center) + 1,
        )

        global_min = np.min([self_min, other_min], axis=0)
        global_max = np.max([self_max, other_max], axis=0)

        new_shape = np.round((global_max - global_min) / self.pixel_size, DECIMALS).astype(int) + 2

        new_self = np.zeros(new_shape, dtype=np.bool_)
        new_other = np.zeros(new_shape, dtype=np.bool_)

        self_start = np.round((self_min - global_min) / self.pixel_size, DECIMALS).astype(int)
        other_start = np.round((other_min - global_min) / self.pixel_size, DECIMALS).astype(int)

        new_self[
            self_start[0]: self_start[0] + self.matrix.shape[0],
            self_start[1]: self_start[1] + self.matrix.shape[1],
        ] = self.matrix

        new_other[
            other_start[0]: other_start[0] + other.matrix.shape[0],
            other_start[1]: other_start[1] + other.matrix.shape[1],
        ] = other.matrix

        result_matrix = logical_operation(new_self, new_other)

        return self.__class__(result_matrix, tuple(global_min), self.pixel_size).crop_matrix()

    def crop_matrix(self) -> "MatrixBasedPixelization":
        """
        Crop the pixel matrix to the smallest possible size.

        :return: The MatrixBasedPixelization with cropped pixel matrix.
        :rtype: MatrixBasedPixelization
        """
        # Find the indices of the True pixels
        true_pixels = np.argwhere(self.matrix)

        if len(true_pixels) == 0:
            # Can't crop a matrix of False values
            return self

        # Find the minimum and maximum indices along each axis
        min_pixel_coords, max_pixel_coords = np.round(np.min(true_pixels, axis=0), DECIMALS), np.round(
            np.max(true_pixels, axis=0), 6
        )

        # Crop the matrix to the smallest possible size
        cropped_matrix = self.matrix[
            min_pixel_coords[0]: max_pixel_coords[0] + 1,
            min_pixel_coords[1]: max_pixel_coords[1] + 1,
        ]

        # Calculate new matrix_origin_center
        new_origin_center = np.round(self.min_grid_center + min_pixel_coords * self.pixel_size, DECIMALS)

        return self.__class__(cropped_matrix, tuple(new_origin_center), self.pixel_size, self.name)
