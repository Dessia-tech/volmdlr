"""
Class for discrete representations of volmdlr models (voxelization for 3D geometries, pixelization for 2D geometries).
"""
import itertools
import math
import warnings
from typing import Any, Dict, Iterable, List, Set, Tuple, TypeVar, Union

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
    round_point_3d_to_digits,
    round_to_digits,
    triangle_intersects_voxel,
    triangles_to_voxel_matrix,
    voxel_triangular_faces,
)
from volmdlr.edges import LineSegment2D
from volmdlr.faces import Face3D, Triangle3D
from volmdlr.shells import ClosedTriangleShell3D, DisplayTriangleShell3D, Shell3D
from volmdlr.wires import ClosedPolygon2D

# pylint: disable=no-name-in-module,too-many-lines


# CUSTOM TYPES
_Point3D = Tuple[float, float, float]
_Triangle3D = Tuple[_Point3D, _Point3D, _Point3D]

_Point2D = Tuple[float, float]
_Segment2D = Tuple[_Point2D, _Point2D]

Octree = Union[int, List["Octree"]]


# GLOBAL VARIABLE
DECIMALS = 9  # Used to round numbers and avoid floating point arithmetic imprecision


class DiscreteRepresentation:
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

    def __eq__(self, other: DiscreteRepresentationType) -> bool:
        """
        Check if two discrete representations are equal.

        :param other: Another discrete representation to compare with.
        :type other: DiscreteRepresentationType

        :return: True if the discrete representations are equal, False otherwise.
        :rtype: bool
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def __len__(self) -> int:
        """
        Get the number of elements in the discrete representation.

        :return: The number of elements in the discrete representation.
        :rtype: int
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def _get_element_centers(self) -> Set[Tuple]:
        """
        Get the center point of each element.

        :return: The center point of each element.
        :rtype: set[tuple[float, ...]]
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    @property
    def min_grid_center(self) -> Tuple:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, ...]
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    @property
    def max_grid_center(self) -> Tuple:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, ...]
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

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

    def union(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform a union operation with another discrete representation.

        :param other: The discrete representation to perform the union with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the union operation.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def difference(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform a difference operation with another discrete representation.

        :param other: The discrete representation to perform the difference with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the difference operation.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def intersection(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform an intersection operation with another discrete representation.

        :param other: The discrete representation to perform the intersection with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the intersection operation.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def symmetric_difference(self, other: DiscreteRepresentationType) -> DiscreteRepresentationType:
        """
        Perform a symmetric difference operation with another discrete representation.

        :param other: The discrete representation to perform the symmetric difference with.
        :type other: DiscreteRepresentationType

        :return: A new discrete representation resulting from the symmetric difference operation.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def inverse(self) -> DiscreteRepresentationType:
        """
        Compute the inverse of the discrete representation.

        :return: A new discrete representation representing the inverse.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

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
    def flood_fill(self, start, fill_with: bool) -> DiscreteRepresentationType:
        """
        Perform a flood fill operation on the discrete representation.

        :param start: The starting point for the flood fill.
        :param fill_with: The value to fill the elements with during the operation.
        :type fill_with: bool

        :return: A new discrete representation resulting from the flood fill operation.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def _fill_outer_elements(self) -> DiscreteRepresentationType:
        """
        Fill the outer elements of the discrete representation.

        :return: A new discrete representation with outer elements filled.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

    def _fill_enclosed_elements(self) -> DiscreteRepresentationType:
        """
        Fill the enclosed elements of the discrete representation.

        :return: A new discrete representation with enclosed elements filled.
        :rtype: DiscreteRepresentationType
        """
        raise NotImplementedError("DiscreteRepresentation is an abstract class and should not be use directly.")

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
            if not round_to_digits((coord - 0.5 * element_size) / element_size, DECIMALS).is_integer():
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
            decimals = len(str(element_size + 1).split(".")[1])

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
    def from_triangles(cls, triangles: List[_Triangle3D], voxel_size: float, name: str = "") -> VoxelizationType:
        """
        Create a voxelization from a list of triangles.

        :param triangles: The list of triangles to create the voxelization from.
        :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the list of triangles.
        :rtype: VoxelizationType
        """
        raise NotImplementedError("Voxelization is an abstract class and should not be use directly.")

    @classmethod
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
        if isinstance(shell, DisplayTriangleShell3D):
            return cls.from_mesh_data(shell.positions, shell.indices, voxel_size, name)

        return cls.from_triangles(cls._shell_to_triangles(shell), voxel_size, name)

    @classmethod
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
        return cls.from_triangles(cls._volume_model_to_triangles(volume_model), voxel_size, name)

    @classmethod
    def from_mesh_data(
        cls, vertices: Iterable[Iterable[float]], faces: Iterable[Iterable[int]], voxel_size: float, name: str = ""
    ) -> "VoxelizationType":
        """
        Create a voxelization from mesh data.

        :param vertices: The vertices of the mesh.
        :type vertices: Iterable[Iterable[float]]
        :param faces: The faces of the mesh, using vertices indexes.
        :type faces: Iterable[Iterable[int]]
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the voxelization.
        :type name: str

        :return: A voxelization created from the mesh data.
        :rtype: VoxelizationType
        """
        return cls.from_triangles(cls._mesh_data_to_triangles(vertices, faces), voxel_size, name)

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

    def to_display_triangle_shell(self) -> DisplayTriangleShell3D:
        """
        Generate a closed triangle shell representing the voxelization.

        :return: A closed triangle shell representation of the voxelization.
        :rtype: ClosedTriangleShell3D
        """
        # Flatten and round the vertices array
        faces = self.to_triangles()
        vertices = np.array([(face[i][0], face[i][1], face[i][2]) for face in faces for i in range(3)])

        # Get unique vertices and their indices
        vertices, unique_indices = np.unique(vertices, axis=0, return_inverse=True)

        # Create the triangle indices array using NumPy indexing
        flattened_indices = unique_indices.reshape(-1, 3)
        faces = flattened_indices[: len(faces)]

        return DisplayTriangleShell3D(vertices, faces, name=self.name)

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

    @staticmethod
    def _mesh_data_to_triangles(
        vertices: Iterable[Iterable[float]], faces: Iterable[Iterable[int]]
    ) -> List[_Triangle3D]:
        """
        Helper method to convert mesh data to a list of triangles.

        :param vertices: The vertices of the mesh.
        :type vertices: Iterable[Iterable[float]]
        :param faces: The faces of the mesh, using vertices indexes.
        :type faces: Iterable[Iterable[int]]

        :return: The list of triangles extracted from the triangulated primitives of the VolumeModel.
        :rtype: List[Triangle]
        """
        triangles = []

        points = list(vertices)

        for i1, i2, i3 in faces:
            triangles.append((tuple(points[i1]), tuple(points[i2]), tuple(points[i3])))

        return triangles

    @staticmethod
    def voxel_to_bounding_box(voxel_center: _Point3D, voxel_size: float) -> BoundingBox:
        """
        Creates a bounding box from a voxel.

        :param voxel_center: The center point of the voxel.
        :type voxel_center: tuple[float, float, float]
        :param voxel_size: The size of the voxel edge.
        :type voxel_size: float

        :return: The created bounding box.
        :rtype: BoundingBox
        """
        half_size = round_to_digits(voxel_size / 2, DECIMALS)
        min_point = (voxel_center[0] - half_size, voxel_center[1] - half_size, voxel_center[2] - half_size)
        max_point = (voxel_center[0] + half_size, voxel_center[1] + half_size, voxel_center[2] + half_size)

        return BoundingBox(min_point[0], max_point[0], min_point[1], max_point[1], min_point[2], max_point[2])


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
    def from_triangles(
        cls, triangles: List[_Triangle3D], voxel_size: float, name: str = ""
    ) -> "PointBasedVoxelization":
        """
        Create a PointBasedVoxelization from a list of triangles.

        :param triangles: The list of triangles to create the PointBasedVoxelization from.
        :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the PointBasedVoxelization.
        :type name: str

        :return: A PointBasedVoxelization created from the list of triangles.
        :rtype: PointBasedVoxelization
        """
        return cls(MatrixBasedVoxelization.from_triangles(triangles, voxel_size).get_voxel_centers(), voxel_size, name)

    @classmethod
    def from_matrix_based_voxelization(
        cls,
        matrix_based_voxelization: "MatrixBasedVoxelization",
    ) -> "PointBasedVoxelization":
        """
        Create a PointBasedVoxelization object from a MatrixBasedVoxelization.

        :param matrix_based_voxelization: The MatrixBasedVoxelization object representing the voxelization.
        :type matrix_based_voxelization: MatrixBasedVoxelization

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

    # EXPORT METHOD
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

    def to_octree_based_voxelization(self) -> "OctreeBasedVoxelization":
        """
        Convert the PointBasedVoxelization to an OctreeBasedVoxelization.

        :return: The octree based voxelization.
        :rtype: OctreeBasedVoxelization
        """
        return OctreeBasedVoxelization.from_point_based_voxelization(self)

    def to_inner_growing_voxelizations(self, layers_minimal_thickness: float) -> List["PointBasedVoxelization"]:
        """
        Convert the PointBasedVoxelization to multiple PointBasedVoxelization, with different voxel size.

        The more the voxelization is inside, the more its voxel size become bigger.

        :param layers_minimal_thickness: The minimal thickness of each layer.
        :type layers_minimal_thickness: float

        :return: A list of PointBasedVoxelization representing the inner growing voxelizations.
        :rtype: List[PointBasedVoxelization]
        """
        i = 2
        inner_growing_voxel_centers = self.to_octree_based_voxelization().get_inner_growing_voxel_centers(
            layers_minimal_thickness, self.to_matrix_based_voxelization().layers_thickness_by_voxel_centers()
        )

        while len(inner_growing_voxel_centers.keys()) == i:
            # Still making the voxelization inner growing
            max_voxel_size = max(inner_growing_voxel_centers.keys())
            i += 1

            new_point_based_voxelization = PointBasedVoxelization(
                inner_growing_voxel_centers[max_voxel_size], max_voxel_size
            )
            for voxel_size, voxel_centers in (
                new_point_based_voxelization.to_octree_based_voxelization()
                .get_inner_growing_voxel_centers(
                    layers_minimal_thickness,
                    new_point_based_voxelization.to_matrix_based_voxelization().layers_thickness_by_voxel_centers(),
                )
                .items()
            ):
                inner_growing_voxel_centers[voxel_size] = voxel_centers

        point_based_voxelizations = []

        # Create a list of PointBasedVoxelization using the inner growing voxel centers
        for voxel_size, voxel_centers in inner_growing_voxel_centers.items():
            point_based_voxelizations.append(PointBasedVoxelization(voxel_centers, voxel_size))

        return point_based_voxelizations

    # HELPER METHOD
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

    def _voxel_centers_distances_to_faces(self) -> Dict[_Point3D, float]:
        """
        Helper method to compute the minimal distance between the voxel centers and the surface of the voxelization.

        :return: The distance to surface for each voxel center.
        :rtype: dict[tuple(float, float, float): float]
        """
        # pylint: disable=import-outside-toplevel
        from igl import signed_distance

        points_coords = np.array(list(self.voxel_centers))

        display_triangle_shell = self.to_display_triangle_shell()
        vertices, faces = display_triangle_shell.positions, display_triangle_shell.indices
        distances_array = signed_distance(points_coords, vertices, faces.astype(int), sign_type=3)[0]

        if len(self) == 1:
            distances_array = np.array([distances_array])

        # Creating a dictionary to map each point to its distance
        distances_dict = dict(zip(self.voxel_centers, distances_array))

        return distances_dict


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
        :param name: The name of the voxelization.
        :type name: str, optional
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
    def from_triangles(
        cls, triangles: List[_Triangle3D], voxel_size: float, name: str = ""
    ) -> "MatrixBasedVoxelization":
        """
        Create a MatrixBasedVoxelization from a list of triangles.

        :param triangles: The list of triangles to create the MatrixBasedVoxelization from.
        :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the MatrixBasedVoxelization.
        :type name: str

        :return: A MatrixBasedVoxelization created from the list of triangles.
        :rtype: MatrixBasedVoxelization
        """
        matrix, min_grid_center = triangles_to_voxel_matrix(triangles, voxel_size)

        return cls(matrix, min_grid_center, voxel_size, name).crop_matrix()

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

    # EXPORT METHODS
    def to_point_based_voxelization(self) -> "PointBasedVoxelization":
        """
        Convert the MatrixBasedVoxelization to a PointBasedVoxelization.

        :return: A PointBasedVoxelization representation of the current voxelization.
        :rtype: PointBasedVoxelization
        """
        return PointBasedVoxelization(self.get_voxel_centers(), self.voxel_size, self.name)

    def to_octree_based_voxelization(self) -> "OctreeBasedVoxelization":
        """
        Convert the MatrixBasedVoxelization to an OctreeBasedVoxelization.

        :return: The octree based voxelization.
        :rtype: OctreeBasedVoxelization
        """
        return OctreeBasedVoxelization.from_point_based_voxelization(self.to_point_based_voxelization())

    def to_inner_growing_voxelizations(self, layers_minimal_thickness: float) -> List["MatrixBasedVoxelization"]:
        """
        Convert the MatrixBasedVoxelization to multiple MatrixBasedVoxelization, with different voxel size.

        The more the voxelization is inside, the more its voxel size become bigger.

        :param layers_minimal_thickness: The minimal thickness of each layer.
        :type layers_minimal_thickness: float

        :return: A list of MatrixBasedVoxelization representing the inner growing voxelizations.
        :rtype: List[MatrixBasedVoxelization]
        """
        i = 2
        inner_growing_voxel_centers = self.to_octree_based_voxelization().get_inner_growing_voxel_centers(
            layers_minimal_thickness, self.layers_thickness_by_voxel_centers()
        )

        while len(inner_growing_voxel_centers.keys()) == i:
            # Still making the voxelization inner growing
            max_voxel_size = max(inner_growing_voxel_centers.keys())
            i += 1

            new_point_based_voxelization = PointBasedVoxelization(
                inner_growing_voxel_centers[max_voxel_size], max_voxel_size
            )
            for voxel_size, voxel_centers in (
                new_point_based_voxelization.to_octree_based_voxelization()
                .get_inner_growing_voxel_centers(
                    layers_minimal_thickness,
                    new_point_based_voxelization.to_matrix_based_voxelization().layers_thickness_by_voxel_centers(),
                )
                .items()
            ):
                inner_growing_voxel_centers[voxel_size] = voxel_centers

        matrix_based_voxelizations = []

        # Create a list of MatrixBasedVoxelization using the inner growing voxel centers
        for voxel_size, voxel_centers in inner_growing_voxel_centers.items():
            matrix_based_voxelizations.append(
                PointBasedVoxelization(voxel_centers, voxel_size).to_matrix_based_voxelization()
            )

        return matrix_based_voxelizations

    # HELPER METHODS
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

        self_min, self_max = np.array(self.min_grid_center), np.array(self.max_grid_center) + 1
        other_min, other_max = (
            np.array(other.min_grid_center),
            np.array(other.max_grid_center) + 1,
        )

        global_min = np.min([self_min, other_min], axis=0)
        global_max = np.max([self_max, other_max], axis=0)

        new_shape = np.round((global_max - global_min) / self.voxel_size, DECIMALS).astype(int)

        new_self = np.zeros(new_shape, dtype=np.bool_)
        new_other = np.zeros(new_shape, dtype=np.bool_)

        self_start = np.round((self_min - global_min) / self.voxel_size, DECIMALS).astype(int)
        other_start = np.round((other_min - global_min) / self.voxel_size, DECIMALS).astype(int)

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
        new_origin_center = np.round(self.min_grid_center + min_voxel_coords * self.voxel_size, DECIMALS)

        return self.__class__(cropped_matrix, tuple(new_origin_center), self.voxel_size, self.name)

    def layers_thickness_by_voxel_centers(self) -> Dict[_Point3D, float]:
        """
        Get a dictionary with voxel centers as keys and their layers thickness to false or edge values as values.

        :return: Dictionary with voxel centers and their corresponding layers thickness to false or edge values.
        :rtype: dict[tuple[float, float, float], float]
        """

        # Get layer values using layers_to_false_or_edge function
        layer_values = self._layers_matrix()

        # Get the indices of True voxels
        indices = np.argwhere(self.matrix)

        # Calculate voxel centers based on indices
        voxel_centers = self.min_grid_center + indices * self.voxel_size
        rounded_voxel_centers = map(tuple, np.round(voxel_centers, DECIMALS))

        # Fetch the layer value for each voxel center using indices and create the dictionary
        layers_dict = {
            center: round_to_digits(layer_values[tuple(ind)] * self.voxel_size, DECIMALS)
            for center, ind in zip(rounded_voxel_centers, indices)
        }

        return layers_dict

    def _layers_matrix(self) -> NDArray[int]:
        """
        Compute for each voxel of the matrix the number of layer there is until an edge or a False value is reached.

        The layers are cubic, i.e. for one voxel, we check the 26 surrounding voxels.

        :return: The computer layer matrix.
        :rtype: np.ndarray[float, np.ndim == 3]
        """
        result = np.zeros_like(self.matrix, dtype=int)

        for x in range(self.matrix.shape[0]):
            for y in range(self.matrix.shape[1]):
                for z in range(self.matrix.shape[2]):
                    if not self.matrix[x][y][z]:
                        result[x, y, z] = 0
                        continue

                    count = 1
                    size = 1

                    while (
                        0 <= x - count
                        and x + count < self.matrix.shape[0]
                        and 0 <= y - count
                        and y + count < self.matrix.shape[1]
                        and 0 <= z - count
                        and z + count < self.matrix.shape[2]
                    ):
                        cube_slice = self.matrix[
                            slice(x - count, x + count + 1),
                            slice(y - count, y + count + 1),
                            slice(z - count, z + count + 1),
                        ]

                        # If any False is found in the current cube, break the loop
                        if not np.all(cube_slice):
                            break

                        count += 1
                        size += 2

                    result[x, y, z] = count - 1

        return result


class OctreeBasedVoxelization(Voxelization):
    """Voxelization implemented as an octree."""

    # pylint: disable=protected-access,too-many-arguments,too-many-locals,too-many-nested-blocks,too-many-branches

    def __init__(
        self,
        octree: Octree,
        root_center: _Point3D,
        octree_depth: int,
        voxel_size: float,
        triangles: List[_Triangle3D] = None,
        name: str = "",
    ):
        """
        Initialize an OctreeBasedVoxelization.

        :param octree: The octree graph represented using lists.
        :type octree: list[list...list[int]]
        :param root_center: The position of the octree root center.
        :type root_center: tuple[float, float, float]
        :param octree_depth: The depth of the octree.
        :type octree_depth: int
        :param voxel_size: The size of the voxel edge.
        :type voxel_size: float
        :param triangles: The list of triangles used to create the voxelization.
        :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        :param name: The name of the voxelization.
        :type name: str, optional
        """
        self._check_element_size_number_of_decimals(voxel_size)

        self._octree = octree
        self._root_center = root_center
        self._octree_depth = octree_depth

        if triangles:
            self._triangles = triangles
        else:
            self._triangles = []

        Voxelization.__init__(self, voxel_size=voxel_size, name=name)

    def _get_element_centers(self) -> Set[_Point3D]:
        """
        Get the center point of each voxel.

        :return: The center point of each voxel.
        :rtype: set[tuple[float, float, float]]
        """
        return self._get_homogeneous_leaf_centers(0, self._root_voxel_size, self._root_center, self._octree)

    @property
    def _root_voxel_size(self) -> float:
        """
        Get the edge size of the root voxel.

        :return: The edge size of the root voxel.
        :rtype: float
        """
        return round_to_digits(self.voxel_size * 2**self._octree_depth, DECIMALS)

    def __eq__(self, other: "OctreeBasedVoxelization") -> bool:
        """
        Check if two OctreeBasedVoxelization are equal.

        :param other: Another OctreeBasedVoxelization to compare with.
        :type other: OctreeBasedVoxelization

        :return: True if the OctreeBasedVoxelization are equal, False otherwise.
        :rtype: bool
        """
        return self.get_voxel_centers() == other.get_voxel_centers()

    def __len__(self) -> int:
        """
        Get the number of voxels in the voxelization (i.e. the number of True value in the 3D voxel matrix).

        :return: The number of voxels in the voxelization.
        :rtype: int
        """
        return len(self.get_voxel_centers())

    @property
    def min_grid_center(self) -> _Point3D:
        """
        Get the minimum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the minimum center in each direction (X, Y, Z).

        :return: The minimum center point.
        :rtype: tuple[float, float, float]
        """
        return self.to_point_based_voxelization().min_grid_center

    @property
    def max_grid_center(self) -> _Point3D:
        """
        Get the maximum center point from the set of voxel centers, in the voxel 3D grid.

        This point may not be a voxel of the voxelization, because it is the maximum center in each direction (X, Y, Z).

        :return: The maximum center point.
        :rtype: tuple[float, float, float]
        """
        return self.to_point_based_voxelization().max_grid_center

    # CLASS METHODS
    @classmethod
    def from_triangles(
        cls, triangles: List[_Triangle3D], voxel_size: float, name: str = ""
    ) -> "OctreeBasedVoxelization":
        """
        Create a OctreeBasedVoxelization from a list of triangles.

        :param triangles: The list of triangles to create the OctreeBasedVoxelization from.
        :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        :param voxel_size: The size of each voxel.
        :type voxel_size: float
        :param name: Optional name for the OctreeBasedVoxelization.
        :type name: str

        :return: A OctreeBasedVoxelization created from the list of triangles.
        :rtype: OctreeBasedVoxelization
        """
        triangles_np = np.array(triangles)
        min_corner = np.min(np.min(triangles_np, axis=1), axis=0)
        max_corner = np.max(np.max(triangles_np, axis=1), axis=0)

        # Compute the corners in the implicit grid defined by the voxel size
        min_corner = (np.floor_divide(min_corner, voxel_size) - 2) * voxel_size
        max_corner = (np.floor_divide(max_corner, voxel_size) + 2) * voxel_size

        root_size = round_to_digits(np.max(np.maximum(np.abs(min_corner), np.abs(max_corner))) * 2, DECIMALS)

        # Compute the max depth corresponding the voxel_size
        max_depth = math.ceil(math.log2(root_size // voxel_size))
        center = (0.0, 0.0, 0.0)

        sizes = [round_to_digits(voxel_size * 2**i, DECIMALS) for i in range(max_depth, -1, -1)]
        sizes.append(round_to_digits(voxel_size * 1 / 2, DECIMALS))

        octree = cls._subdivide_from_triangles(triangles, list(range(len(triangles))), center, sizes, 0, max_depth)

        return cls(octree, center, max_depth, voxel_size, triangles, name)

    @classmethod
    def from_point_based_voxelization(
        cls, point_based_voxelization: "PointBasedVoxelization"
    ) -> "OctreeBasedVoxelization":
        """
        Create a OctreeBasedVoxelization from a PointBasedVoxelization.

        :param point_based_voxelization: The PointBasedVoxelization to create the OctreeBasedVoxelization from.
        :type point_based_voxelization: PointBasedVoxelization

        :return: A OctreeBasedVoxelization created from the PointBasedVoxelization.
        :rtype: OctreeBasedVoxelization
        """
        min_corner = np.min(np.array(point_based_voxelization.min_grid_center), axis=0)
        max_corner = np.max(np.array(point_based_voxelization.max_grid_center), axis=0)

        voxel_size = point_based_voxelization.voxel_size

        # Compute the corners in the implicit grid defined by the voxel size
        min_corner = (np.floor_divide(min_corner, voxel_size) - 2) * voxel_size
        max_corner = (np.floor_divide(max_corner, voxel_size) + 2) * voxel_size

        root_size = round_to_digits(np.max(np.maximum(np.abs(min_corner), np.abs(max_corner))) * 2, DECIMALS)

        # Compute the max depth corresponding the voxel_size
        max_depth = math.ceil(math.log2(root_size // voxel_size))
        center = (0.0, 0.0, 0.0)

        sizes = [round_to_digits(voxel_size * 2**i, DECIMALS) for i in range(max_depth, -1, -1)]
        sizes.append(round_to_digits(voxel_size * 1 / 2, DECIMALS))

        octree = cls._subdivide_from_points(
            list(point_based_voxelization.voxel_centers),
            center,
            sizes,
            0,
            max_depth,
        )

        return cls(octree, center, max_depth, voxel_size)

    # BOOLEAN OPERATIONS
    def is_intersecting(self, other: "OctreeBasedVoxelization") -> bool:
        """
        Check is two OctreeBasedVoxelization are intersecting.

        :param other: The other voxelization to check if there is an intersection with.
        :type other: OctreeBasedVoxelization

        :return: True if the voxelizations are intersecting, False otherwise.
        :rtype: bool
        """
        self_sizes = [round_to_digits(self.voxel_size * 2**i, DECIMALS) for i in range(self._octree_depth, -1, -1)]
        self_sizes.append(round_to_digits(self.voxel_size * 1 / 2, DECIMALS))
        other_sizes = [round_to_digits(other.voxel_size * 2**i, DECIMALS) for i in range(other._octree_depth, -1, -1)]
        other_sizes.append(round_to_digits(other.voxel_size * 1 / 2, DECIMALS))

        self_stack = [(0, self._root_center, self._octree)]
        other_stack = [(0, other._root_center, other._octree)]

        while self_stack and other_stack:
            # Check for intersection until were sure there is / there isn't one
            self_current_depth, self_current_center, self_current_octree = self_stack.pop()
            other_current_depth, other_current_center, other_current_octree = other_stack.pop()

            if not self._check_voxel_intersection(
                self_current_center,
                self_sizes[self_current_depth + 1],
                other_current_center,
                other_sizes[other_current_depth + 1],
            ):
                # If these two voxels are not intersecting, we don't need to subdivide further
                continue

            if self_current_depth == self._octree_depth and other_current_depth == other._octree_depth:
                # If the voxel are intersecting and are leaves voxel, we are sure the voxelizations are intersecting
                return True

            self_new_stack = []
            other_new_stack = []

            if self_current_depth == self._octree_depth:
                # If it is a leaf voxel, we can't subdivide further, so we re-add it to the stack.
                self_new_stack.append((self_current_depth, self_current_center, self_current_octree))

            else:
                # We subdivide further
                half_size = self_sizes[self_current_depth + 1]

                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            if self_current_octree[i * 4 + j * 2 + k]:
                                # Check for child voxels
                                sub_voxel_center = round_point_3d_to_digits(
                                    (
                                        self_current_center[0] + (i - 0.5) * half_size,
                                        self_current_center[1] + (j - 0.5) * half_size,
                                        self_current_center[2] + (k - 0.5) * half_size,
                                    ),
                                    DECIMALS,
                                )

                                self_new_stack.append(
                                    (
                                        self_current_depth + 1,
                                        sub_voxel_center,
                                        self_current_octree[i * 4 + j * 2 + k],
                                    )
                                )

            if other_current_depth == other._octree_depth:
                # If it is a leaf voxel, we can't subdivide further, so we re-add it to the stack.
                other_new_stack.append((other_current_depth, other_current_center, other_current_octree))

            else:
                # We subdivide further
                half_size = other_sizes[other_current_depth + 1]

                for i in range(2):
                    for j in range(2):
                        for k in range(2):
                            if other_current_octree[i * 4 + j * 2 + k]:
                                sub_voxel_center = round_point_3d_to_digits(
                                    (
                                        other_current_center[0] + (i - 0.5) * half_size,
                                        other_current_center[1] + (j - 0.5) * half_size,
                                        other_current_center[2] + (k - 0.5) * half_size,
                                    ),
                                    DECIMALS,
                                )

                                other_new_stack.append(
                                    (
                                        other_current_depth + 1,
                                        sub_voxel_center,
                                        other_current_octree[i * 4 + j * 2 + k],
                                    )
                                )

            # We need to add all the permutations to the stack, to check there intersection.
            for self_voxel in self_new_stack:
                for other_voxel in other_new_stack:
                    self_stack.append(self_voxel)
                    other_stack.append(other_voxel)

        return False

    @staticmethod
    def _check_voxel_intersection(
        voxel_center_1: _Point3D, half_size_1: float, voxel_center_2: _Point3D, half_size_2: float
    ) -> bool:
        """
        Helper method to check if two voxels intersect.

        :param voxel_center_1: Center of the first voxel as (x, y, z) coordinates.
        :type voxel_center_1: tuple[float, float, float]
        :param half_size_1: The half edge size of the first voxel.
        :type half_size_1: float
        :param voxel_center_2: Center of the second voxel as (x, y, z) coordinates.
        :type voxel_center_2: tuple[float, float, float]
        :param half_size_2: The half edge size of the second voxel.
        :type half_size_2: float

        :return: True if the cubes intersect, False otherwise.
        :rtype: bool
        """
        # Calculate the minimum and maximum coordinates of each cube along each axis
        min_x1, max_x1 = voxel_center_1[0] - half_size_1, voxel_center_1[0] + half_size_1
        min_y1, max_y1 = voxel_center_1[1] - half_size_1, voxel_center_1[1] + half_size_1
        min_z1, max_z1 = voxel_center_1[2] - half_size_1, voxel_center_1[2] + half_size_1

        min_x2, max_x2 = voxel_center_2[0] - half_size_2, voxel_center_2[0] + half_size_2
        min_y2, max_y2 = voxel_center_2[1] - half_size_2, voxel_center_2[1] + half_size_2
        min_z2, max_z2 = voxel_center_2[2] - half_size_2, voxel_center_2[2] + half_size_2

        # Check for intersection along each axis
        x_intersect = max(min_x1, min_x2) < min(max_x1, max_x2)
        y_intersect = max(min_y1, min_y2) < min(max_y1, max_y2)
        z_intersect = max(min_z1, min_z2) < min(max_z1, max_z2)

        # The cubes intersect if they intersect along all three axes
        return x_intersect and y_intersect and z_intersect

    def union(self, other: "OctreeBasedVoxelization") -> "OctreeBasedVoxelization":
        """
        Perform a union operation with another OctreeBasedVoxelization.

        :param other: The OctreeBasedVoxelization to perform the union with.
        :type other: OctreeBasedVoxelization

        :return: A new OctreeBasedVoxelization resulting from the union operation.
        :rtype: OctreeBasedVoxelization
        """
        octree_1, octree_2 = self._make_octrees_same_depth(self, other)

        octree_union = self._recursive_union(
            0, octree_1._octree_depth, len(octree_1._triangles), octree_1._octree, octree_2._octree
        )
        triangles = self._triangles + other._triangles

        return self.__class__(octree_union, (0.0, 0.0, 0.0), octree_1._octree_depth, octree_1.voxel_size, triangles)

    @staticmethod
    def _recursive_union(
        current_depth: int,
        max_depth: int,
        n_triangles_1: int,
        current_octree_1,
        current_octree_2,
    ):
        """Recursive method to perform a union between two octree based voxelizations."""
        if current_depth == max_depth:  # if _octree_depth reached, it is a leaf node
            return current_octree_1 + [i_triangle + n_triangles_1 for i_triangle in current_octree_2]

        sub_voxels = []

        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if current_octree_1[i * 4 + j * 2 + k] and current_octree_2[i * 4 + j * 2 + k]:
                        # if it is in both octrees
                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_union(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                current_octree_1[i * 4 + j * 2 + k],
                                current_octree_2[i * 4 + j * 2 + k],
                            )
                        )

                    elif current_octree_1[i * 4 + j * 2 + k]:
                        # if it is in first octree only

                        if current_depth + 1 == max_depth:
                            _current_octree_2 = []
                        else:
                            _current_octree_2 = [[], [], [], [], [], [], [], []]

                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_union(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                current_octree_1[i * 4 + j * 2 + k],
                                _current_octree_2,
                            )
                        )

                    elif current_octree_2[i * 4 + j * 2 + k]:
                        # if it is in second octree only

                        if current_depth + 1 == max_depth:
                            _current_octree_1 = []
                        else:
                            _current_octree_1 = [[], [], [], [], [], [], [], []]

                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_union(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                _current_octree_1,
                                current_octree_2[i * 4 + j * 2 + k],
                            )
                        )

                    else:
                        sub_voxels.append([])

        return sub_voxels

    def difference(self, other: "OctreeBasedVoxelization") -> "OctreeBasedVoxelization":
        """
        Perform a difference operation with another voxelization.

        :param other: The OctreeBasedVoxelization to perform the difference with.
        :type other: OctreeBasedVoxelization

        :return: A new OctreeBasedVoxelization resulting from the difference operation.
        :rtype: OctreeBasedVoxelization
        """
        octree_1, octree_2 = self._make_octrees_same_depth(self, other)

        octree_difference = self._recursive_difference(
            0, octree_1._octree_depth, len(octree_1._triangles), octree_1._octree, octree_2._octree
        )
        triangles = self._triangles

        return self.__class__(
            octree_difference, (0.0, 0.0, 0.0), octree_1._octree_depth, octree_1.voxel_size, triangles
        )

    @staticmethod
    def _recursive_difference(
        current_depth: int,
        max_depth: int,
        n_triangles_1: int,
        current_octree_1,
        current_octree_2,
    ):
        """Recursive method to perform a difference between two octree based voxelizations."""
        if current_depth == max_depth:  # if _octree_depth reached, it is a leaf node
            if not current_octree_2:
                return current_octree_1
            return []

        sub_voxels = []

        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if current_octree_1[i * 4 + j * 2 + k] and current_octree_2[i * 4 + j * 2 + k]:
                        # if it is in both octrees
                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_difference(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                current_octree_1[i * 4 + j * 2 + k],
                                current_octree_2[i * 4 + j * 2 + k],
                            )
                        )

                    elif current_octree_1[i * 4 + j * 2 + k]:
                        # if it is in first octree only

                        if current_depth + 1 == max_depth:
                            _current_octree_2 = []
                        else:
                            _current_octree_2 = [[], [], [], [], [], [], [], []]

                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_difference(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                current_octree_1[i * 4 + j * 2 + k],
                                _current_octree_2,
                            )
                        )

                    else:
                        sub_voxels.append([])

        return sub_voxels

    def intersection(self, other: "OctreeBasedVoxelization") -> "OctreeBasedVoxelization":
        """
        Perform an intersection operation with another OctreeBasedVoxelization.

        :param other: The OctreeBasedVoxelization to perform the intersection with.
        :type other: OctreeBasedVoxelization

        :return: A new OctreeBasedVoxelization resulting from the intersection operation.
        :rtype: OctreeBasedVoxelization
        """
        octree_1, octree_2 = self._make_octrees_same_depth(self, other)

        octree_intersection = self._recursive_intersection(
            0, octree_1._octree_depth, len(octree_1._triangles), octree_1._octree, octree_2._octree
        )
        triangles = self._triangles + other._triangles

        return self.__class__(
            octree_intersection, (0.0, 0.0, 0.0), octree_1._octree_depth, octree_1.voxel_size, triangles
        )

    @staticmethod
    def _recursive_intersection(
        current_depth: int,
        max_depth: int,
        n_triangles_1: int,
        current_octree_1,
        current_octree_2,
    ):
        """Recursive method to perform a union between two octree based voxelizations."""
        if current_depth == max_depth:  # if _octree_depth reached, it is a leaf node
            return current_octree_1 + [i_triangle + n_triangles_1 for i_triangle in current_octree_2]

        sub_voxels = []

        for i in range(2):
            for j in range(2):
                for k in range(2):
                    # if it is in both octrees
                    if current_octree_1[i * 4 + j * 2 + k] and current_octree_2[i * 4 + j * 2 + k]:
                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_intersection(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                current_octree_1[i * 4 + j * 2 + k],
                                current_octree_2[i * 4 + j * 2 + k],
                            )
                        )

                    else:
                        sub_voxels.append([])

        return sub_voxels

    def symmetric_difference(self, other: "OctreeBasedVoxelization") -> "OctreeBasedVoxelization":
        """
        Perform a symmetric difference operation with another OctreeBasedVoxelization.

        :param other: The OctreeBasedVoxelization to perform the symmetric difference with.
        :type other: OctreeBasedVoxelization

        :return: A new OctreeBasedVoxelization resulting from the symmetric difference operation.
        :rtype: OctreeBasedVoxelization
        """
        octree_1, octree_2 = self._make_octrees_same_depth(self, other)

        octree_symmetric_difference = self._recursive_symmetric_difference(
            0, octree_1._octree_depth, len(octree_1._triangles), octree_1._octree, octree_2._octree
        )
        triangles = self._triangles + other._triangles

        return self.__class__(
            octree_symmetric_difference, (0.0, 0.0, 0.0), octree_1._octree_depth, octree_1.voxel_size, triangles
        )

    @staticmethod
    def _recursive_symmetric_difference(
        current_depth: int,
        max_depth: int,
        n_triangles_1: int,
        current_octree_1,
        current_octree_2,
    ):
        """Recursive method to perform a symmetric difference between two octree based voxelizations."""
        if current_depth == max_depth:  # if _octree_depth reached, it is a leaf node
            if not current_octree_2:
                return current_octree_1
            if not current_octree_1:
                return [i_triangle + n_triangles_1 for i_triangle in current_octree_2]
            return []

        sub_voxels = []

        for i in range(2):
            for j in range(2):
                for k in range(2):
                    if current_octree_1[i * 4 + j * 2 + k] and current_octree_2[i * 4 + j * 2 + k]:
                        # if it is in both octrees
                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_symmetric_difference(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                current_octree_1[i * 4 + j * 2 + k],
                                current_octree_2[i * 4 + j * 2 + k],
                            )
                        )

                    elif current_octree_1[i * 4 + j * 2 + k]:
                        # if it is in first octree only

                        if current_depth + 1 == max_depth:
                            _current_octree_2 = []
                        else:
                            _current_octree_2 = [[], [], [], [], [], [], [], []]

                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_symmetric_difference(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                current_octree_1[i * 4 + j * 2 + k],
                                _current_octree_2,
                            )
                        )

                    elif current_octree_2[i * 4 + j * 2 + k]:
                        # if it is in second octree only

                        if current_depth + 1 == max_depth:
                            _current_octree_1 = []
                        else:
                            _current_octree_1 = [[], [], [], [], [], [], [], []]

                        sub_voxels.append(
                            OctreeBasedVoxelization._recursive_symmetric_difference(
                                current_depth + 1,
                                max_depth,
                                n_triangles_1,
                                _current_octree_1,
                                current_octree_2[i * 4 + j * 2 + k],
                            )
                        )

                    else:
                        sub_voxels.append([])

        return sub_voxels

    def inverse(self) -> "OctreeBasedVoxelization":
        """
        Compute the inverse of the voxelization.

        :return: A new voxelization representing the inverse.
        :rtype: OctreeBasedVoxelization
        """
        return self.to_matrix_based_voxelization().inverse().to_octree_based_voxelization()

    def flood_fill(self, start: Tuple[int, int, int], fill_with: bool) -> "OctreeBasedVoxelization":
        """
        Perform a flood fill operation on the voxelization.

        :param start: The indexes of the starting voxel in the 3D matrix for the flood fill.
        :type start: tuple[int, int, int]
        :param fill_with: The value to fill the voxels with during the operation.
        :type fill_with: bool

        :return: A new voxelization resulting from the flood fill operation.
        :rtype: OctreeBasedVoxelization
        """
        return self.to_matrix_based_voxelization().flood_fill(start, fill_with).to_octree_based_voxelization()

    def _fill_outer_elements(self) -> "OctreeBasedVoxelization":
        """
        Fill the outer voxels of the voxelization.

        :return: A new voxelization with outer voxels filled.
        :rtype: OctreeBasedVoxelization
        """
        return self.to_matrix_based_voxelization().fill_outer_voxels().to_octree_based_voxelization()

    def _fill_enclosed_elements(self) -> "OctreeBasedVoxelization":
        """
        Fill the enclosed voxels of the voxelization.

        :return: A new voxelization with enclosed voxels filled.
        :rtype: OctreeBasedVoxelization
        """
        return self.to_matrix_based_voxelization().fill_enclosed_voxels().to_octree_based_voxelization()

    # EXPORT METHODS
    def to_point_based_voxelization(self) -> "PointBasedVoxelization":
        """
        Convert the OctreeBasedVoxelization to a PointBasedVoxelization.

        :return: A PointBasedVoxelization representation of the current voxelization.
        :rtype: PointBasedVoxelization
        """
        return PointBasedVoxelization(self.get_voxel_centers(), self.voxel_size, self.name)

    def to_matrix_based_voxelization(self) -> "MatrixBasedVoxelization":
        """
        Convert the OctreeBasedVoxelization to a MatrixBasedVoxelization.

        :return: A MatrixBasedVoxelization representation of the current voxelization.
        :rtype: MatrixBasedVoxelization
        """
        return self.to_point_based_voxelization().to_matrix_based_voxelization()

    def to_non_homogeneous_point_based_voxelizations(self) -> List["PointBasedVoxelization"]:
        """
        Convert the OctreeBasedVoxelization to multiple PointBasedVoxelization, with different size.

        :return: A PointBasedVoxelization representation of the current voxelization.
        :rtype: PointBasedVoxelization
        """
        point_based_voxelizations = []

        for voxel_size, voxel_centers in self._get_non_homogeneous_voxel_centers().items():
            point_based_voxelizations.append(PointBasedVoxelization(voxel_centers, voxel_size))

        return point_based_voxelizations

    def to_inner_growing_voxelizations(self, layers_minimal_thickness: float) -> List["OctreeBasedVoxelization"]:
        """
        Convert the OctreeBasedVoxelization to multiple OctreeBasedVoxelization, with different size.

        The more the voxelization is inside, the more its voxel size become bigger.

        :param layers_minimal_thickness: The minimal thickness of each layer.
        :type layers_minimal_thickness: float

        :return: A list of PointBasedVoxelization representing the inner growing voxelization.
        :rtype: List[PointBasedVoxelization]
        """
        i = 2
        inner_growing_voxel_centers = self.get_inner_growing_voxel_centers(layers_minimal_thickness)

        while len(inner_growing_voxel_centers.keys()) == i:
            # Still making the voxelization inner growing
            max_voxel_size = max(inner_growing_voxel_centers.keys())
            i += 1

            for voxel_size, voxel_centers in (
                OctreeBasedVoxelization.from_point_based_voxelization(
                    PointBasedVoxelization(inner_growing_voxel_centers[max_voxel_size], max_voxel_size)
                )
                .get_inner_growing_voxel_centers(layers_minimal_thickness)
                .items()
            ):
                inner_growing_voxel_centers[voxel_size] = voxel_centers

        octree_based_voxelizations = []

        for voxel_size, voxel_centers in inner_growing_voxel_centers.items():
            octree_based_voxelizations.append(
                PointBasedVoxelization(voxel_centers, voxel_size).to_octree_based_voxelization()
            )

        return octree_based_voxelizations

    # HELPER CREATION METHODS
    @staticmethod
    def _subdivide_from_triangles(
        triangles: List[_Triangle3D],
        intersecting_indices: List[int],
        center: _Point3D,
        sizes: List[float],
        depth: int,
        max_depth: int,
    ) -> Octree:
        """Recursive method to create an OctreeBasedVoxelization from a list of triangles."""
        if depth < max_depth:  # not yet reached max depth
            half_size = sizes[depth + 1]
            quarter_size = sizes[depth + 2]

            sub_voxels = []

            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        # calculate the center of the sub-voxel
                        sub_voxel_center = round_point_3d_to_digits(
                            (
                                center[0] + (i - 0.5) * half_size,
                                center[1] + (j - 0.5) * half_size,
                                center[2] + (k - 0.5) * half_size,
                            ),
                            DECIMALS,
                        )

                        # check for intersecting triangle with the sub-voxel
                        sub_voxel_intersecting_indices = [
                            i
                            for i in intersecting_indices
                            if triangle_intersects_voxel(
                                triangles[i], sub_voxel_center, (quarter_size, quarter_size, quarter_size)
                            )
                        ]

                        # Recursive process
                        if not sub_voxel_intersecting_indices:
                            # If sub-voxel not intersecting
                            sub_voxels.append([])

                        else:
                            # If sub-voxel intersecting
                            sub_voxels.append(
                                OctreeBasedVoxelization._subdivide_from_triangles(
                                    triangles=triangles,
                                    intersecting_indices=sub_voxel_intersecting_indices,
                                    center=sub_voxel_center,
                                    sizes=sizes,
                                    depth=depth + 1,
                                    max_depth=max_depth,
                                )
                            )

            if all(not sub_voxel for sub_voxel in sub_voxels):
                return []
            return sub_voxels

        # reached max depth
        return intersecting_indices

    @staticmethod
    def _subdivide_from_points(
        points: List[_Point3D],
        center: _Point3D,
        sizes: List[float],
        depth: int,
        max_depth: int,
    ) -> Octree:
        """Recursive method to create an OctreeBasedVoxelization from a list of points."""

        if depth < max_depth:  # not yet reached max depth
            half_size = sizes[depth + 1]

            sub_voxels = []

            # Initialize lists for sub-voxel points
            sub_voxel_points = [[] for _ in range(8)]

            # Check each point and determine which sub-voxel it belongs to
            for point in points:
                # Relative position to the center
                idx = 0
                if point[0] > center[0]:
                    idx |= 4
                if point[1] > center[1]:
                    idx |= 2
                if point[2] > center[2]:
                    idx |= 1

                sub_voxel_points[idx].append(point)

            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        # calculate the center of the sub-voxel
                        sub_voxel_center = round_point_3d_to_digits(
                            (
                                center[0] + (i - 0.5) * half_size,
                                center[1] + (j - 0.5) * half_size,
                                center[2] + (k - 0.5) * half_size,
                            ),
                            DECIMALS,
                        )

                        idx = (i << 2) + (j << 1) + k

                        if not sub_voxel_points[idx]:
                            sub_voxels.append([])

                        else:
                            sub_voxels.append(
                                OctreeBasedVoxelization._subdivide_from_points(
                                    points=sub_voxel_points[idx],
                                    center=sub_voxel_center,
                                    sizes=sizes,
                                    depth=depth + 1,
                                    max_depth=max_depth,
                                )
                            )

            if all(not sub_voxel for sub_voxel in sub_voxels):
                return []
            return sub_voxels

        # reached max depth
        return [1]

    # HELPER EXPORT METHODS
    def _get_homogeneous_leaf_centers(
        self, current_depth: int, current_size: float, current_center: _Point3D, current_octree
    ) -> Set[_Point3D]:
        """Recursive method to extract all the leaf voxel center (voxels of minimal size)."""
        if current_depth == self._octree_depth:  # if _octree_depth reached, it is a leaf node
            return {current_center}

        centers = set()
        half_size = round_to_digits(current_size / 2, DECIMALS)

        for i in range(2):
            for j in range(2):
                for k in range(2):
                    # if it is the octree
                    if current_octree[i * 4 + j * 2 + k]:
                        # calculate the center of the sub-voxel
                        sub_voxel_center = round_point_3d_to_digits(
                            (
                                current_center[0] + (i - 0.5) * half_size,
                                current_center[1] + (j - 0.5) * half_size,
                                current_center[2] + (k - 0.5) * half_size,
                            ),
                            DECIMALS,
                        )

                        centers = centers.union(
                            self._get_homogeneous_leaf_centers(
                                current_depth + 1, half_size, sub_voxel_center, current_octree[i * 4 + j * 2 + k]
                            )
                        )

        return centers

    def _get_non_homogeneous_voxel_centers(self) -> Dict[float, Set[_Point3D]]:
        """
        Get the center points of non-homogeneous voxels and organize them by voxel size.

        This method maximize the voxel size without any less of information.

        :return: A dictionary where the keys are voxel sizes and the values are sets of voxel centers.
        :rtype: dict[float, set[tuple[float, float, float]]]
        """
        return self._get_non_homogeneous_leaf_centers(0, self._root_voxel_size, self._root_center, self._octree)

    def _get_non_homogeneous_leaf_centers(
        self, current_depth: int, current_size: float, current_center: _Point3D, current_octree
    ) -> Dict[float, Set[_Point3D]]:
        """
        Recursive method to extract all the non-homogeneous voxel centers.

        This method maximize the voxel size without any less of information.

        :return: A dictionary where the keys are voxel sizes and the values are sets of voxel centers.
        :rtype: dict[float, set[tuple[float, float, float]]]
        """
        centers_by_voxel_size = {}

        if current_depth == self._octree_depth:
            # Return leaf nodes
            centers_by_voxel_size[current_size] = {current_center}

        else:
            half_size = round_to_digits(current_size / 2, DECIMALS)

            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        # if it is the octree
                        if current_octree[i * 4 + j * 2 + k]:
                            # calculate the center of the sub-voxel
                            sub_voxel_center = round_point_3d_to_digits(
                                (
                                    current_center[0] + (i - 0.5) * half_size,
                                    current_center[1] + (j - 0.5) * half_size,
                                    current_center[2] + (k - 0.5) * half_size,
                                ),
                                DECIMALS,
                            )
                            # Recursive process
                            sub_centers = self._get_non_homogeneous_leaf_centers(
                                current_depth + 1, half_size, sub_voxel_center, current_octree[i * 4 + j * 2 + k]
                            )

                            # Merge sub-centers into the result dictionary, keeping track of voxel sizes
                            for size, sub_voxel_centers in sub_centers.items():
                                if size not in centers_by_voxel_size:
                                    centers_by_voxel_size[size] = set()
                                centers_by_voxel_size[size].update(sub_voxel_centers)

            if len(centers_by_voxel_size.get(half_size, [])) == 8:
                # Merge voxels if possible
                del centers_by_voxel_size[half_size]
                centers_by_voxel_size[current_size] = {current_center}

        return centers_by_voxel_size

    def get_inner_growing_voxel_centers(
        self, layers_minimal_thickness: float, layer_dict: Dict[_Point3D, float] = None
    ) -> Dict[float, Set[_Point3D]]:
        """
        Get the center points of inner growing voxels and organize them by voxel size.

        :param layers_minimal_thickness: The minimal thickness of each layer.
        :type layers_minimal_thickness: float
        :param layer_dict: A dict with voxel centers and their corresponding layers thickness to false or edge values.
        :type layer_dict: dict[tuple[float, float, float], float]

        :return: A dictionary where the keys are voxel sizes and the values are sets of voxel centers.
        :rtype: dict[float, set[tuple[float, float, float]]]
        """
        if not layer_dict:
            layer_dict = self.to_matrix_based_voxelization().layers_thickness_by_voxel_centers()

        return self._get_inner_growing_leaf_centers(
            0,
            self._root_voxel_size,
            self._root_center,
            self._octree,
            layer_dict,
            layers_minimal_thickness,
        )

    def _get_inner_growing_leaf_centers(
        self,
        current_depth: int,
        current_size: float,
        current_center: _Point3D,
        current_octree,
        layer_dict: Dict[_Point3D, float],
        min_layer_thickness: float,
    ) -> Dict[float, Set[_Point3D]]:
        """Recursive method to extract inner growing voxel centers."""
        centers_by_voxel_size = {}

        if current_depth == self._octree_depth:
            # Return leaf nodes
            centers_by_voxel_size[current_size] = {current_center}

        else:
            half_size = round_to_digits(current_size / 2, DECIMALS)

            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        # if it is the octree
                        if current_octree[i * 4 + j * 2 + k]:
                            # calculate the center of the sub-voxel
                            sub_voxel_center = round_point_3d_to_digits(
                                (
                                    current_center[0] + (i - 0.5) * half_size,
                                    current_center[1] + (j - 0.5) * half_size,
                                    current_center[2] + (k - 0.5) * half_size,
                                ),
                                DECIMALS,
                            )

                            # Recursive process
                            sub_centers = self._get_inner_growing_leaf_centers(
                                current_depth + 1,
                                half_size,
                                sub_voxel_center,
                                current_octree[i * 4 + j * 2 + k],
                                layer_dict,
                                min_layer_thickness,
                            )

                            # Merge sub-centers into the result dictionary, keeping track of voxel sizes
                            for size, sub_voxel_centers in sub_centers.items():
                                if size not in centers_by_voxel_size:
                                    centers_by_voxel_size[size] = set()
                                centers_by_voxel_size[size].update(sub_voxel_centers)

            if len(centers_by_voxel_size.get(half_size, [])) == 8:
                if current_depth + 1 == self._octree_depth and all(
                    layer_dict[voxel_center] >= min_layer_thickness for voxel_center in centers_by_voxel_size[half_size]
                ):
                    # Merge voxels
                    centers_by_voxel_size[current_size] = {current_center}
                    del centers_by_voxel_size[half_size]

        return centers_by_voxel_size

    # FACES INTERSECIONS METHODS
    @classmethod
    def intersecting_faces_combinations(
        cls, shell_1: Shell3D, shell_2: Shell3D, voxel_size: float
    ) -> List[Tuple[Tuple[Face3D, Face3D], PointBasedVoxelization]]:
        """
        Compute the intersecting faces combinations and where the faces are located, as a PointBasedVoxelization.

        :param shell_1: The first shell to find the intersecting faces with.
        :type shell_1: Shell3D
        :param shell_1: The second shell to find the intersecting faces with.
        :type shell_2: Shell3D
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: The possibly intersecting face combinations.
        :rtype: list[tuple[tuple[Face3D, Face3D], PointBasedVoxelization]]
        """
        face_idx_by_triangle_1, shell_triangles_1 = cls._shell_to_face_idx_by_triangle(shell_1)
        face_idx_by_triangle_2, shell_triangles_2 = cls._shell_to_face_idx_by_triangle(shell_2)

        voxelization_1 = cls.from_triangles(shell_triangles_1, voxel_size)
        voxelization_2 = cls.from_triangles(shell_triangles_2, voxel_size)

        intersection = voxelization_1.intersection(voxelization_2)
        triangle_combinations = intersection._get_intersections_voxel_centers(len(voxelization_1._triangles))

        face_id_combinations = {}

        for (i, j), voxel_centers in triangle_combinations.items():
            face_1 = face_idx_by_triangle_1[shell_triangles_1[i]]
            face_2 = face_idx_by_triangle_2[shell_triangles_2[j - len(shell_triangles_1)]]

            if (face_1, face_2) not in face_id_combinations:
                face_id_combinations[(face_1, face_2)] = set()
            face_id_combinations[(face_1, face_2)].update(voxel_centers)

        face_combinations = []

        for (i, j), voxel_centers in face_id_combinations.items():
            face_combinations.append(
                ((shell_1.faces[i], shell_2.faces[j]), PointBasedVoxelization(voxel_centers, voxel_size))
            )

        return face_combinations

    def _get_intersections_voxel_centers(self, threshold: int = None) -> Dict[Tuple[int, int], Set[_Point3D]]:
        """
        Get the center points of non-homogeneous voxels and organize them by voxel size.

        :return: A dictionary where the keys are voxel sizes and the values are sets of voxel centers.
        :rtype: dict[float, set[tuple[float, float, float]]]
        """
        return self._get_intersections_leaf_centers(
            0, self._root_voxel_size, self._root_center, self._octree, threshold
        )

    def _get_intersections_leaf_centers(
        self,
        current_depth: int,
        current_size: float,
        current_center: _Point3D,
        current_octree,
        threshold: int = None,
    ) -> Dict[Tuple[int, int], Set[_Point3D]]:
        """Recursive method to extract all the leaf voxel center (voxels of minimal size)."""
        intersections_locations = {}

        if current_depth == self._octree_depth:  # if _octree_depth reached, it is a leaf node
            for combination in itertools.combinations(current_octree, 2):
                combination = tuple(sorted(combination))

                if threshold:
                    if not combination[0] < threshold <= combination[1]:
                        continue

                intersections_locations[combination] = {current_center}

        else:
            half_size = round_to_digits(current_size / 2, DECIMALS)

            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        # if it is the octree
                        if current_octree[i * 4 + j * 2 + k]:
                            # calculate the center of the sub-voxel
                            sub_voxel_center = round_point_3d_to_digits(
                                (
                                    current_center[0] + (i - 0.5) * half_size,
                                    current_center[1] + (j - 0.5) * half_size,
                                    current_center[2] + (k - 0.5) * half_size,
                                ),
                                DECIMALS,
                            )

                            # Recursive process
                            sub_intersections_locations = self._get_intersections_leaf_centers(
                                current_depth + 1,
                                half_size,
                                sub_voxel_center,
                                current_octree[i * 4 + j * 2 + k],
                                threshold,
                            )

                            # Merge sub-centers into the result dictionary, keeping track of voxel sizes
                            for intersection, location in sub_intersections_locations.items():
                                if intersection not in intersections_locations:
                                    intersections_locations[intersection] = set()
                                intersections_locations[intersection].update(location)

        return intersections_locations

    @staticmethod
    def _shell_to_face_idx_by_triangle(shell: Shell3D):
        """Helper method to triangulate a Shell3D while keeping its faces correspondence with the triangles."""
        face_idx_by_triangle = {}
        shell_triangles = []

        for i, face in enumerate(shell.faces):
            try:
                triangulation = face.triangulation()

                face_triangles = [
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

                for triangle in face_triangles:
                    face_idx_by_triangle[triangle] = i

                shell_triangles.extend(face_triangles)

            except RuntimeError:
                warnings.warn(f"Failed triangulation of {face}.")

        return face_idx_by_triangle, shell_triangles

    # HELPER METHODS
    def _increment_octree_depth(self) -> "OctreeBasedVoxelization":
        """
        Increment the octree depth by doubling the size of the root voxel.

        :return: The modified octree based voxelization.
        :rtype: OctreeBasedVoxelization
        """
        if self._root_center != (0.0, 0.0, 0.0):
            raise NotImplementedError

        new_octree = [
            [[], [], [], [], [], [], [], self._octree[0]],
            [[], [], [], [], [], [], self._octree[1], []],
            [[], [], [], [], [], self._octree[2], [], []],
            [[], [], [], [], self._octree[3], [], [], []],
            [[], [], [], self._octree[4], [], [], [], []],
            [[], [], self._octree[5], [], [], [], [], []],
            [[], self._octree[6], [], [], [], [], [], []],
            [self._octree[7], [], [], [], [], [], [], []],
        ]

        return self.__class__(
            new_octree, self._root_center, self._octree_depth + 1, self.voxel_size, self._triangles, self.name
        )

    @staticmethod
    def _make_octrees_same_depth(
        octree_1: "OctreeBasedVoxelization", octree_2: "OctreeBasedVoxelization"
    ) -> Tuple["OctreeBasedVoxelization", "OctreeBasedVoxelization"]:
        """
        Helper method to make two OctreeBasedVoxelization having the same octree depth.

        The OctreeBasedVoxelization must have same voxel size and octree root.

        :param octree_1: The first OctreeBasedVoxelization.
        :type octree_1: OctreeBasedVoxelization
        :param octree_2: The second OctreeBasedVoxelization.
        :type octree_2: OctreeBasedVoxelization

        :return: The two OctreeBasedVoxelization modified to have the same depth.
        :rtype: tuple[OctreeBasedVoxelization, OctreeBasedVoxelization]
        """
        octree_1._check_other_type(octree_2)
        octree_1._check_other_element_size(octree_2)

        # Make the octree have the same depth
        while octree_1._octree_depth < octree_2._octree_depth:
            octree_1 = octree_1._increment_octree_depth()

        while octree_2._octree_depth < octree_1._octree_depth:
            octree_2 = octree_2._increment_octree_depth()

        return octree_1, octree_2

    # SERIALIZATION
    def to_dict(
        self, use_pointers: bool = True, memo=None, path: str = "#", id_method=True, id_memo=None
    ) -> JsonSerializable:
        """Specific 'to_dict' method."""
        dict_ = self.base_dict()

        dict_["octree"] = self._octree
        dict_["root_center"] = self._root_center
        dict_["octree_depth"] = self._octree_depth
        dict_["voxel_size"] = self.voxel_size
        dict_["triangles"] = self._triangles
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
    ) -> "OctreeBasedVoxelization":
        """Specific 'dict_to_object' method."""

        octree = dict_["octree"]
        root_center = dict_["root_center"]
        octree_depth = dict_["octree_depth"]
        voxel_size = dict_["voxel_size"]
        triangles = dict_["triangles"]
        name = dict_["name"]

        return cls(octree, root_center, octree_depth, voxel_size, triangles, name)


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
        raise NotImplementedError("Pixelization is an abstract class and should not be use directly.")

    @classmethod
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
        raise NotImplementedError("Pixelization is an abstract class and should not be use directly.")

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

        new_shape = np.round((global_max - global_min) / self.pixel_size, DECIMALS).astype(int)

        new_self = np.zeros(new_shape, dtype=np.bool_)
        new_other = np.zeros(new_shape, dtype=np.bool_)

        self_start = np.round((self_min - global_min) / self.pixel_size, DECIMALS).astype(int)
        other_start = np.round((other_min - global_min) / self.pixel_size, DECIMALS).astype(int)

        new_self[
            self_start[0] : self_start[0] + self.matrix.shape[0],
            self_start[1] : self_start[1] + self.matrix.shape[1],
        ] = self.matrix

        new_other[
            other_start[0] : other_start[0] + other.matrix.shape[0],
            other_start[1] : other_start[1] + other.matrix.shape[1],
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
            min_pixel_coords[0] : max_pixel_coords[0] + 1,
            min_pixel_coords[1] : max_pixel_coords[1] + 1,
        ]

        # Calculate new matrix_origin_center
        new_origin_center = np.round(self.min_grid_center + min_pixel_coords * self.pixel_size, DECIMALS)

        return self.__class__(cropped_matrix, tuple(new_origin_center), self.pixel_size, self.name)
