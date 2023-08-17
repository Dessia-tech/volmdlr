# cython: language_level=3
# distutils: language = c++
"""
Pure python module to define cython function.
This module needs to be compiled!
"""
from typing import List, Set, Tuple

import cython
import cython.cimports.libc.math as math_c
import numpy
import numpy as np
from cython.cimports.libcpp import bool as bool_C
from cython.cimports.libcpp.stack import stack
from cython.cimports.libcpp.vector import vector

# TODO: refactor, add docstrings

# CUSTOM PYTHON TYPES

Point = Tuple[float, float, float]
Triangle = Tuple[Point, Point, Point]

# PYTHON FUNCTIONS


def triangles_to_voxels(triangles: List[Triangle], voxel_size: float) -> Set[Point]:
    """
    Helper function to compute all the voxels intersecting with a given list of triangles.

    :param triangles: The triangles to compute the intersecting voxels.
    :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
    :param voxel_size: The voxel edges size.
    :type voxel_size: float

    :return: The centers of the voxels that intersect with the triangles.
    :rtype: set[tuple[float, float, float]]
    """
    # TODO: replace this method with triangles_to_voxel_matrix
    voxel_centers = set()

    for triangle in triangles:
        min_point = tuple(min(p[i] for p in triangle) for i in range(3))
        max_point = tuple(max(p[i] for p in triangle) for i in range(3))

        for bbox_center in _aabb_intersecting_boxes(
            min_point,
            max_point,
            voxel_size,
        ):
            bbox_center = tuple(bbox_center)
            if bbox_center not in voxel_centers:
                if _triangle_intersects_voxel(
                    triangle,
                    bbox_center,
                    (0.5 * voxel_size, 0.5 * voxel_size, 0.5 * voxel_size),
                ):
                    voxel_centers.add(bbox_center)

    return voxel_centers


def triangles_to_voxel_matrix(
    triangles: List[Tuple[Tuple[float, float, float]], Tuple[float, float, float], Tuple[float, float, float]],
    voxel_size: float,
) -> Tuple[np.ndarray[np.bool_, np.ndim == 3], Tuple[float, float, float]]:
    """
    Helper function to compute the voxel matrix of all the voxels intersecting with a given list of triangles.

    :param triangles: The triangles to compute the intersecting voxels.
    :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
    :param voxel_size: The voxel edges size.
    :type voxel_size: float

    :return: The voxel matrix and the origin center of the matrix.
    :rtype: tuple[np.ndarray[np.bool_, np.ndim == 3], tuple[float, float, float]]
    """
    # compute the size of the matrix and min matrix origin center
    min_point, max_point = _triangles_min_max_points(triangles)
    shape = (
        int(max_point[0] // voxel_size + 1) - int(min_point[0] // voxel_size) + 2,
        int(max_point[1] // voxel_size + 1) - int(min_point[1] // voxel_size) + 2,
        int(max_point[2] // voxel_size + 1) - int(min_point[2] // voxel_size) + 2,
    )
    matrix = numpy.zeros(shape, dtype=np.bool_)
    matrix_origin_center = (
        round((min_point[0] // voxel_size - 0.5) * voxel_size, 6),
        round((min_point[1] // voxel_size - 0.5) * voxel_size, 6),
        round((min_point[2] // voxel_size - 0.5) * voxel_size, 6),
    )

    # compute the intersecting voxel
    return (
        np.asarray(
            _triangles_to_voxel_matrix(
                np.array(triangles),
                len(triangles),
                voxel_size,
                matrix,
                matrix_origin_center,
            ),
            np.bool_,
        ),
        matrix_origin_center,
    )


def flood_fill_matrix_2d(
    matrix: np.ndarray[np.bool_, np.ndim == 2], start: Tuple[int, int], fill_with: bool
) -> np.ndarray[np.bool_, np.ndim == 2]:
    """
    Perform 2D flood fill on the given matrix starting from the specified position.

    :param matrix: The matrix to perform flood fill on.
    :type matrix: np.ndarray[np.bool_, np.ndim == 2]
    :param start: The starting position for flood fill.
    :type start: Tuple[int, int]
    :param fill_with: The value to fill the connected component with.
    :type fill_with: bool

    :return: The matrix after performing flood fill.
    :rtype: np.ndarray[np.bool_, np.ndim == 2]
    """
    return np.asarray(
        _flood_fill_matrix_2d(
            matrix.astype(np.bool_),
            (start[0], start[1]),
            fill_with,
            (matrix.shape[0], matrix.shape[1]),
        ),
        dtype=np.bool_,
    )


def flood_fill_matrix_3d(
    matrix: np.ndarray[np.bool_, np.ndim == 3], start: Tuple[int, int, int], fill_with: bool
) -> np.ndarray[np.bool_, np.ndim == 3]:
    """
    Perform 3D flood fill on the given matrix starting from the specified position.

    :param matrix: The matrix to perform flood fill on.
    :type matrix: np.ndarray[np.bool_, np.ndim == 3]
    :param start: The starting position for flood fill.
    :type start: tuple[int, int, int]
    :param fill_with: The value to fill the connected component with.
    :type fill_with: bool

    :return: The matrix after performing flood fill.
    :rtype: np.ndarray[np.bool_, np.ndim == 3]
    """
    return np.asarray(
        _flood_fill_matrix_3d(
            matrix.astype(np.bool_),
            (start[0], start[1], start[2]),
            fill_with,
            (matrix.shape[0], matrix.shape[1], matrix.shape[2]),
        ),
        dtype=np.bool_,
    )


def line_segments_to_pixels(
    line_segments: List[Tuple[Tuple[float, float], Tuple[float, float]]], pixel_size: float
) -> Set[Tuple[float, float]]:
    """
    Convert a list of line segments to a set of pixel coordinates.

    This function takes a list of line segments, each defined by two endpoints, and converts them into a set of pixel
    coordinates based on the specified pixel size. The pixel coordinates represent the pixels that the line segments
    intersect with.

    :param line_segments: List of line segments to convert.
    :type line_segments: list[tuple[tuple[float, float], tuple[float, float]]]
    :param pixel_size: The size of each pixel.
    :type pixel_size: float

    :return: A set of pixel coordinates representing the intersection points of the line segments.
    :rtype: set[tuple[float, float]]
    """
    return set(_line_segments_to_pixels(line_segments, pixel_size))


# CYTHON FUNCTIONS


@cython.cfunc
@cython.cdivision(True)
@cython.exceptval(check=False)
def _round_to_digits(num: cython.double, digits: cython.int) -> cython.double:
    """Round the given number to the specified number of digits after the decimal point."""
    multiplier: cython.double = math_c.pow(10.0, digits)
    return math_c.round(num * multiplier) / multiplier


@cython.cfunc
@cython.exceptval(check=False)
def _is_integer(value: cython.double) -> bool_C:
    """Check if the given value is an integer (has no fractional part)."""
    return cython.cast(cython.int, value) == value


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def _triangle_intersects_voxel(
    triangle: Tuple[
        Tuple[cython.double, cython.double, cython.double],
        Tuple[cython.double, cython.double, cython.double],
        Tuple[cython.double, cython.double, cython.double],
    ],
    voxel_center: Tuple[cython.double, cython.double, cython.double],
    voxel_extents: Tuple[cython.double, cython.double, cython.double],
) -> bool_C:
    """Check if a 3D triangle intersects with a voxel defined by its center and extents."""

    # Ported from https://gist.github.com/zvonicek/fe73ba9903f49d57314cf7e8e0f05dcf

    v0: cython.double[3]
    v1: cython.double[3]
    v2: cython.double[3]
    f0: cython.double[3]
    f1: cython.double[3]
    f2: cython.double[3]
    box_center: cython.double[3]
    box_extents: cython.double[3]
    plane_normal: cython.double[3]
    plane_distance: cython.double
    r: cython.double
    a00: cython.double[3]
    a01: cython.double[3]
    a02: cython.double[3]
    a10: cython.double[3]
    a11: cython.double[3]
    a12: cython.double[3]
    a20: cython.double[3]
    a21: cython.double[3]
    a22: cython.double[3]
    p0: cython.double
    p1: cython.double
    p2: cython.double

    # Translate triangle as conceptually moving AABB to origin
    v0[0] = triangle[0][0] - voxel_center[0]
    v0[1] = triangle[0][1] - voxel_center[1]
    v0[2] = triangle[0][2] - voxel_center[2]

    v1[0] = triangle[1][0] - voxel_center[0]
    v1[1] = triangle[1][1] - voxel_center[1]
    v1[2] = triangle[1][2] - voxel_center[2]

    v2[0] = triangle[2][0] - voxel_center[0]
    v2[1] = triangle[2][1] - voxel_center[1]
    v2[2] = triangle[2][2] - voxel_center[2]

    # Compute edge vectors for triangle
    f0[0] = triangle[1][0] - triangle[0][0]
    f0[1] = triangle[1][1] - triangle[0][1]
    f0[2] = triangle[1][2] - triangle[0][2]

    f1[0] = triangle[2][0] - triangle[1][0]
    f1[1] = triangle[2][1] - triangle[1][1]
    f1[2] = triangle[2][2] - triangle[1][2]

    f2[0] = triangle[0][0] - triangle[2][0]
    f2[1] = triangle[0][1] - triangle[2][1]
    f2[2] = triangle[0][2] - triangle[2][2]

    # REGION TEST THE THREE AXES CORRESPONDING TO THE FACE NORMALS OF AABB B (CATEGORY 1)

    # Exit if...
    # ... [-extents.X, extents.X] and [min(v0.X,v1.X,v2.X), max(v0.X,v1.X,v2.X)] do not overlap
    if max(v0[0], v1[0], v2[0]) < -voxel_extents[0] or min(v0[0], v1[0], v2[0]) > voxel_extents[0]:
        return False

    # ... [-extents.Y, extents.Y] and [min(v0.Y,v1.Y,v2.Y), max(v0.Y,v1.Y,v2.Y)] do not overlap
    if max(v0[1], v1[1], v2[1]) < -voxel_extents[1] or min(v0[1], v1[1], v2[1]) > voxel_extents[1]:
        return False

    # ... [-extents.Z, extents.Z] and [min(v0.Z,v1.Z,v2.Z), max(v0.Z,v1.Z,v2.Z)] do not overlap
    if max(v0[2], v1[2], v2[2]) < -voxel_extents[2] or min(v0[2], v1[2], v2[2]) > voxel_extents[2]:
        return False

    # ENDREGION

    # REGION TEST SEPARATING AXIS CORRESPONDING TO TRIANGLE FACE NORMAL (CATEGORY 2)

    plane_normal[0] = f0[1] * f1[2] - f0[2] * f1[1]
    plane_normal[1] = f0[2] * f1[0] - f0[0] * f1[2]
    plane_normal[2] = f0[0] * f1[1] - f0[1] * f1[0]

    plane_distance = math_c.fabs(plane_normal[0] * v0[0] + plane_normal[1] * v0[1] + plane_normal[2] * v0[2])

    # Compute the projection interval radius of b onto L(t) = b.c + t * p.n
    r = (
        voxel_extents[0] * math_c.fabs(plane_normal[0])
        + voxel_extents[1] * math_c.fabs(plane_normal[1])
        + voxel_extents[2] * math_c.fabs(plane_normal[2])
    )

    # Intersection occurs when plane distance falls within [-r,+r] interval
    if plane_distance > r:
        return False

    # ENDREGION

    # REGION TEST AXES a00..a22 (CATEGORY 3)

    # Test axis a00
    a00[0] = 0
    a00[1] = -f0[2]
    a00[2] = f0[1]
    if _calculate_axis_values(v0, v1, v2, a00, f0, voxel_extents):
        return False

    # Test axis a01
    a01[0] = 0
    a01[1] = -f1[2]
    a01[2] = f1[1]
    if _calculate_axis_values(v0, v1, v2, a01, f1, voxel_extents):
        return False

    # Test axis a02
    a02[0] = 0
    a02[1] = -f2[2]
    a02[2] = f2[1]
    if _calculate_axis_values(v0, v1, v2, a02, f2, voxel_extents):
        return False

    # Test axis a10
    a10[0] = f0[2]
    a10[1] = 0
    a10[2] = -f0[0]
    if _calculate_axis_values(v0, v1, v2, a10, f0, voxel_extents):
        return False

    # Test axis a11
    a11[0] = f1[2]
    a11[1] = 0
    a11[2] = -f1[0]
    if _calculate_axis_values(v0, v1, v2, a11, f1, voxel_extents):
        return False

    # Test axis a12
    a12[0] = f2[2]
    a12[1] = 0
    a12[2] = -f2[0]
    if _calculate_axis_values(v0, v1, v2, a12, f2, voxel_extents):
        return False

    # Test axis a20
    a20[0] = -f0[1]
    a20[1] = f0[0]
    a20[2] = 0
    if _calculate_axis_values(v0, v1, v2, a20, f0, voxel_extents):
        return False

    # Test axis a21
    a21[0] = -f1[1]
    a21[1] = f1[0]
    a21[2] = 0
    if _calculate_axis_values(v0, v1, v2, a21, f1, voxel_extents):
        return False

    # Test axis a22
    a22[0] = -f2[1]
    a22[1] = f2[0]
    a22[2] = 0
    if _calculate_axis_values(v0, v1, v2, a22, f2, voxel_extents):
        return False

    # ENDREGION

    return True


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def _calculate_axis_values(
    v0: cython.double[3],
    v1: cython.double[3],
    v2: cython.double[3],
    ax: cython.double[3],
    f: cython.double[3],
    voxel_extents: Tuple[cython.double, cython.double, cython.double],
) -> bool_C:
    """Calculate axis values used in triangle intersection tests with an axis-aligned box."""

    p0 = v0[0] * ax[0] + v0[1] * ax[1] + v0[2] * ax[2]
    p1 = v1[0] * ax[0] + v1[1] * ax[1] + v1[2] * ax[2]
    p2 = v2[0] * ax[0] + v2[1] * ax[1] + v2[2] * ax[2]
    r = (
        voxel_extents[0] * math_c.fabs(f[2])
        + voxel_extents[1] * math_c.fabs(f[0])
        + voxel_extents[2] * math_c.fabs(f[1])
    )

    return max(-max(p0, p1, p2), min(p0, p1, p2)) > r


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _aabb_intersecting_boxes(
    min_point: Tuple[cython.double, cython.double, cython.double],
    max_point: Tuple[cython.double, cython.double, cython.double],
    voxel_size: cython.double,
) -> vector[Tuple[cython.double, cython.double, cython.double]]:
    """Calculate the voxel centers that intestects with a given axis-aligned boxes."""

    x_start: cython.int
    x_end: cython.int
    y_start: cython.int
    y_end: cython.int
    z_start: cython.int
    z_end: cython.int
    x: cython.int
    y: cython.int
    z: cython.int
    num_centers: cython.int

    x_start = cython.cast(cython.int, (min_point[0] / voxel_size) - 1)
    x_end = cython.cast(cython.int, (max_point[0] / voxel_size) + 1)
    y_start = cython.cast(cython.int, (min_point[1] / voxel_size) - 1)
    y_end = cython.cast(cython.int, (max_point[1] / voxel_size) + 1)
    z_start = cython.cast(cython.int, (min_point[2] / voxel_size) - 1)
    z_end = cython.cast(cython.int, (max_point[2] / voxel_size) + 1)

    num_centers = (x_end - x_start) * (y_end - y_start) * (z_end - z_start)
    centers: vector[Tuple[cython.double, cython.double, cython.double]]
    centers.resize(num_centers)

    num_centers = 0
    for x in range(x_start, x_end):
        for y in range(y_start, y_end):
            for z in range(z_start, z_end):
                centers[num_centers] = (
                    _round_to_digits((x + 0.5) * voxel_size, 6),
                    _round_to_digits((y + 0.5) * voxel_size, 6),
                    _round_to_digits((z + 0.5) * voxel_size, 6),
                )
                num_centers += 1

    return centers


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def _flood_fill_matrix_2d(
    matrix: bool_C[:, :],
    start: Tuple[cython.int, cython.int],
    fill_with: bool_C,
    shape: Tuple[cython.int, cython.int],
) -> bool_C[:, :]:
    """Apply a flood fill algorithm to a 2D boolean matrix from a given starting point."""

    dx: cython.int[4] = [0, 0, -1, 1]
    dy: cython.int[4] = [-1, 1, 0, 0]
    nx: cython.int
    ny: cython.int
    x: cython.int
    y: cython.int
    sx: cython.int = shape[0]
    sy: cython.int = shape[1]

    old_value: cython.int = matrix[start[0], start[1]]

    if old_value == fill_with:
        return matrix

    fill_stack: stack[Tuple[cython.int, cython.int]]
    fill_stack.push((start[0], start[1]))

    while not fill_stack.empty():
        x, y = fill_stack.top()
        fill_stack.pop()
        matrix[x, y] = fill_with

        for i in range(4):
            nx, ny = x + dx[i], y + dy[i]

            if 0 <= nx < sx and 0 <= ny < sy and matrix[nx, ny] == old_value:
                fill_stack.push((nx, ny))

    return matrix


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def _flood_fill_matrix_3d(
    matrix: bool_C[:, :, :],
    start: Tuple[cython.int, cython.int, cython.int],
    fill_with: bool_C,
    shape: Tuple[cython.int, cython.int, cython.int],
) -> bool_C[:, :, :]:
    """Apply a flood fill algorithm to a 3D boolean matrix from a given starting point."""

    dx: cython.int[6] = [0, 0, -1, 1, 0, 0]
    dy: cython.int[6] = [-1, 1, 0, 0, 0, 0]
    dz: cython.int[6] = [0, 0, 0, 0, -1, 1]
    nx: cython.int
    ny: cython.int
    nz: cython.int
    x: cython.int
    y: cython.int
    z: cython.int
    sx: cython.int = shape[0]
    sy: cython.int = shape[1]
    sz: cython.int = shape[2]

    old_value: cython.int = matrix[start[0], start[1], start[2]]

    if old_value == fill_with:
        return matrix

    fill_stack: stack[Tuple[cython.int, cython.int, cython.int]]
    fill_stack.push((start[0], start[1], start[2]))

    while not fill_stack.empty():
        x, y, z = fill_stack.top()
        fill_stack.pop()
        matrix[x, y, z] = fill_with

        for i in range(6):
            nx, ny, nz = x + dx[i], y + dy[i], z + dz[i]

            if 0 <= nx < sx and 0 <= ny < sy and 0 <= nz < sz and matrix[nx, ny, nz] == old_value:
                fill_stack.push((nx, ny, nz))

    return matrix


@cython.cfunc
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def _line_segment_intersects_pixel(
    x1: cython.double,
    y1: cython.double,
    x2: cython.double,
    y2: cython.double,
    pixel_center_x: cython.double,
    pixel_center_y: cython.double,
    pixel_size: cython.double,
) -> bool_C:
    """Check if a line segment intersects with a pixel defined by its center and size."""

    # Determine the coordinates of lower-left and upper-right of rectangle
    xmin, xmax = pixel_center_x - pixel_size / 2, pixel_center_x + pixel_size / 2
    ymin, ymax = pixel_center_y - pixel_size / 2, pixel_center_y + pixel_size / 2

    # Compute the line equation for a point

    line_eq1 = _round_to_digits((y2 - y1) * xmin + (x1 - x2) * ymin + (x2 * y1 - x1 * y2), 6)
    line_eq2 = _round_to_digits((y2 - y1) * xmin + (x1 - x2) * ymax + (x2 * y1 - x1 * y2), 6)
    line_eq3 = _round_to_digits((y2 - y1) * xmax + (x1 - x2) * ymin + (x2 * y1 - x1 * y2), 6)
    line_eq4 = _round_to_digits((y2 - y1) * xmax + (x1 - x2) * ymax + (x2 * y1 - x1 * y2), 6)

    # Check if all corners are on the same side of the line
    miss: bool_C = (
        (line_eq1 >= 0 and line_eq2 >= 0 and line_eq3 >= 0 and line_eq4 >= 0)
        or (line_eq1 < 0 and line_eq2 < 0 and line_eq3 < 0 and line_eq4 < 0)
    ) and (
        (line_eq1 > 0 and line_eq2 > 0 and line_eq3 > 0 and line_eq4 > 0)
        or (line_eq1 <= 0 and line_eq2 <= 0 and line_eq3 <= 0 and line_eq4 <= 0)
    )

    # Does it miss based on the shadow intersection test?
    shadow_miss: bool_C = (
        (x1 > xmax and x2 > xmax) or (x1 < xmin and x2 < xmin) or (y1 > ymax and y2 > ymax) or (y1 < ymin and y2 < ymin)
    )

    # A hit is if it doesn't miss on both tests!
    return not (miss or shadow_miss)


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _line_segments_to_pixels(
    line_segments: vector[Tuple[Tuple[cython.double, cython.double], Tuple[cython.double, cython.double]]],
    pixel_size: cython.double,
) -> vector[Tuple[cython.double, cython.double]]:
    """Convert line segments to a list of pixel centers they intersect with."""

    pixel_centers: vector[Tuple[cython.double, cython.double]]

    for i in range(line_segments.size()):
        x1: cython.double = line_segments[i][0][0]
        y1: cython.double = line_segments[i][0][1]
        x2: cython.double = line_segments[i][1][0]
        y2: cython.double = line_segments[i][1][1]

        # Calculate the bounding box of the line segment
        xmin = min(x1, x2)
        ymin = min(y1, y2)
        xmax = max(x1, x2)
        ymax = max(y1, y2)

        # Calculate the indices of the box that intersect with the bounding box of the line segment
        x_start = cython.cast(cython.int, (xmin / pixel_size) - 2)
        x_end = cython.cast(cython.int, (xmax / pixel_size) + 2)
        y_start = cython.cast(cython.int, (ymin / pixel_size) - 2)
        y_end = cython.cast(cython.int, (ymax / pixel_size) + 2)

        # Create a list of the centers of all the intersecting voxels
        for x in range(x_start, x_end):
            for y in range(y_start, y_end):
                x_coord: cython.double = (cython.cast(cython.double, x) + 0.5) * pixel_size
                y_coord: cython.double = (cython.cast(cython.double, y) + 0.5) * pixel_size
                center: Tuple[cython.double, cython.double] = (
                    _round_to_digits(x_coord, 6),
                    _round_to_digits(y_coord, 6),
                )

                if _line_segment_intersects_pixel(x1, y1, x2, y2, center[0], center[1], pixel_size):
                    pixel_centers.push_back(center)

    return pixel_centers


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _triangle_2d_to_pixels(
    triangle_2d: Tuple[
        Tuple[cython.double, cython.double], Tuple[cython.double, cython.double], Tuple[cython.double, cython.double]
    ],
    pixel_size: cython.double,
) -> vector[Tuple[cython.double, cython.double]]:
    """Convert line segments to a list of pixel centers they intersect with."""

    line_segments: vector[Tuple[Tuple[cython.double, cython.double], Tuple[cython.double, cython.double]]]
    line_segments.push_back((triangle_2d[0], triangle_2d[1]))
    line_segments.push_back((triangle_2d[1], triangle_2d[2]))
    line_segments.push_back((triangle_2d[2], triangle_2d[0]))

    return _line_segments_to_pixels(line_segments, pixel_size)


@cython.cfunc
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def _triangles_to_voxel_matrix(
    triangles: vector[
        Tuple[
            Tuple[cython.double, cython.double, cython.double],
            Tuple[cython.double, cython.double, cython.double],
            Tuple[cython.double, cython.double, cython.double],
        ]
    ],
    n_triangles: cython.int,
    voxel_size: cython.double,
    matrix: bool_C[:, :, :],
    matrix_origin_center: Tuple[cython.double, cython.double, cython.double],
) -> bool_C[:, :, :]:
    """Convert 3D triangles to a voxel matrix representation using intersection tests."""

    # Check interface voxels
    for i in range(n_triangles):
        x_abscissa = _round_to_digits(triangles[i][0][0], 6)
        y_abscissa = _round_to_digits(triangles[i][0][1], 6)
        z_abscissa = _round_to_digits(triangles[i][0][2], 6)

        # Check if two points of the triangle are equal
        if _check_triangle_equal_point(triangles[i]):
            pass

        # Check if the triangle is in the YZ plane at the interface between voxels
        elif (
            _round_to_digits(triangles[i][0][0], 6)
            == _round_to_digits(triangles[i][1][0], 6)
            == _round_to_digits(triangles[i][2][0], 6)
        ) and _is_integer(_round_to_digits(x_abscissa / voxel_size, 6)):
            # Define the 3D triangle in 2D
            p0: Tuple[cython.double, cython.double] = (triangles[i][0][1], triangles[i][0][2])
            p1: Tuple[cython.double, cython.double] = (triangles[i][1][1], triangles[i][1][2])
            p2: Tuple[cython.double, cython.double] = (triangles[i][2][1], triangles[i][2][2])

            triangle_2d: Tuple[
                Tuple[cython.double, cython.double],
                Tuple[cython.double, cython.double],
                Tuple[cython.double, cython.double],
            ] = (p0, p1, p2)

            # Compute intersecting pixels
            pixels = _triangle_2d_to_pixels(triangle_2d, voxel_size)

            min_center: Tuple[cython.double, cython.double] = _get_min_pixel_grid_center(pixels)
            max_center: Tuple[cython.double, cython.double] = _get_max_pixel_grid_center(pixels)

            dim_x = cython.cast(cython.int, math_c.round((max_center[0] - min_center[0]) / voxel_size + 1))
            dim_y = cython.cast(cython.int, math_c.round((max_center[1] - min_center[1]) / voxel_size + 1))

            pixel_matrix = _pixel_centers_to_outer_filled_pixel_matrix(pixels, voxel_size, (dim_x, dim_y), min_center)

            # Put the corresponding voxels at True, using the indices
            ix1 = cython.cast(
                cython.int, _round_to_digits((x_abscissa - matrix_origin_center[0] - voxel_size / 2) / voxel_size, 6)
            )
            ix2 = cython.cast(
                cython.int, _round_to_digits((x_abscissa - matrix_origin_center[0] + voxel_size / 2) / voxel_size, 6)
            )

            dx: cython.int
            dy: cython.int
            for dx in range(dim_x):
                for dy in range(dim_y):
                    if pixel_matrix[dx + 1, dy + 1]:
                        iy = cython.cast(
                            cython.int,
                            _round_to_digits(
                                ((min_center[0] + dx * voxel_size) - matrix_origin_center[1]) / voxel_size, 6
                            ),
                        )
                        iz = cython.cast(
                            cython.int,
                            _round_to_digits(
                                ((min_center[1] + dy * voxel_size) - matrix_origin_center[2]) / voxel_size, 6
                            ),
                        )

                        matrix[ix1, iy, iz] = True
                        matrix[ix2, iy, iz] = True

        # Check if the triangle is in the XZ plane at the interface between voxels
        elif (
            _round_to_digits(triangles[i][0][1], 6)
            == _round_to_digits(triangles[i][1][1], 6)
            == _round_to_digits(triangles[i][2][1], 6)
        ) and _is_integer(_round_to_digits(y_abscissa / voxel_size, 6)):
            # Define the 3D triangle in 2D
            p0: Tuple[cython.double, cython.double] = (triangles[i][0][0], triangles[i][0][2])
            p1: Tuple[cython.double, cython.double] = (triangles[i][1][0], triangles[i][1][2])
            p2: Tuple[cython.double, cython.double] = (triangles[i][2][0], triangles[i][2][2])

            triangle_2d: Tuple[
                Tuple[cython.double, cython.double],
                Tuple[cython.double, cython.double],
                Tuple[cython.double, cython.double],
            ] = (p0, p1, p2)

            # Compute intersecting pixels
            pixels = _triangle_2d_to_pixels(triangle_2d, voxel_size)

            min_center: Tuple[cython.double, cython.double] = _get_min_pixel_grid_center(pixels)
            max_center: Tuple[cython.double, cython.double] = _get_max_pixel_grid_center(pixels)

            dim_x = cython.cast(cython.int, math_c.round((max_center[0] - min_center[0]) / voxel_size + 1))
            dim_y = cython.cast(cython.int, math_c.round((max_center[1] - min_center[1]) / voxel_size + 1))

            pixel_matrix = _pixel_centers_to_outer_filled_pixel_matrix(pixels, voxel_size, (dim_x, dim_y), min_center)

            # Put the corresponding voxels at True, using the indices
            iy1 = cython.cast(
                cython.int, _round_to_digits((y_abscissa - matrix_origin_center[1] - voxel_size / 2) / voxel_size, 6)
            )
            iy2 = cython.cast(
                cython.int, _round_to_digits((y_abscissa - matrix_origin_center[1] + voxel_size / 2) / voxel_size, 6)
            )

            dx: cython.int
            dy: cython.int
            for dx in range(dim_x):
                for dy in range(dim_y):
                    if pixel_matrix[dx + 1, dy + 1]:
                        ix = cython.cast(
                            cython.int,
                            _round_to_digits(
                                ((min_center[0] + dx * voxel_size) - matrix_origin_center[0]) / voxel_size, 6
                            ),
                        )
                        iz = cython.cast(
                            cython.int,
                            _round_to_digits(
                                ((min_center[1] + dy * voxel_size) - matrix_origin_center[2]) / voxel_size, 6
                            ),
                        )

                        matrix[ix, iy1, iz] = True
                        matrix[ix, iy2, iz] = True

        # Check if the triangle is in the XY plane at the interface between voxels
        elif (
            _round_to_digits(triangles[i][0][2], 6)
            == _round_to_digits(triangles[i][1][2], 6)
            == _round_to_digits(triangles[i][2][2], 6)
        ) and _is_integer(_round_to_digits(z_abscissa / voxel_size, 6)):
            # Define the 3D triangle in 2D
            p0: Tuple[cython.double, cython.double] = (triangles[i][0][0], triangles[i][0][1])
            p1: Tuple[cython.double, cython.double] = (triangles[i][1][0], triangles[i][1][1])
            p2: Tuple[cython.double, cython.double] = (triangles[i][2][0], triangles[i][2][1])

            triangle_2d: Tuple[
                Tuple[cython.double, cython.double],
                Tuple[cython.double, cython.double],
                Tuple[cython.double, cython.double],
            ] = (p0, p1, p2)

            # Compute intersecting pixels
            pixels = _triangle_2d_to_pixels(triangle_2d, voxel_size)

            min_center: Tuple[cython.double, cython.double] = _get_min_pixel_grid_center(pixels)
            max_center: Tuple[cython.double, cython.double] = _get_max_pixel_grid_center(pixels)

            dim_x = cython.cast(cython.int, math_c.round((max_center[0] - min_center[0]) / voxel_size + 1))
            dim_y = cython.cast(cython.int, math_c.round((max_center[1] - min_center[1]) / voxel_size + 1))

            pixel_matrix = _pixel_centers_to_outer_filled_pixel_matrix(pixels, voxel_size, (dim_x, dim_y), min_center)

            # Put the corresponding voxels at True, using the indices
            iz1 = cython.cast(
                cython.int, _round_to_digits((z_abscissa - matrix_origin_center[2] - voxel_size / 2) / voxel_size, 6)
            )
            iz2 = cython.cast(
                cython.int, _round_to_digits((z_abscissa - matrix_origin_center[2] + voxel_size / 2) / voxel_size, 6)
            )

            dx: cython.int
            dy: cython.int
            for dx in range(dim_x):
                for dy in range(dim_y):
                    if pixel_matrix[dx + 1, dy + 1]:
                        ix = cython.cast(
                            cython.int,
                            _round_to_digits(
                                ((min_center[0] + dx * voxel_size) - matrix_origin_center[0]) / voxel_size, 6
                            ),
                        )
                        iy = cython.cast(
                            cython.int,
                            _round_to_digits(
                                ((min_center[1] + dy * voxel_size) - matrix_origin_center[1]) / voxel_size, 6
                            ),
                        )

                        matrix[ix, iy, iz1] = True
                        matrix[ix, iy, iz2] = True

        # Check intersecting voxels
        else:
            triangle: Tuple[
                Tuple[cython.double, cython.double, cython.double],
                Tuple[cython.double, cython.double, cython.double],
                Tuple[cython.double, cython.double, cython.double],
            ] = triangles[i]
            min_point: Tuple[cython.double, cython.double, cython.double]
            max_point: Tuple[cython.double, cython.double, cython.double]
            min_point, max_point = _triangle_min_max_points(triangle)

            voxel_centers = _aabb_intersecting_boxes(min_point, max_point, voxel_size)
            for j in range(voxel_centers.size()):
                ix = cython.cast(
                    cython.int,
                    _round_to_digits((voxel_centers[j][0] - matrix_origin_center[0]) / voxel_size, 6),
                )
                iy = cython.cast(
                    cython.int,
                    _round_to_digits((voxel_centers[j][1] - matrix_origin_center[1]) / voxel_size, 6),
                )
                iz = cython.cast(
                    cython.int,
                    _round_to_digits((voxel_centers[j][2] - matrix_origin_center[2]) / voxel_size, 6),
                )

                if not matrix[ix, iy, iz]:
                    if _triangle_intersects_voxel(
                        triangle,
                        voxel_centers[j],
                        (0.5 * voxel_size, 0.5 * voxel_size, 0.5 * voxel_size),
                    ):
                        matrix[ix, iy, iz] = True

    return matrix


@cython.cfunc
@cython.exceptval(check=False)
def _get_min_pixel_grid_center(
    pixel_centers: vector[Tuple[cython.double, cython.double]]
) -> Tuple[cython.double, cython.double]:
    """Calculate and return the minimum x and y coordinates among a collection of pixel centers."""

    min_x = min_y = math_c.INFINITY
    for i in range(pixel_centers.size()):
        min_x = min(min_x, pixel_centers[i][0])
        min_y = min(min_y, pixel_centers[i][1])

    return min_x, min_y


@cython.cfunc
@cython.exceptval(check=False)
def _get_max_pixel_grid_center(
    pixel_centers: vector[Tuple[cython.double, cython.double]]
) -> Tuple[cython.double, cython.double]:
    """Calculate and return the maximum x and y coordinates among a collection of pixel centers."""

    max_x = max_y = -math_c.INFINITY
    for i in range(pixel_centers.size()):
        max_x = max(max_x, pixel_centers[i][0])
        max_y = max(max_y, pixel_centers[i][1])

    return max_x, max_y


@cython.cfunc
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def _pixel_centers_to_outer_filled_pixel_matrix(
    pixel_centers: vector[Tuple[cython.double, cython.double]],
    pixel_size: cython.double,
    shape: Tuple[cython.int, cython.int],
    min_center: Tuple[cython.double, cython.double],
) -> bool_C[:, :]:
    """
    Convert pixel centers to a boolean matrix where pixels are marked as filled and their outer boundary is identified.
    """

    matrix: bool_C[:, :] = np.zeros((shape[0] + 2, shape[1] + 2), dtype=np.bool_)

    for i in range(pixel_centers.size()):
        ix = cython.cast(cython.int, _round_to_digits((pixel_centers[i][0] - min_center[0]) / pixel_size, 6)) + 1
        iy = cython.cast(cython.int, _round_to_digits((pixel_centers[i][1] - min_center[1]) / pixel_size, 6)) + 1
        matrix[ix, iy] = True

    matrix_outer_filled = _flood_fill_matrix_2d(matrix, (0, 0), True, shape)

    for ix in range(shape[0]):
        for iy in range(shape[1]):
            if not matrix_outer_filled[ix, iy]:
                matrix[ix, iy] = True

    return matrix


@cython.cfunc
@cython.exceptval(check=False)
def _check_triangle_equal_point(
    triangle: Tuple[
        Tuple[cython.double, cython.double, cython.double],
        Tuple[cython.double, cython.double, cython.double],
        Tuple[cython.double, cython.double, cython.double],
    ]
) -> bool_C:
    """Check if any two points of a triangle are equal, implying degeneracy."""

    return (
        (triangle[0][0] == triangle[1][0] and triangle[0][1] == triangle[1][1] and triangle[0][2] == triangle[1][2])
        or (triangle[0][0] == triangle[2][0] and triangle[0][1] == triangle[2][1] and triangle[0][2] == triangle[2][2])
        or (triangle[1][0] == triangle[2][0] and triangle[1][1] == triangle[2][1] and triangle[1][2] == triangle[2][2])
    )


@cython.cfunc
@cython.exceptval(check=False)
def _triangle_min_max_points(
    triangle: Tuple[
        Tuple[cython.double, cython.double, cython.double],
        Tuple[cython.double, cython.double, cython.double],
        Tuple[cython.double, cython.double, cython.double],
    ]
) -> Tuple[Tuple[cython.double, cython.double, cython.double], Tuple[cython.double, cython.double, cython.double]]:
    """Calculate and return the minimum and maximum coordinates of a 3D triangle."""

    min_x: cython.double = math_c.INFINITY
    min_y: cython.double = math_c.INFINITY
    min_z: cython.double = math_c.INFINITY
    max_x: cython.double = -math_c.INFINITY
    max_y: cython.double = -math_c.INFINITY
    max_z: cython.double = -math_c.INFINITY

    min_x = min(min_x, triangle[0][0])
    min_x = min(min_x, triangle[1][0])
    min_x = min(min_x, triangle[2][0])

    min_y = min(min_y, triangle[0][1])
    min_y = min(min_y, triangle[1][1])
    min_y = min(min_y, triangle[2][1])

    min_z = min(min_z, triangle[0][2])
    min_z = min(min_z, triangle[1][2])
    min_z = min(min_z, triangle[2][2])

    max_x = max(max_x, triangle[0][0])
    max_x = max(max_x, triangle[1][0])
    max_x = max(max_x, triangle[2][0])

    max_y = max(max_y, triangle[0][1])
    max_y = max(max_y, triangle[1][1])
    max_y = max(max_y, triangle[2][1])

    max_z = max(max_z, triangle[0][2])
    max_z = max(max_z, triangle[1][2])
    max_z = max(max_z, triangle[2][2])

    return (min_x, min_y, min_z), (max_x, max_y, max_z)


@cython.cfunc
def _triangles_min_max_points(
    triangles: vector[
        Tuple[
            Tuple[cython.double, cython.double, cython.double],
            Tuple[cython.double, cython.double, cython.double],
            Tuple[cython.double, cython.double, cython.double],
        ]
    ]
) -> Tuple[Tuple[cython.double, cython.double, cython.double], Tuple[cython.double, cython.double, cython.double]]:
    """Calculate and return the minimum and maximum coordinates across a collection of 3D triangles."""

    min_x: cython.double = math_c.INFINITY
    min_y: cython.double = math_c.INFINITY
    min_z: cython.double = math_c.INFINITY
    max_x: cython.double = -math_c.INFINITY
    max_y: cython.double = -math_c.INFINITY
    max_z: cython.double = -math_c.INFINITY

    for i in range(triangles.size()):
        min_x = min(min_x, triangles[i][0][0])
        min_x = min(min_x, triangles[i][1][0])
        min_x = min(min_x, triangles[i][2][0])

        min_y = min(min_y, triangles[i][0][1])
        min_y = min(min_y, triangles[i][1][1])
        min_y = min(min_y, triangles[i][2][1])

        min_z = min(min_z, triangles[i][0][2])
        min_z = min(min_z, triangles[i][1][2])
        min_z = min(min_z, triangles[i][2][2])

        max_x = max(max_x, triangles[i][0][0])
        max_x = max(max_x, triangles[i][1][0])
        max_x = max(max_x, triangles[i][2][0])

        max_y = max(max_y, triangles[i][0][1])
        max_y = max(max_y, triangles[i][1][1])
        max_y = max(max_y, triangles[i][2][1])

        max_z = max(max_z, triangles[i][0][2])
        max_z = max(max_z, triangles[i][1][2])
        max_z = max(max_z, triangles[i][2][2])

    return (min_x, min_y, min_z), (max_x, max_y, max_z)
