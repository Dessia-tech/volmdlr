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
from cython.cimports.libcpp.stack import stack
from cython.cimports.libcpp.vector import vector

# CUSTOM PYTHON TYPES

Point = Tuple[float, ...]
Triangle = Tuple[Point, ...]


# TODO: refactor, add docstrings


@cython.cfunc
@cython.cdivision(True)
@cython.exceptval(check=False)
def _round_to_digits(num: cython.double, digits: cython.int) -> cython.double:
    multiplier: cython.double = math_c.pow(10.0, digits)
    return math_c.round(num * multiplier) / multiplier


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def _triangle_intersects_voxel(
    triangle: cython.double[3][3],
    voxel_center: cython.double[3],
    voxel_extents: cython.double[3],
) -> cython.bint:
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
    voxel_extents: cython.double[3],
) -> cython.bint:
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
    min_point: cython.double[3], max_point: cython.double[3], voxel_size: cython.double
) -> vector[Tuple[cython.double, cython.double, cython.double]]:
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


def triangles_to_voxels(triangles: List[Triangle], voxel_size: float) -> Set[Point]:
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

    for triangle in triangles:
        min_point = tuple(min(p[i] for p in triangle) for i in range(3))
        max_point = tuple(max(p[i] for p in triangle) for i in range(3))

        for bbox_center in _aabb_intersecting_boxes(
            [min_point[0], min_point[1], min_point[2]],
            [max_point[0], max_point[1], max_point[2]],
            voxel_size,
        ):
            bbox_center = tuple(bbox_center)
            if bbox_center not in voxel_centers:
                if _triangle_intersects_voxel(
                    [
                        [triangle[0][0], triangle[0][1], triangle[0][2]],
                        [triangle[1][0], triangle[1][1], triangle[1][2]],
                        [triangle[2][0], triangle[2][1], triangle[2][2]],
                    ],
                    [bbox_center[0], bbox_center[1], bbox_center[2]],
                    [0.5 * voxel_size, 0.5 * voxel_size, 0.5 * voxel_size],
                ):
                    voxel_centers.add(bbox_center)

    return voxel_centers


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def _flood_fill_matrix_3d(
    matrix: cython.bint[:, :, :], start: cython.int[3], fill_with: cython.bint, shape: cython.int[3]
) -> cython.bint[:, :, :]:
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


def flood_fill_matrix_3d(
    matrix: np.ndarray[np.bool_, np.ndim == 3], start: Tuple[int, int, int], fill_with: bool
) -> np.ndarray[np.bool_, np.ndim == 3]:
    return np.asarray(
        _flood_fill_matrix_3d(
            matrix.astype(np.int32),
            [start[0], start[1], start[2]],
            fill_with,
            [matrix.shape[0], matrix.shape[1], matrix.shape[2]],
        ),
        dtype=np.bool_,
    )


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def _flood_fill_matrix_2d(
    matrix: cython.bint[:, :], start: cython.int[2], fill_with: cython.bint, shape: cython.int[2]
) -> cython.bint[:, :]:
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


def flood_fill_matrix_2d(
    matrix: np.ndarray[np.bool_, np.ndim == 2], start: Tuple[int, int], fill_with: bool
) -> np.ndarray[np.bool_, np.ndim == 2]:
    return np.asarray(
        _flood_fill_matrix_2d(
            matrix.astype(np.int32),
            [start[0], start[1]],
            fill_with,
            [matrix.shape[0], matrix.shape[1]],
        ),
        dtype=np.bool_,
    )


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
) -> cython.bint:
    # Determine the coordinates of lower-left and upper-right of rectangle
    xmin, xmax = pixel_center_x - pixel_size / 2, pixel_center_x + pixel_size / 2
    ymin, ymax = pixel_center_y - pixel_size / 2, pixel_center_y + pixel_size / 2

    # Compute the line equation for a point

    line_eq1 = (y2 - y1) * xmin + (x1 - x2) * ymin + (x2 * y1 - x1 * y2)
    line_eq2 = (y2 - y1) * xmin + (x1 - x2) * ymax + (x2 * y1 - x1 * y2)
    line_eq3 = (y2 - y1) * xmax + (x1 - x2) * ymin + (x2 * y1 - x1 * y2)
    line_eq4 = (y2 - y1) * xmax + (x1 - x2) * ymax + (x2 * y1 - x1 * y2)

    # Check if all corners are on the same side of the line
    miss: cython.bint = (line_eq1 <= 0 and line_eq2 <= 0 and line_eq3 <= 0 and line_eq4 <= 0) or (
        line_eq1 >= 0 and line_eq2 >= 0 and line_eq3 >= 0 and line_eq4 >= 0
    )

    # Does it miss based on the shadow intersection test?
    shadow_miss: cython.bint = (
        (x1 > xmax and x2 > xmax) or (x1 < xmin and x2 < xmin) or (y1 > ymax and y2 > ymax) or (y1 < ymin and y2 < ymin)
    )

    # A hit is if it doesn't miss on both tests!
    return not (miss or shadow_miss)


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _line_segments_to_pixels(
    line_segments: vector[Tuple[Tuple[cython.double, cython.double], Tuple[cython.double, cython.double]]], pixel_size: cython.double
) -> vector[Tuple[cython.double, cython.double]]:
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
        x_start = cython.cast(cython.int, (xmin / pixel_size) - 1)
        x_end = cython.cast(cython.int, (xmax / pixel_size) + 1)
        y_start = cython.cast(cython.int, (ymin / pixel_size) - 1)
        y_end = cython.cast(cython.int, (ymax / pixel_size) + 1)

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


def line_segments_to_pixels(
    line_segments: List[Tuple[Tuple[float, float], Tuple[float, float]]], pixel_size: float
) -> Set[Tuple[float, float]]:
    return set(_line_segments_to_pixels(line_segments, pixel_size))


def triangles_to_voxel_matrix(
    triangles: List[Tuple[Tuple[float, float, float]], Tuple[float, float, float], Tuple[float, float, float]],
    voxel_size: float,
) -> Tuple[np.ndarray[np.bool_, np.ndim == 3], Tuple[float, float, float]]:
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
    matrix = np.asarray(_triangles_to_voxel_matrix(triangles, voxel_size, matrix, matrix_origin_center), dtype=np.bool_)

    return matrix, matrix_origin_center


@cython.cfunc
@cython.cdivision(True)
def _triangles_to_voxel_matrix(
    triangles: vector[
        Tuple[
            Tuple[cython.double, cython.double, cython.double],
            Tuple[cython.double, cython.double, cython.double],
            Tuple[cython.double, cython.double, cython.double],
        ]
    ],
    voxel_size: cython.double,
    matrix: cython.bint[:, :, :],
    matrix_origin_center: Tuple[cython.double, cython.double, cython.double],
) -> cython.bint[:, :, :]:
    # Check interface voxels
    for i in range(triangles.size()):
        # Check if the triangle is in the YZ plane
        if (
            _round_to_digits(triangles[i][0][0], 6)
            == _round_to_digits(triangles[i][1][0], 6)
            == _round_to_digits(triangles[i][2][0], 6)
        ):
            # Check if this plane is defined is at the interface between voxels
            abscissa = _round_to_digits(triangles[i][0][0], 6)
            if _is_integer(_round_to_digits(abscissa / voxel_size, 6)):
                # Define the 3D triangle in 2D
                p0: Tuple[cython.double, cython.double] = (triangles[i][0][1], triangles[i][0][2])
                p1: Tuple[cython.double, cython.double] = (triangles[i][1][1], triangles[i][1][2])
                p2: Tuple[cython.double, cython.double] = (triangles[i][2][1], triangles[i][2][2])

                line_segments: vector[
                    Tuple[Tuple[cython.double, cython.double], Tuple[cython.double, cython.double]]
                ]
                line_segments.push_back((p0, p1))
                line_segments.push_back((p1, p2))
                line_segments.push_back((p2, p0))

                pixels = _line_segments_to_pixels(line_segments, voxel_size)

                pass

        # Check if the triangle is in the XZ plane

        # Check if the triangle is in the XY plane

        # Check intersecting voxels
        else:
            # Compute intersectings voxel
            # for these which are false
            # triangle intersects voxel
            pass

    return matrix


@cython.cfunc
@cython.exceptval(check=False)
def _is_integer(value: cython.double) -> cython.bint:
    return cython.cast(cython.int, value) == value


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
