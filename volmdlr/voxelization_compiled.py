# cython: language_level=3
# distutils: language = c++
"""
Pure python module to define cython function.
This module needs to be compiled!
"""
import cython
import cython.cimports.libc.math as math_c
from cython.cimports.libcpp.vector import vector

from typing import List, Tuple, Set
import numpy as np

# CUSTOM PYTHON TYPES

Point = Tuple[float, ...]
Triangle = Tuple[Point, ...]


@cython.cfunc
@cython.cdivision(True)
def round_to_digits(num: cython.double, digits: cython.int) -> cython.double:
    multiplier: cython.double = math_c.pow(10.0, digits)
    return math_c.round(num * multiplier) / multiplier


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
def triangle_intersects_voxel(
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
    if calculate_axis_values(v0, v1, v2, a00, f0, voxel_extents):
        return False

    # Test axis a01
    a01[0] = 0
    a01[1] = -f1[2]
    a01[2] = f1[1]
    if calculate_axis_values(v0, v1, v2, a01, f1, voxel_extents):
        return False

    # Test axis a02
    a02[0] = 0
    a02[1] = -f2[2]
    a02[2] = f2[1]
    if calculate_axis_values(v0, v1, v2, a02, f2, voxel_extents):
        return False

    # Test axis a10
    a10[0] = f0[2]
    a10[1] = 0
    a10[2] = -f0[0]
    if calculate_axis_values(v0, v1, v2, a10, f0, voxel_extents):
        return False

    # Test axis a11
    a11[0] = f1[2]
    a11[1] = 0
    a11[2] = -f1[0]
    if calculate_axis_values(v0, v1, v2, a11, f1, voxel_extents):
        return False

    # Test axis a12
    a12[0] = f2[2]
    a12[1] = 0
    a12[2] = -f2[0]
    if calculate_axis_values(v0, v1, v2, a12, f2, voxel_extents):
        return False

    # Test axis a20
    a20[0] = -f0[1]
    a20[1] = f0[0]
    a20[2] = 0
    if calculate_axis_values(v0, v1, v2, a20, f0, voxel_extents):
        return False

    # Test axis a21
    a21[0] = -f1[1]
    a21[1] = f1[0]
    a21[2] = 0
    if calculate_axis_values(v0, v1, v2, a21, f1, voxel_extents):
        return False

    # Test axis a22
    a22[0] = -f2[1]
    a22[1] = f2[0]
    a22[2] = 0
    if calculate_axis_values(v0, v1, v2, a22, f2, voxel_extents):
        return False

    # ENDREGION

    return True


@cython.cfunc
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.exceptval(check=False)
def calculate_axis_values(
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
def aabb_intersecting_boxes(
    min_point: cython.double[3], max_point: cython.double[3], voxel_size: cython.double
) -> List[cython.double[3]]:
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
    center: cython.double[3]

    x_start = cython.cast(cython.int, (min_point[0] / voxel_size) - 1)
    x_end = cython.cast(cython.int, (max_point[0] / voxel_size) + 1)
    y_start = cython.cast(cython.int, (min_point[1] / voxel_size) - 1)
    y_end = cython.cast(cython.int, (max_point[1] / voxel_size) + 1)
    z_start = cython.cast(cython.int, (min_point[2] / voxel_size) - 1)
    z_end = cython.cast(cython.int, (max_point[2] / voxel_size) + 1)

    num_centers = (x_end - x_start) * (y_end - y_start) * (z_end - z_start)
    centers: list = num_centers * [center]  # TODO: use cpp vector

    num_centers = 0
    for x in range(x_start, x_end):
        for y in range(y_start, y_end):
            for z in range(z_start, z_end):
                center[0] = round_to_digits((x + 0.5) * voxel_size, 6)
                center[1] = round_to_digits((y + 0.5) * voxel_size, 6)
                center[2] = round_to_digits((z + 0.5) * voxel_size, 6)
                centers[num_centers] = center
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
        # Check if the triangle is at the interface of two voxels, and add them
        # voxel_centers = voxel_centers.union(Voxelization._triangle_interface_voxels(triangle, voxel_size))

        min_point = tuple(min(p[i] for p in triangle) for i in range(3))
        max_point = tuple(max(p[i] for p in triangle) for i in range(3))

        for bbox_center in aabb_intersecting_boxes(
            [min_point[0], min_point[1], min_point[2]],
            [max_point[0], max_point[1], max_point[2]],
            voxel_size,
        ):
            bbox_center = tuple(bbox_center)
            if bbox_center not in voxel_centers:
                if triangle_intersects_voxel(
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
def flood_fill_matrix(matrix: vector[vector[vector[cython.int]]], start: cython.int[0], fill_with: cython.int) -> vector[vector[vector[cython.int]]]:
    dx: cython.int[6] = [0, 0, -1, 1, 0, 0]
    dy: cython.int[6] = [-1, 1, 0, 0, 0, 0]
    dz: cython.int[6] = [0, 0, 0, 0, -1, 1]
    nx: cython.int
    ny: cython.int
    nz: cython.int
    x: cython.int
    y: cython.int
    z: cython.int
    sx: cython.int = matrix.size()
    sy: cython.int = matrix[0].size()
    sz: cython.int = matrix[0][0].size()

    old_value = matrix[start[0]][start[1]][start[2]]

    if old_value == fill_with:
        return matrix

    stack: vector[cython.pint]
    stack.push_back(start)

    while stack.size() > 0:
        x, y, z = stack.pop_back()

        matrix[x][y][z] = fill_with

        for i in range(6):
            nx, ny, nz = x + dx[i], y + dy[i], z + dz[i]

            if (
                0 <= nx < sx
                and 0 <= ny < sy
                and 0 <= nz < sz
                and matrix[nx][ny][nz] == old_value
            ):
                stack.push_back([nx, ny, nz])

    return matrix