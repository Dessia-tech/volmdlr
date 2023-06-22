# cython: language_level=3

from cpython cimport bool
from libc.stdlib cimport malloc, free
from libc.math cimport fabs, sqrt
from cython.view cimport array as cvarray
from cython.operator cimport dereference as deref
from typing import List, Tuple

cdef struct Point:
    double x
    double y
    double z

cdef struct Triangle:
    Point[3] points

cdef struct VoxelCenter:
    double x
    double y
    double z

cdef struct VoxelExtents:
    double x
    double y
    double z

cdef bint triangle_voxel_intersection_c(Triangle triangle, VoxelCenter voxel_center, VoxelExtents voxel_extents):
    # Method ported from https://gist.github.com/zvonicek/fe73ba9903f49d57314cf7e8e0f05dcf

    cdef Point v0, v1, v2
    cdef Point f0, f1, f2
    cdef Point box_center
    cdef Point box_extents
    cdef Point plane_normal
    cdef double plane_distance, r
    cdef Point a00, a01, a02, a10, a11, a12, a20, a21, a22
    cdef double p0, p1, p2

    # Translate triangle as conceptually moving AABB to origin
    v0.x = triangle.points[0].x - voxel_center.x
    v0.y = triangle.points[0].y - voxel_center.y
    v0.z = triangle.points[0].z - voxel_center.z

    v1.x = triangle.points[1].x - voxel_center.x
    v1.y = triangle.points[1].y - voxel_center.y
    v1.z = triangle.points[1].z - voxel_center.z

    v2.x = triangle.points[2].x - voxel_center.x
    v2.y = triangle.points[2].y - voxel_center.y
    v2.z = triangle.points[2].z - voxel_center.z

    # Compute edge vectors for triangle
    f0.x = triangle.points[1].x - triangle.points[0].x
    f0.y = triangle.points[1].y - triangle.points[0].y
    f0.z = triangle.points[1].z - triangle.points[0].z

    f1.x = triangle.points[2].x - triangle.points[1].x
    f1.y = triangle.points[2].y - triangle.points[1].y
    f1.z = triangle.points[2].z - triangle.points[1].z

    f2.x = triangle.points[0].x - triangle.points[2].x
    f2.y = triangle.points[0].y - triangle.points[2].y
    f2.z = triangle.points[0].z - triangle.points[2].z

    # REGION TEST THE THREE AXES CORRESPONDING TO THE FACE NORMALS OF AABB B (CATEGORY 1)

    # Exit if...
    # ... [-extents.X, extents.X] and [min(v0.X,v1.X,v2.X), max(v0.X,v1.X,v2.X)] do not overlap
    if max(v0.x, v1.x, v2.x) < -voxel_extents.x or min(v0.x, v1.x, v2.x) > voxel_extents.x:
        return False

    # ... [-extents.Y, extents.Y] and [min(v0.Y,v1.Y,v2.Y), max(v0.Y,v1.Y,v2.Y)] do not overlap
    if max(v0.y, v1.y, v2.y) < -voxel_extents.y or min(v0.y, v1.y, v2.y) > voxel_extents.y:
        return False

    # ... [-extents.Z, extents.Z] and [min(v0.Z,v1.Z,v2.Z), max(v0.Z,v1.Z,v2.Z)] do not overlap
    if max(v0.z, v1.z, v2.z) < -voxel_extents.z or min(v0.z, v1.z, v2.z) > voxel_extents.z:
        return False

    # ENDREGION

    # REGION TEST SEPARATING AXIS CORRESPONDING TO TRIANGLE FACE NORMAL (CATEGORY 2)

    plane_normal.x = f0.y * f1.z - f0.z * f1.y
    plane_normal.y = f0.z * f1.x - f0.x * f1.z
    plane_normal.z = f0.x * f1.y - f0.y * f1.x

    plane_distance = fabs(plane_normal.x * v0.x + plane_normal.y * v0.y + plane_normal.z * v0.z)

    # Compute the projection interval radius of b onto L(t) = b.c + t * p.n
    r = (
        voxel_extents.x * fabs(plane_normal.x)
        + voxel_extents.y * fabs(plane_normal.y)
        + voxel_extents.z * fabs(plane_normal.z)
    )

    # Intersection occurs when plane distance falls within [-r,+r] interval
    if plane_distance > r:
        return False

    # ENDREGION

    # REGION TEST AXES a00..a22 (CATEGORY 3)

    # Test axis a00
    a00.x = 0
    a00.y = -f0.z
    a00.z = f0.y
    p0 = v0.x * a00.x + v0.y * a00.y + v0.z * a00.z
    p1 = v1.x * a00.x + v1.y * a00.y + v1.z * a00.z
    p2 = v2.x * a00.x + v2.y * a00.y + v2.z * a00.z
    r = voxel_extents.y * fabs(f0.z) + voxel_extents.z * fabs(f0.y)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a01
    a01.x = 0
    a01.y = -f1.z
    a01.z = f1.y
    p0 = v0.x * a01.x + v0.y * a01.y + v0.z * a01.z
    p1 = v1.x * a01.x + v1.y * a01.y + v1.z * a01.z
    p2 = v2.x * a01.x + v2.y * a01.y + v2.z * a01.z
    r = voxel_extents.y * fabs(f1.z) + voxel_extents.z * fabs(f1.y)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a02
    a02.x = 0
    a02.y = -f2.z
    a02.z = f2.y
    p0 = v0.x * a02.x + v0.y * a02.y + v0.z * a02.z
    p1 = v1.x * a02.x + v1.y * a02.y + v1.z * a02.z
    p2 = v2.x * a02.x + v2.y * a02.y + v2.z * a02.z
    r = voxel_extents.y * fabs(f2.z) + voxel_extents.z * fabs(f2.y)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a10
    a10.x = f0.z
    a10.y = 0
    a10.z = -f0.x
    p0 = v0.x * a10.x + v0.y * a10.y + v0.z * a10.z
    p1 = v1.x * a10.x + v1.y * a10.y + v1.z * a10.z
    p2 = v2.x * a10.x + v2.y * a10.y + v2.z * a10.z
    r = voxel_extents.x * fabs(f0.z) + voxel_extents.z * fabs(f0.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a11
    a11.x = f1.z
    a11.y = 0
    a11.z = -f1.x
    p0 = v0.x * a11.x + v0.y * a11.y + v0.z * a11.z
    p1 = v1.x * a11.x + v1.y * a11.y + v1.z * a11.z
    p2 = v2.x * a11.x + v2.y * a11.y + v2.z * a11.z
    r = voxel_extents.x * fabs(f1.z) + voxel_extents.z * fabs(f1.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a12
    a11.x = f2.z
    a11.y = 0
    a11.z = -f2.x
    p0 = v0.x * a11.x + v0.y * a11.y + v0.z * a11.z
    p1 = v1.x * a11.x + v1.y * a11.y + v1.z * a11.z
    p2 = v2.x * a11.x + v2.y * a11.y + v2.z * a11.z
    r = voxel_extents.x * fabs(f2.z) + voxel_extents.z * fabs(f2.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a20
    a20.x = -f0.y
    a20.y = f0.x
    a20.z = 0
    p0 = v0.x * a20.x + v0.y * a20.y + v0.z * a20.z
    p1 = v1.x * a20.x + v1.y * a20.y + v1.z * a20.z
    p2 = v2.x * a20.x + v2.y * a20.y + v2.z * a20.z
    r = voxel_extents.x * fabs(f0.y) + voxel_extents.y * fabs(f0.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a21
    a21.x = -f1.y
    a21.y = f1.x
    a21.z = 0
    p0 = v0.x * a21.x + v0.y * a21.y + v0.z * a21.z
    p1 = v1.x * a21.x + v1.y * a21.y + v1.z * a21.z
    p2 = v2.x * a21.x + v2.y * a21.y + v2.z * a21.z
    r = voxel_extents.x * fabs(f1.y) + voxel_extents.y * fabs(f1.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # Test axis a22
    a22.x = -f2.y
    a22.y = f2.x
    a22.z = 0
    p0 = v0.x * a22.x + v0.y * a22.y + v0.z * a22.z
    p1 = v1.x * a22.x + v1.y * a22.y + v1.z * a22.z
    p2 = v2.x * a22.x + v2.y * a22.y + v2.z * a22.z
    r = voxel_extents.x * fabs(f2.y) + voxel_extents.y * fabs(f2.x)
    if max(-max(p0, p1, p2), min(p0, p1, p2)) > r:
        return False

    # ENDREGION

    return True


def triangle_intersects_voxel(
    triangle: Tuple[Tuple[float, ...]], voxel_center: Tuple[float, ...], voxel_extents: List[float]
) -> bool:
    """
    Helper method to compute if there is an intersection between a 3D triangle and a voxel.
    This method uses the "Separating Axis Theorem".

    :param triangle: The triangle to check if it intersects with the voxel.
    :type: triangle: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
    :param voxel_center: The center point of the voxel.
    :type voxel_center: tuple[float, float, float]
    :param voxel_extents: The extents of the voxel in each direction (half-size of the voxel size).
    :type voxel_extents: list[float, float, float]

    :return: True if there is an intersection, False otherwise.
    :rtype: bool
    """
    cdef Triangle tri
    cdef VoxelCenter center
    cdef VoxelExtents extents

    tri.points[0].x = triangle[0][0]
    tri.points[0].y = triangle[0][1]
    tri.points[0].z = triangle[0][2]

    tri.points[1].x = triangle[1][0]
    tri.points[1].y = triangle[1][1]
    tri.points[1].z = triangle[1][2]

    tri.points[2].x = triangle[2][0]
    tri.points[2].y = triangle[2][1]
    tri.points[2].z = triangle[2][2]

    center.x = voxel_center[0]
    center.y = voxel_center[1]
    center.z = voxel_center[2]

    extents.x = voxel_extents[0]
    extents.y = voxel_extents[1]
    extents.z = voxel_extents[2]

    return triangle_voxel_intersection_c(tri, center, extents)
