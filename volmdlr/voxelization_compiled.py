"""
Pure python module to define cython function.
This module needs to be compiled!
"""
import cython
from typing import Tuple
from cython.cimports.libc.math import fabs, sqrt

# C TYPES DEFINITION


class Vector3D:
    x: cython.double
    y: cython.double
    z: cython.double


class Triangle3D:
    points: Tuple[Vector3D, Vector3D, Vector3D]


@cython.ccall
@cython.boundscheck(False)
def triangle_intersects_voxel_c(
    triangle: Triangle3D,
    voxel_center: Vector3D,
    voxel_extents: Vector3D,
) -> cython.bint:
    v0: Vector3D
    v1: Vector3D
    v2: Vector3D
    f0: Vector3D
    f1: Vector3D
    f2: Vector3D
    box_center: Vector3D
    box_extents: Vector3D
    plane_normal: Vector3D
    plane_distance: cython.double
    r: cython.double
    a00: Vector3D
    a01: Vector3D
    a02: Vector3D
    a10: Vector3D
    a11: Vector3D
    a12: Vector3D
    a20: Vector3D
    a21: Vector3D
    a22: Vector3D
    p0: cython.double
    p1: cython.double
    p2: cython.double

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
    if not calculate_axis_values(v0, v1, v2, a00, f0, voxel_extents):
        return False

    # Test axis a01
    a01.x = 0
    a01.y = -f1.z
    a01.z = f1.y
    if not calculate_axis_values(v0, v1, v2, a01, f1, voxel_extents):
        return False

    # Test axis a02
    a02.x = 0
    a02.y = -f2.z
    a02.z = f2.y
    if not calculate_axis_values(v0, v1, v2, a02, f2, voxel_extents):
        return False

    # Test axis a10
    a10.x = f0.z
    a10.y = 0
    a10.z = -f0.x
    if not calculate_axis_values(v0, v1, v2, a10, f0, voxel_extents):
        return False

    # Test axis a11
    a11.x = f1.z
    a11.y = 0
    a11.z = -f1.x
    if not calculate_axis_values(v0, v1, v2, a11, f1, voxel_extents):
        return False

    # Test axis a12
    a12.x = f2.z
    a12.y = 0
    a12.z = -f2.x
    if not calculate_axis_values(v0, v1, v2, a12, f2, voxel_extents):
        return False

    # Test axis a20
    a20.x = -f0.y
    a20.y = f0.x
    a20.z = 0
    if not calculate_axis_values(v0, v1, v2, a20, f0, voxel_extents):
        return False

    # Test axis a21
    a21.x = -f1.y
    a21.y = f1.x
    a21.z = 0
    if not calculate_axis_values(v0, v1, v2, a21, f1, voxel_extents):
        return False

    # Test axis a22
    a22.x = -f2.y
    a22.y = f2.x
    a22.z = 0
    if not calculate_axis_values(v0, v1, v2, a22, f2, voxel_extents):
        return False

    # ENDREGION

    return True


@cython.ccall
@cython.boundscheck(False)
def calculate_axis_values(
    v0: Vector3D, v1: Vector3D, v2: Vector3D, ax: Vector3D, f: Vector3D, voxel_extents: Vector3D
) -> cython.bint:
    p0 = v0.x * ax.x + v0.y * ax.y + v0.z * ax.z
    p1 = v1.x * ax.x + v1.y * ax.y + v1.z * ax.z
    p2 = v2.x * ax.x + v2.y * ax.y + v2.z * ax.z
    r = voxel_extents.x * fabs(f.z) + voxel_extents.y * fabs(f.x) + voxel_extents.z * fabs(f.y)

    return max(-max(p0, p1, p2), min(p0, p1, p2)) <= r
