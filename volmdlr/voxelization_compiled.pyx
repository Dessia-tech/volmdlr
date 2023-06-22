# cython: language_level=3
# cython: embedsignature=True
# cython: boundscheck=False
# cython: wraparound=False

cdef float Max(float a, float b, float c):
    return max(max(a, b), c)

cdef float Min(float a, float b, float c):
    return min(min(a, b), c)

cdef bint BoxIntersectsTriangle(tuple triangle, tuple voxel_center, list voxel_extents) except *:
    cdef tuple v0 = tuple(triangle[0][i] - voxel_center[i] for i in range(3))
    cdef tuple v1 = tuple(triangle[1][i] - voxel_center[i] for i in range(3))
    cdef tuple v2 = tuple(triangle[2][i] - voxel_center[i] for i in range(3))

    # Compute edge vectors for triangle
    cdef tuple f0 = tuple(triangle[1][i] - triangle[0][i] for i in range(3))
    cdef tuple f1 = tuple(triangle[2][i] - triangle[1][i] for i in range(3))
    cdef tuple f2 = tuple(triangle[0][i] - triangle[2][i] for i in range(3))

    # Test axes a00..a22 (category 3)
    cdef tuple a00 = (0, -f0[2], f0[1])
    if not test_axis(a00, v0, v1, v2, f0, voxel_extents):
        return False

    cdef tuple a01 = (0, -f1[2], f1[1])
    if not test_axis(a01, v0, v1, v2, f1, voxel_extents):
        return False

    cdef tuple a02 = (0, -f2[2], f2[1])
    if not test_axis(a02, v0, v1, v2, f2, voxel_extents):
        return False

    cdef tuple a10 = (f0[2], 0, -f0[0])
    if not test_axis(a10, v0, v1, v2, f0, voxel_extents):
        return False

    cdef tuple a11 = (f1[2], 0, -f1[0])
    if not test_axis(a11, v0, v1, v2, f1, voxel_extents):
        return False

    cdef tuple a12 = (f2[2], 0, -f2[0])
    if not test_axis(a12, v0, v1, v2, f2, voxel_extents):
        return False

    cdef tuple a20 = (-f0[1], f0[0], 0)
    if not test_axis(a20, v0, v1, v2, f0, voxel_extents):
        return False

    cdef tuple a21 = (-f1[1], f1[0], 0)
    if not test_axis(a21, v0, v1, v2, f1, voxel_extents):
        return False

    cdef tuple a22 = (-f2[1], f2[0], 0)
    if not test_axis(a22, v0, v1, v2, f2, voxel_extents):
        return False

    # Test the three axes corresponding to the face normals of AABB b (category 1)
    if not test_axis_face(v0, v1, v2, voxel_extents):
        return False

    # Test separating axis corresponding to triangle face normal (category 2)
    cdef tuple plane_normal = (f0[1] * f1[2] - f0[2] * f1[1], f0[2] * f1[0] - f0[0] * f1[2], f0[0] * f1[1] - f0[1] * f1[0])
    if not test_axis_triangle(plane_normal, v0, v1, v2, voxel_extents):
        return False

    return True

cdef bint test_axis(tuple axis, tuple v0, tuple v1, tuple v2, tuple f, list voxel_extents) except *:
    cdef float p0 = v0[0] * axis[0] + v0[1] * axis[1] + v0[2] * axis[2]
    cdef float p1 = v1[0] * axis[0] + v1[1] * axis[1] + v1[2] * axis[2]
    cdef float p2 = v2[0] * axis[0] + v2[1] * axis[1] + v2[2] * axis[2]
    cdef float r = voxel_extents[1] * abs(f[2]) + voxel_extents[2] * abs(f[1])
    return max(-Max(p0, p1, p2), Min(p0, p1, p2)) <= r

cdef bint test_axis_face(tuple v0, tuple v1, tuple v2, list voxel_extents) except *:
    if Max(v0[0], v1[0], v2[0]) < -voxel_extents[0] or Min(v0[0], v1[0], v2[0]) > voxel_extents[0]:
        return False
    if Max(v0[1], v1[1], v2[1]) < -voxel_extents[1] or Min(v0[1], v1[1], v2[1]) > voxel_extents[1]:
        return False
    if Max(v0[2], v1[2], v2[2]) < -voxel_extents[2] or Min(v0[2], v1[2], v2[2]) > voxel_extents[2]:
        return False
    return True

cdef bint test_axis_triangle(tuple plane_normal, tuple v0, tuple v1, tuple v2, list voxel_extents) except *:
    cdef float plane_distance = abs(plane_normal[0] * v0[0] + plane_normal[1] * v0[1] + plane_normal[2] * v0[2])
    cdef float r = voxel_extents[0] * abs(plane_normal[0]) + voxel_extents[1] * abs(plane_normal[1]) + voxel_extents[2] * abs(plane_normal[2])
    return plane_distance <= r
