# distutils: language = c++
cdef extern from *:
    """
    #include <cmath>
    #include <algorithm>

    struct Vector3 {
        float x, y, z;
        Vector3(float x, float y, float z) : x(x), y(y), z(z) {}
        Vector3 operator-(const Vector3& other) const {
            return Vector3(x - other.x, y - other.y, z - other.z);
        }
    };

    struct Bounds {
        Vector3 center, extents;
        Bounds(Vector3 center, Vector3 extents) : center(center), extents(extents) {}
    };

    inline float Max(float a, float b, float c) {
        return std::max({a, b, c});
    }

    inline float Min(float a, float b, float c) {
        return std::min({a, b, c});
    }

    inline Vector3 Cross(const Vector3& v1, const Vector3& v2) {
        return Vector3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }

    inline float Dot(const Vector3& v1, const Vector3& v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    inline float Abs(float x) {
        return std::abs(x);
    }

    bool BoxIntersectsTriangle(const Bounds& bounds, const Vector3& triangle0, const Vector3& triangle1, const Vector3& triangle2) {
        // The implementation goes here
    }
    """


    cdef cppclass Vector3:
        Vector3(float, float, float)
        Vector3 operator-(Vector3)
        float x
        float y
        float z

    cdef cppclass Bounds:
        Bounds(Vector3, Vector3)
        Vector3 center
        Vector3 extents

    float Max(float, float, float)
    float Min(float, float, float)
    Vector3 Cross(Vector3, Vector3)
    float Dot(Vector3, Vector3)
    float Abs(float)
    bool BoxIntersectsTriangle(Bounds, Vector3, Vector3, Vector3)

cpdef bool box_intersects_triangle((float,float,float) voxel_center, (float,float,float) voxel_extents, ((float,float,float),(float,float,float),(float,float,float)) triangle):
    cdef:
        Vector3 center = Vector3(voxel_center[0], voxel_center[1], voxel_center[2])
        Vector3 extents = Vector3(voxel_extents[0], voxel_extents[1], voxel_extents[2])
        Bounds b = Bounds(center, extents)
        Vector3 t0 = Vector3(triangle[0][0], triangle[0][1], triangle[0][2])
        Vector3 t1 = Vector3(triangle[1][0], triangle[1][1], triangle[1][2])
        Vector3 t2 = Vector3(triangle[2][0], triangle[2][1], triangle[2][2])
    return BoxIntersectsTriangle(b, t0, t1, t2)
