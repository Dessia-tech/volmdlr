cdef class Vector:
    cdef public str name

cdef class Vector2D(Vector):
    cdef public double x, y

cdef class Point2D(Vector2D):
    pass

cdef class Vector3D(Vector):
    cdef public double x, y, z

cdef class Point3D(Vector3D):
    pass
