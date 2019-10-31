#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Cython functions

"""
# =============================================================================

cdef (double, double) Csub2D(double u1, double u2,
                             double v1, double v2):
    return (u1-v1, u2-v2)

def sub2D(vector1, vector2):
    return Csub2D(vector1[0], vector1[1],
                  vector2[0], vector2[1])

# =============================================================================

cdef (double, double) Cadd2D(double u1, double u2,
                             double v1, double v2,):
    return (u1+v1, u2+v2)

def add2D(vector1, vector2):
    return Cadd2D(vector1[0], vector1[1],
                  vector2[0], vector2[1])

# =============================================================================

cdef (double, double) Cmul2D(double u1, double u2, double value):
    return (u1*value, u2*value)

def mul2D(vector, value):
    return Cmul2D(vector[0], vector[1], value)

# =============================================================================

cdef double CVector2DDot(double u1, double u2,
                         double v1, double v2):
    return u1*v1 + u2*v2

def Vector2DDot(vector1, vector2):
    return CVector2DDot(vector1[0], vector1[1],
                        vector2[0], vector2[1])

# =============================================================================

cdef double CVector2DNorm(double u1, double u2):
    return (u1*u1 + u2*u2)**0.5

def Vector2DNorm(vector):
    return CVector2DNorm(vector[0], vector[1])

# =============================================================================

cdef (double, double, double) Csub3D(double u1, double u2, double u3,
                                     double v1, double v2, double v3):
    return (u1-v1, u2-v2, u3-v3)

def sub3D(vector1, vector2):
    return Csub3D(vector1[0], vector1[1], vector1[2],
                  vector2[0], vector2[1], vector2[2])

# =============================================================================

cdef (double, double, double) Cadd3D(double u1, double u2, double u3,
                                     double v1, double v2, double v3):
    return (u1+v1, u2+v2, u3+v3)

def add3D(vector1, vector2):
    return Cadd3D(vector1[0], vector1[1], vector1[2],
                  vector2[0], vector2[1], vector2[2])

# =============================================================================

cdef (double, double, double) Cmul3D(double u1, double u2, double u3,
                                     double value):
    return (u1*value, u2*value, u3*value)

def mul3D(vector, value):
    return Cmul3D(vector[0], vector[1], vector[2], value)

# =============================================================================

cdef double CVector3DDot(double u1, double u2, double u3,
                         double v1, double v2, double v3):
    return u1*v1 + u2*v2 + u3*v3

def Vector3DDot(vector1, vector2):
    return CVector3DDot(vector1[0], vector1[1], vector1[2],
                        vector2[0], vector2[1], vector2[2])

# =============================================================================

cdef double CVector3DNorm(double u1, double u2, double u3):
    return (u1*u1 + u2*u2 + u3*u3)**0.5

def Vector3DNorm(vector):
    return CVector3DNorm(vector[0], vector[1], vector[2])

# =============================================================================

def PolygonPointBelongs(point, points):

    cdef int i
    cdef int n = len(points)
    cdef bint inside = False
    cdef float x, y, p1x, p1y, p2x, p2y, xints
    x,y=point
    p1x,p1y = points[0]

    for i in range(n+1):
        p2x, p2y = points[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

# =============================================================================

cdef (double, (double, double)) CLineSegment2DPointDistance((double, double) p1, (double, double) p2, (double, double) point):
    cdef double t
    cdef (double, double) u, projection

    u = sub2D(p2, p1)
    t = max(0, min(1, Vector2DDot(sub2D(point, p1), u) / Vector2DNorm(u)**2))
    projection = add2D(p1, mul2D(u, t))
    return Vector2DNorm(sub2D(projection, point)), projection

def LineSegment2DPointDistance(points, point):
    return CLineSegment2DPointDistance(tuple(points[0]), tuple(points[1]), tuple(point))
