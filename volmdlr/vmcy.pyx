#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Cython functions

"""


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


#cdef float CVector3DDot(float u1, float u2, float u3,
#                        float v1, float v2, float v3):
#    return u1*v1 + u2*v2 + u3*v3
#
#def Vector3DDot(vector1, vector2):
#    return CVector3DDot(vector1[0], vector1[1], vector1[2],
#                        vector2[0], vector2[1], vector2[2])
#    
#
#cdef (float, float, float) Csub(float u1, float u2, float u3,
#                                float v1, float v2, float v3):
#    return (u1-v1, u2-v2, u3-v3)
#
#def sub(vector1, vector2):
#    return Csub(vector1[0], vector1[1], vector1[2],
#                vector2[0], vector2[1], vector2[2])
    

cdef double CVector3DDot(double u1, double u2, double u3,
                        double v1, double v2, double v3):
    return u1*v1 + u2*v2 + u3*v3

def Vector3DDot(vector1, vector2):
    return CVector3DDot(vector1[0], vector1[1], vector1[2],
                        vector2[0], vector2[1], vector2[2])
    

cdef (double, double, double) Csub(double u1, double u2, double u3,
                                double v1, double v2, double v3):
    return (u1-v1, u2-v2, u3-v3)

def sub(vector1, vector2):
    return Csub(vector1[0], vector1[1], vector1[2],
                vector2[0], vector2[1], vector2[2])