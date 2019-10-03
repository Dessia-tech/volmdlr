#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 16:26:30 2017

@author: steven
"""
#from libcpp cimport bool

def PolygonPointBelongs(point,points):
    
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

def Vector3DDot(vector1, vector2):
    
    cdef float u1, u2, u3, v1, v2, v3
    
    u1, u2, u3 = vector1
    v1, v2, v3 = vector2
    return u1*v1 + u2*v2 + u3*v3