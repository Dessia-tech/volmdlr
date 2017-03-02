#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:14:32 2017

@author: steven
"""
from numpy import dot,cross
from volmdlr.primitives2D import Point2D

def PointProjectionPlane(point,plane_origin,plane_normal):
    """
    Plane is defined by (plane_point,plane_normal)
    :returns : coordinates in global coordinate system
    """
    return Point2D(point-dot(point-plane_origin,plane_normal)*plane_normal)

def PointLocalProjectionPlane(point,plane_origin,x1,x2):
    """
    Plane is defined by (plane_point,x1,x2)
    :returns : coordinates in local coordinate system
    """
    x1=x1.vector
    x2=x2.vector
    point=point.vector
    plane_origin=plane_origin.vector
    plane_normal=cross(x1,x2)
    xp=point-dot(point-plane_origin,plane_normal)*plane_normal-plane_origin# projeted point
    return Point2D((dot(xp,x1),dot(xp,x2)))