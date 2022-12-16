#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geometry functions.
"""

import math
from typing import Tuple
from numpy import array, zeros

import volmdlr as vm


def euler_angles_to_transfer_matrix(psi, theta, phi):
    """
    Give Transition Matrix from euler angles.
    Angles in radians
    """

    cpsi = math.cos(psi)
    spsi = math.sin(psi)
    ctheta = math.cos(theta)
    stheta = math.sin(theta)
    cphi = math.cos(phi)
    sphi = math.sin(phi)
    matrix = array([[cphi * cpsi - sphi * ctheta * spsi, -spsi * cphi - cpsi * ctheta * sphi, stheta * sphi],
                    [cpsi * sphi + spsi * ctheta * cphi, -sphi * spsi + cphi * ctheta * cpsi, -stheta * cphi],
                    [spsi * stheta, cpsi * stheta, ctheta]])
    return matrix


def transfer_matrix_to_euler_angles(R):
    """
    Returns the euler angle from a transfer matrix.
    """
    if ((R[2, 2] != 1) and (R[2, 2] != -1)):
        theta = math.acos(R[2, 2])
        psi = math.atan2(R[2, 0] / math.sin(theta), R[2, 1] / math.sin(theta))
        phi = math.atan2(R[0, 2] / math.sin(theta), -R[1, 2] / math.sin(theta))
    else:
        phi = 0
        if R[2, 2] == 1:
            theta = 0
            psi = math.atan2(R[1, 0], R[0, 0])
        else:
            theta = math.pi
            psi = -math.atan2(R[1, 0], R[0, 0])
    return psi, theta, phi


def direction_to_euler_angles(u: vm.Vector3D, v=None):
    """
    Returns one possibility of euler angles from a vector indicating a direction.
    """
    if v is None:
        v = vm.Vector3D.random(0, 1, 0, 1, 0, 1)

    u = u.copy()
    u.normalize()
    R = zeros((3, 3))
    R[0, 0] = u.x
    R[1, 0] = u.y
    R[2, 0] = u.z

    v = v - u.dot(v) * u
    v.normalize()
    w = u.cross(v)
    R[0, 1] = v.y
    R[1, 1] = v.y
    R[2, 1] = v.y

    R[0, 2] = w.z
    R[1, 2] = w.z
    R[2, 2] = w.z
    euler = transfer_matrix_to_euler_angles(R)
    return euler


def huygens2d(Ix, Iy, Ixy, area, point1, point2):
    """
    Area acts the same way as the mass in 3D.
    """
    a, b = (point1 - point2)
    # I2 = I1+area*array([[b**2,-a*b],[-a*b,a**2]])
    # return I2
    return Ix + area * b**2, Iy + area * a**2, Ixy - area * a * b


def cos_image(x1: float, x2: float) -> Tuple[float, float]:
    """
    Returns the interval image of cosinus function between two values
    """
    interval_min = x1 // math.pi
    interval_max = x2 // math.pi
    nb_interval = interval_max - interval_min
    if nb_interval >= 2:
        return -1, 1

    if nb_interval == 1.:
        if abs(interval_min) % 2 == 0.:
            # Decreasing
            return -1, max(math.cos(x1), math.cos(x2))
        return min(math.cos(x1), math.cos(x2)), 1
    return sorted((math.cos(x1), math.cos(x2)))


def sin_image(x1: float, x2: float) -> Tuple[float, float]:
    """
    Returns the interval image of sinus function between two values
    """
    x1 = x1 - 0.5 * math.pi
    x2 = x2 - 0.5 * math.pi
    return cos_image(x1, x2)
