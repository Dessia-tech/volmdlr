#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

import math
from typing import Tuple
from numpy import dot, cross, array, zeros, random
from scipy.linalg import norm


def euler_angles_to_transfer_matrix(psi, theta, phi):
    """
    Give Transition Matrix from euler angles
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


def direction_to_euler_angles(u, v=random.random(3)):
    #    u=npy.array([ux,uy,uz])
    u = u / norm(u)
    R = zeros((3, 3))
    R[:, 0] = u
    v = v - dot(u, v) * u
    v = v / norm(v)
    w = cross(u, v)
    R[:, 1] = v
    R[:, 2] = w
    euler = transfer_matrix_to_euler_angles(R)
    return euler


def huygens2d(Ix, Iy, Ixy, area, point1, point2):
    """
    area acts the same way as the mass in 3D
    """
    a, b = (point1 - point2)
    # I2 = I1+area*array([[b**2,-a*b],[-a*b,a**2]])
    # return I2
    return Ix + area * b**2, Iy + area * a**2, Ixy - area * a * b


def cos_image(x1: float, x2: float) -> Tuple[float, float]:
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
    x1 = x1 - 0.5 * math.pi
    x2 = x2 - 0.5 * math.pi
    return cos_image(x1, x2)
