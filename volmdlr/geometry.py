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


def transfer_matrix_to_euler_angles(r_matrix):
    """Returns the Euler angle from a transfer matrix."""
    if (r_matrix[2, 2] != 1) and (r_matrix[2, 2] != -1):
        theta = math.acos(r_matrix[2, 2])
        psi = math.atan2(r_matrix[2, 0] / math.sin(theta), r_matrix[2, 1] / math.sin(theta))
        phi = math.atan2(r_matrix[0, 2] / math.sin(theta), -r_matrix[1, 2] / math.sin(theta))
    else:
        phi = 0
        if r_matrix[2, 2] == 1:
            theta = 0
            psi = math.atan2(r_matrix[1, 0], r_matrix[0, 0])
        else:
            theta = math.pi
            psi = -math.atan2(r_matrix[1, 0], r_matrix[0, 0])
    return psi, theta, phi


def direction_to_euler_angles(u: vm.Vector3D, v=None):
    """
    Returns one possibility of euler angles from a vector indicating a direction.
    """
    if v is None:
        v = vm.Vector3D.random(0, 1, 0, 1, 0, 1)

    u = u.copy()
    u = u.unit_vector()
    matrix_r = zeros((3, 3))
    matrix_r[0, 0] = u.x
    matrix_r[1, 0] = u.y
    matrix_r[2, 0] = u.z

    v = v - u.dot(v) * u
    v = v.unit_vector()
    w = u.cross(v)
    matrix_r[0, 1] = v.y
    matrix_r[1, 1] = v.y
    matrix_r[2, 1] = v.y

    matrix_r[0, 2] = w.z
    matrix_r[1, 2] = w.z
    matrix_r[2, 2] = w.z
    euler = transfer_matrix_to_euler_angles(matrix_r)
    return euler


def huygens2d(Ix, Iy, Ixy, area, point1, point2):
    """
    Area acts the same way as the mass in 3D.
    """
    a, b = point1 - point2
    return Ix + area * b**2, Iy + area * a**2, Ixy - area * a * b


def cos_image(x1: float, x2: float) -> Tuple[float, float]:
    """
    Returns the interval image of cosinus function between two values.

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
    Returns the interval image of sinus function between two values.

    """
    x1 = x1 - 0.5 * math.pi
    x2 = x2 - 0.5 * math.pi
    return cos_image(x1, x2)


def vectors3d_angle(vector1, vector2):
    """
    Computes the angle between two 3 dimensional vectors.

    :param vector1: The fist 3 dimensional vector
    :type vector1: :class:`volmdlr.Vector3D`
    :param vector2: The second 3 dimensional vectors
    :type vector2: :class:`volmdlr.Vector3D`
    :return: The angle between the two vectors
    :rtype: float
    """
    dot_v1v2 = vector1.dot(vector2)
    ratio = dot_v1v2 / (vector1.norm() * vector2.norm())
    ratio = max(min(ratio, 1.0), -1.0)
    theta = math.acos(ratio)

    return theta


def sin_cos_angle(u1, u2):
    """
    Returns an angle between 0 and 2*PI verifying cos(theta)=u1, sin(theta)=u2.

    :param u1: The value of the cosinus of the returned angle
    :type u1: float
    :param u2: The value of the sinus of the returned angle
    :type u2: float
    :return: The angle verifying the two equations
    :rtype: float
    """
    u1 = max(-1.0, min(u1, 1.0))
    u2 = max(-1.0, min(u2, 1.0))

    if u1 > 0:
        if u2 >= 0:
            theta = math.acos(u1)
        else:
            theta = vm.TWO_PI + math.asin(u2)
    else:
        if u2 >= 0:
            theta = math.acos(u1)
        else:
            theta = vm.TWO_PI - math.acos(u1)
    if math.isclose(theta, vm.TWO_PI, abs_tol=1e-9):
        return 0.
    return theta


def clockwise_interior_from_circle3d(start, end, circle):
    """
    Returns the clockwise interior point between start and end on the circle.
    """
    start2d = start.to_2d(plane_origin=circle.frame.origin,
                          x=circle.frame.u, y=circle.frame.v)
    end2d = end.to_2d(plane_origin=circle.frame.origin,
                      x=circle.frame.u, y=circle.frame.v)

    # p1 angle
    u1, u2 = start2d.x / circle.radius, start2d.y / circle.radius
    theta1 = sin_cos_angle(u1, u2)
    # p2 angle
    u3, u4 = end2d.x / circle.radius, end2d.y / circle.radius
    theta2 = sin_cos_angle(u3, u4)

    if theta1 > theta2:
        theta3 = (theta1 + theta2) / 2
    elif theta2 >= theta1:
        theta3 = (theta1 + theta2) / 2 + vm.TWO_PI / 2
    else:
        raise NotImplementedError

    if theta3 > vm.TWO_PI:
        theta3 -= vm.TWO_PI

    interior2d = vm.Point2D(circle.radius * math.cos(theta3),
                            circle.radius * math.sin(theta3))
    interior3d = interior2d.to_3d(plane_origin=circle.frame.origin,
                                  vx=circle.frame.u, vy=circle.frame.v)
    return interior3d


def offset_angle(trigo, angle_start, angle_end):
    """
    Calculates the offset and angle.

    Seems unused
    """
    if trigo:
        offset = angle_start
    else:
        offset = angle_end
    if angle_start > angle_end:
        angle = angle_start - angle_end
    else:
        angle = angle_end - angle_start
    return offset, angle


def angle_principal_measure(angle, min_angle=-math.pi):
    """
    Returns angle between O and 2 pi.
    """
    max_angle = min_angle + vm.TWO_PI
    angle = angle % vm.TWO_PI

    if math.isclose(angle, min_angle, abs_tol=1e-9):
        return min_angle
    if math.isclose(angle, max_angle, abs_tol=1e-9):
        return max_angle

    return angle


def clockwise_angle(vector1, vector2):
    """
    Return the clockwise angle in radians between vector1 and vector2.
    """
    vector0 = vm.O2D
    if vector0 in (vector1, vector2):
        return 0

    dot_v1v2 = vector1.dot(vector2)
    norm_vec_1 = vector1.norm()
    norm_vec_2 = vector2.norm()
    sol = dot_v1v2 / (norm_vec_1 * norm_vec_2)
    cross_v1v2 = vector1.cross(vector2)
    sol = max(-1.0, min(1.0, sol))
    inner_angle = math.acos(sol)

    if cross_v1v2 < 0:
        return inner_angle

    return vm.TWO_PI - inner_angle
