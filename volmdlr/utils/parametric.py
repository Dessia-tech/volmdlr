"""
volmdlr utils for calculating 3D to surface parametric domain operationa
"""
import math
import bisect

import volmdlr
import volmdlr.edges as vme


def repair_singularity(primitive, last_primitive):
    """
    Repairs the Contour2D of SphericalSurface3D and ConicalSurface3D parametric face representations.

    Used when transforming from spatial to parametric coordinates when the surface contains a sigularity
    """
    v1 = primitive.unit_direction_vector()
    v2 = last_primitive.unit_direction_vector()
    dot = v1.dot(volmdlr.X2D)
    cross = v1.cross(v2)
    new_primitives = []
    if cross == 0 and dot == 0:
        if primitive.start.x == math.pi:
            primitive = primitive.translation(volmdlr.Vector2D(-2 * math.pi, 0))
            new = vme.LineSegment2D(last_primitive.end, primitive.start)
        elif primitive.start.x == -math.pi:
            primitive = primitive.translation(volmdlr.Vector2D(2 * math.pi, 0))
            new = vme.LineSegment2D(last_primitive.end, primitive.start)
        else:
            new = vme.LineSegment2D(last_primitive.end, primitive.start)

        new_primitives.append(new)
        new_primitives.append(primitive)
    else:
        delta = last_primitive.end - primitive.start
        new_primitives.append(primitive.translation(delta))
    return new_primitives


def repair_start_end_angle_periodicity(angle, ref_angle):
    """
    Repairs start and end angles in parametric coordinates.

    Uses ref_angle (angle just after start angle, if repairing start angle or angle just before end if repairing end
        angle).

    :param angle: Angle to repair.
    :type angle: float
    :param ref_angle: Angle of reference.
    :return: Returns the repaired angle.
    :rtype: float
    """
    # Verify if theta1 or theta2 point should be -pi because atan2() -> ]-pi, pi]
    if math.isclose(angle, math.pi, abs_tol=1e-6) and ref_angle < 0:
        angle = -math.pi
    elif math.isclose(angle, -math.pi, abs_tol=1e-6) and ref_angle > 0:
        angle = math.pi
    return angle


def repair_arc3d_angle_continuity(angle_start, angle_after_start, angle_end, angle3d):
    """
    Repairs Arc3D continuity after convertion of points to parametric 2D space.
    """
    ref_low = angle_start - angle3d
    ref_up = angle_start + angle3d

    # angle_after_start < angle_start --> angle coordinate axis going clockwise
    # ref_low < -math.pi -> crossing lower bound of atan2  [-math.pi, math.pi]
    if angle_after_start < angle_start and ref_low < -math.pi:
        angle_end = ref_low

    # angle_after_start > angle_start --> angle coordinate axis going trigowise
    #  ref_up > math.pi -> crossing lower bound of atan2  [-math.pi, math.pi]
    elif angle_after_start > angle_start and ref_up > math.pi:
        angle_end = ref_up

    if angle_start > 0 > angle_after_start:
        angle_start -= 2 * math.pi
    elif angle_start < 0 < angle_after_start:
        angle_start += 2 * math.pi

    return angle_start, angle_end


def arc3d_to_cylindrical_verification(start, end, angle3d, theta3, theta4):
    """
    Verifies theta from start and end of an Arc3D after transformation from spatial to parametric coordinates.
    """
    theta1, z1 = start
    theta2, z2 = end

    if abs(theta1) == math.pi:
        theta1 = repair_start_end_angle_periodicity(theta1, theta3)
    if abs(theta2) == math.pi:
        theta2 = repair_start_end_angle_periodicity(theta2, theta4)

    theta1, theta2 = repair_arc3d_angle_continuity(theta1, theta3, theta2, angle3d)

    start = volmdlr.Point2D(theta1, z1)
    end = volmdlr.Point2D(theta2, z2)
    return [start, end]


def arc3d_to_spherical_verification(start, end, angle3d, point_after_start, point_before_end):
    """
    Verifies theta and phi from start and end of an arc3d after transformation from spatial to parametric coordinates.
    """
    theta1, phi1 = start
    theta2, phi2 = end
    theta3, phi3 = point_after_start
    theta4, phi4 = point_before_end
    # Verify if theta1 or theta2 point should be -pi or pi because atan2() -> ]-pi, pi]
    if abs(theta1) == math.pi:
        theta1 = repair_start_end_angle_periodicity(theta1, theta3)
    if abs(theta2) == math.pi:
        theta2 = repair_start_end_angle_periodicity(theta2, theta4)

    # Verify if phi1 or phi2 point should be -pi or pi because phi -> ]-pi, pi]
    if abs(phi1) == math.pi:
        phi1 = repair_start_end_angle_periodicity(phi1, phi3)
    if abs(phi2) == math.pi:
        phi2 = repair_start_end_angle_periodicity(phi2, phi4)

    if math.isclose(phi1, phi2, abs_tol=1e-4):
        theta1, theta2 = repair_arc3d_angle_continuity(theta1, theta3, theta2, angle3d)

    if math.isclose(theta1, theta2, abs_tol=1e-4):
        phi1, phi2 = repair_arc3d_angle_continuity(phi1, phi3, phi2, angle3d)

    start = volmdlr.Point2D(theta1, phi1)
    end = volmdlr.Point2D(theta2, phi2)

    return start, end


def array_range_search(x, xmin, xmax):
    """
    Find the indices of the elements in the sorted list `x` that fall within the specified range.

    This function use bisect pyhton builtin module, which uses binary search and has a time complexity of O(log(n)).
    Where n is the array length.

    :param x: A sorted list of values.
    :type x: list
    :param xmin: The minimum value in the range.
    :type xmin: float
    :param xmax: The maximum value in the range.
    :type xmax: float
    :return: A python range from the first to the last elements in `x`.
    :rtype: range

    :Example:

    >>> x = [1, 2, 3, 4, 5]
    >>> array_range_search(x, 2, 4)
    range(1, 3)
    """

    left = bisect.bisect_left(x, xmin)
    right = bisect.bisect_right(x, xmax)
    return range(left, right)
