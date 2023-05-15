"""
volmdlr utils for calculating 3D to surface parametric domain operation.

"""
import bisect
import math
import warnings

import volmdlr
import volmdlr.edges as vme


def find_sign_changes(list_of_values):
    """
    Finds the position of sign changes in a list.

    :param list_of_values: The list to search for sign changes.
    :type list_of_values: list
    :returns: A list of indices where the sign changes occur.
    :rtype: list
    """
    sign_changes = []
    for i in range(1, len(list_of_values)):
        if list_of_values[i] * list_of_values[i - 1] < 0:
            sign_changes.append(i)
    return sign_changes


def angle_discontinuity(angle_list):
    """
    Returns True if there is some angle discontinuity in the angle_list.

    This discontinuity can occur when we transform some 3D primitive into the parametric domain of periodical surfaces.

    :param angle_list: List with angle axis values of the primitive in parametric domain.
    :type angle_list: List[float]
    :return: Returns True if there is discontinuity, False otherwise.
    :rtype: bool
    """
    indexes_sign_changes = find_sign_changes(angle_list)
    discontinuity = False
    indexes_angle_discontinuity = []
    if indexes_sign_changes:
        for index in indexes_sign_changes:
            sign = round(angle_list[index - 1] / abs(angle_list[index - 1]), 2)
            delta = abs(angle_list[index] + sign * volmdlr.TWO_PI - angle_list[index - 1])
            if math.isclose(abs(angle_list[index]), math.pi, abs_tol=delta) and \
                    not math.isclose(abs(angle_list[index]), 0, abs_tol=delta):
                indexes_angle_discontinuity.append(index)
                discontinuity = True
    return discontinuity, indexes_angle_discontinuity


def is_undefined_brep_primitive(primitive, periodicity):
    """
    2D primitives on parametric surface domain can be in wrong bound side due to periodic behavior of some surfaces.

    A B-Rep primitive can be undefined when it's not directly provided, but it's found by some 3D to 2D operation.
    This can result in a primitive that is entirely contained on the periodical boundary, and we can only know its
    right position when it's placed on the boundary representation, by analyzing the continuity of the 2D contour.

    :param primitive: primitive to perform verification.
    :type primitive: vme.Edge
    :param periodicity: list with periodicity in x and y direction
    :type periodicity: list
    """
    start = primitive.start
    end = primitive.end

    if periodicity[0] and not periodicity[1]:
        i = 0
    elif not periodicity[0] and periodicity[1]:
        i = 1
    else:
        return False
    #     if ((2 * start.x) % periodicity[0]) == 0 and ((2 * start.y) % periodicity[1]) == 0:
    #         return True
    #     return False
    if ((2 * start[i]) % periodicity[i]) == 0 and end[i] == start[i]:
        return True
    return False


def find_index_defined_brep_primitive_on_periodical_surface(primitives2d, periodicity):
    """
    Search for a primitive on the boundary representation that can be used as reference for reparing the contour.
    """
    pos = 0
    for i, primitive in enumerate(primitives2d):
        if not is_undefined_brep_primitive(primitive, periodicity):
            return i
    return pos


def repair_singularity(primitive, last_primitive):
    """
    Repairs the Contour2D of SphericalSurface3D and ConicalSurface3D parametric face representations.

    Used when transforming from spatial to parametric coordinates when the surface contains a singularity
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
    elif math.isclose(angle, 0.5 * math.pi, abs_tol=1e-6) and ref_angle < 0:
        angle = -0.5 * math.pi
    elif math.isclose(angle, -0.5 * math.pi, abs_tol=1e-6) and ref_angle > 0:
        angle = 0.5 * math.pi
    return angle


def repair_arc3d_angle_continuity(angle_start, angle_after_start, angle_end, angle3d, periodicity):
    """
    Repairs Arc3D continuity after conversion of points to parametric 2D space.
    """
    ref_low = angle_start - angle3d
    ref_up = angle_start + angle3d

    # angle_after_start < angle_start --> angle coordinate axis going clockwise
    # ref_low < -math.pi -> crossing lower bound of atan2  [-math.pi, math.pi]
    if angle_after_start < angle_start and ref_low < -math.pi:
        angle_end = ref_low

    # angle_after_start > angle_start --> angle coordinate axis going trigo wise
    #  ref_up > math.pi -> crossing lower bound of atan2  [-math.pi, math.pi]
    elif angle_after_start > angle_start and ref_up > math.pi:
        angle_end = ref_up
    if not math.isclose(angle_start, 0, abs_tol=1e-4):
        if angle_start > 0 > angle_after_start:
            angle_start -= periodicity
        elif angle_start < 0 < angle_after_start:
            angle_start += periodicity

    return angle_start, angle_end


def arc3d_to_cylindrical_coordinates_verification(start, end, angle3d, theta3, theta4):
    """
    Verifies theta from start and end of an Arc3D after transformation from spatial to parametric coordinates.
    """
    theta1, z1 = start
    theta2, z2 = end

    if abs(theta1) == math.pi or abs(theta1) == 0.5 * math.pi:
        theta1 = repair_start_end_angle_periodicity(theta1, theta3)
    if abs(theta2) == math.pi or abs(theta2) == 0.5 * math.pi:
        theta2 = repair_start_end_angle_periodicity(theta2, theta4)

    theta1, theta2 = repair_arc3d_angle_continuity(theta1, theta3, theta2, angle3d, volmdlr.TWO_PI)

    start = volmdlr.Point2D(theta1, z1)
    end = volmdlr.Point2D(theta2, z2)
    return [start, end]


def fullarc_to_cylindrical_coordinates_verification(start, end, theta3):
    """
    Verifies theta from start and end of a full arc after transformation from spatial to parametric coordinates.
    """
    theta1, _ = start
    _, z2 = end

    if theta1 > theta3:
        end = volmdlr.Point2D(theta1 + volmdlr.TWO_PI, z2)
    elif theta1 < theta3:
        end = volmdlr.Point2D(theta1 - volmdlr.TWO_PI, z2)
    return [start, end]


def toroidal_repair_start_end_angle_periodicity(start, end, point_after_start, point_before_end):
    """
    Verifies theta and phi from start and end of an arc 3D after transformation from spatial to parametric coordinates.
    """
    theta1, phi1 = start
    theta2, phi2 = end
    theta3, phi3 = point_after_start
    theta4, phi4 = point_before_end
    # Verify if theta1 or theta2 point should be -pi or pi because atan2() -> ]-pi, pi]
    if abs(theta1) == math.pi or abs(theta1) == 0.5 * math.pi:
        theta1 = repair_start_end_angle_periodicity(theta1, theta3)
    if abs(theta2) == math.pi or abs(theta2) == 0.5 * math.pi:
        theta2 = repair_start_end_angle_periodicity(theta2, theta4)

    # Verify if phi1 or phi2 point should be -pi or pi because phi -> ]-pi, pi]
    if abs(phi1) == math.pi or abs(phi1) == 0.5 * math.pi:
        phi1 = repair_start_end_angle_periodicity(phi1, phi3)
    if abs(phi2) == math.pi or abs(phi2) == 0.5 * math.pi:
        phi2 = repair_start_end_angle_periodicity(phi2, phi4)
    return volmdlr.Point2D(theta1, phi1), volmdlr.Point2D(theta2, phi2)


def spherical_repair_start_end_angle_periodicity(start, end, point_after_start, point_before_end):
    """
    Verifies theta and phi from start and end of an arc 3D after transformation from spatial to parametric coordinates.
    """
    theta1, phi1 = start
    theta2, phi2 = end
    theta3, _ = point_after_start
    theta4, _ = point_before_end
    # Verify if theta1 or theta2 point should be -pi or pi because atan2() -> ]-pi, pi]
    if abs(theta1) == math.pi or abs(theta1) == 0.5 * math.pi:
        theta1 = repair_start_end_angle_periodicity(theta1, theta3)
    if abs(theta2) == math.pi or abs(theta2) == 0.5 * math.pi:
        theta2 = repair_start_end_angle_periodicity(theta2, theta4)

    return volmdlr.Point2D(theta1, phi1), volmdlr.Point2D(theta2, phi2)


def arc3d_to_toroidal_coordinates_verification(start, end, angle3d, reference_points, periodicity):
    """
    Verifies theta and phi from start and end of an arc 3D after transformation from spatial to parametric coordinates.
    """
    point_after_start = reference_points[0]
    point_before_end = reference_points[1]
    theta3, phi3 = point_after_start
    start, end = toroidal_repair_start_end_angle_periodicity(start, end, point_after_start, point_before_end)
    theta1, phi1 = start
    theta2, phi2 = end
    if math.isclose(phi1, phi2, abs_tol=1e-4):
        theta1, theta2 = repair_arc3d_angle_continuity(theta1, theta3, theta2,
                                                       angle3d, periodicity[0])

    if math.isclose(theta1, theta2, abs_tol=1e-4):
        phi1, phi2 = repair_arc3d_angle_continuity(phi1, phi3, phi2,
                                                   angle3d, periodicity[1])

    return volmdlr.Point2D(theta1, phi1), volmdlr.Point2D(theta2, phi2)


def arc3d_to_spherical_coordinates_verification(start, end, angle3d, reference_points, periodicity):
    """
    Verifies theta and phi from start and end of an arc 3D after transformation from spatial to parametric coordinates.
    """
    point_after_start = reference_points[0]
    point_before_end = reference_points[1]
    theta3, _ = point_after_start
    start, end = spherical_repair_start_end_angle_periodicity(start, end, point_after_start, point_before_end)
    theta1, phi1 = start
    theta2, phi2 = end
    if math.isclose(phi1, phi2, abs_tol=1e-4):
        theta1, theta2 = repair_arc3d_angle_continuity(theta1, theta3, theta2,
                                                       angle3d, periodicity[0])

    return volmdlr.Point2D(theta1, phi1), volmdlr.Point2D(theta2, phi2)


def array_range_search(x, xmin, xmax):
    """
    Find the indices of the elements in the sorted list `x` that fall within the specified range.

    This function use bisect python builtin module, which uses binary search and has a time complexity of O(log(n)).
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


def contour2d_healing(contour2d):
    """
    Heals topologies incoherences on the boundary representation.
    """
    contour2d = contour2d_healing_self_intersection(contour2d)
    contour2d = contour2d_healing_close_gaps(contour2d)
    return contour2d


def contour2d_healing_close_gaps(contour2d):
    """
    Heals topologies incoherences on the boundary representation.
    """
    new_primitives = [contour2d.primitives[0]]
    for i, (prim1, prim2) in enumerate(
            zip(contour2d.primitives, contour2d.primitives[1:] + [contour2d.primitives[0]])):
        if prim1 and prim2:
            if not prim1.end.is_close(prim2.start):
                new_primitives.append(vme.LineSegment2D(prim1.end, prim2.start))
            if i < len(contour2d.primitives) - 1:
                new_primitives.append(prim2)
    contour2d.primitives = new_primitives

    return contour2d


def contour2d_healing_self_intersection(contour2d):
    """
    Heals topologies incoherences on the boundary representation.
    """
    for i, (prim1, prim2) in enumerate(
            zip(contour2d.primitives, contour2d.primitives[1:] + [contour2d.primitives[0]])):
        if not prim1.end.is_close(prim2.start):
            # check intersection
            intersections = prim1.intersections(prim2)
            if intersections:
                if len(intersections) > 1:
                    warnings.warn("More than one crossings found while detecting contour self intersection.")
                    return contour2d
                split_point = intersections[0]
                if prim1.is_point_edge_extremity(split_point):
                    new_prim1 = prim1
                else:
                    new_prim1 = prim1.split(split_point)[0]
                if prim2.is_point_edge_extremity(split_point):
                    new_prim2 = prim2
                else:
                    new_prim2 = prim2.split(split_point)[1]
                contour2d.primitives[i] = new_prim1
                contour2d.primitives[(i + 1) % len(contour2d.primitives)] = new_prim2
    return contour2d
