"""
volmdlr utils for calculating curves intersections.

"""
import math

import volmdlr
from volmdlr.utils.common_operations import get_abscissa_discretization, get_plane_equation_coefficients


def circle_3d_line_intersections(circle_3d, line, abs_tol: float = 1e-6):
    """
    Calculates the intersections between a Circle3D and a Line3D.

    :param circle_3d: Circle3D or Arc3D
    :param line: Line3D to verify intersections
    :param abs_tol: Tolerance.
    :return: list of points intersecting Circle
    """
    intersections = []
    if not math.isclose(abs(circle_3d.frame.w.dot(volmdlr.Z3D)), 1, abs_tol=abs_tol):
        frame_mapped_circle = circle_3d.frame_mapping(circle_3d.frame, 'new')
        frame_mapped_line = line.frame_mapping(circle_3d.frame, 'new')
        circle_linseg_intersections = circle_3d_line_intersections(frame_mapped_circle, frame_mapped_line)
        for inter in circle_linseg_intersections:
            intersections.append(circle_3d.frame.local_to_global_coordinates(inter))
        return intersections
    distance_center_lineseg = line.point_distance(circle_3d.frame.origin)
    if distance_center_lineseg > circle_3d.radius:
        return []
    direction_vector = line.direction_vector()
    if math.isclose(line.point1.z, line.point2.z, abs_tol=1e-6) and \
            math.isclose(line.point2.z, circle_3d.frame.origin.z, abs_tol=abs_tol):
        if line.point1.is_close(circle_3d.center):
            point1 = line.point2
            vec = line.point1 - line.point2
        else:
            point1 = line.point1
            vec = line.point2 - line.point1
        quadratic_equation_a = vec.dot(vec)
        quadratic_equation_b = 2 * vec.dot(point1 - circle_3d.center)
        quadratic_equation_c = point1.dot(point1) + circle_3d.center.dot(circle_3d.center) - 2 * point1.dot(
            circle_3d.center) - circle_3d.radius ** 2

        delta = quadratic_equation_b ** 2 - 4 * quadratic_equation_a * quadratic_equation_c
        if delta < 0:  # No real solutions, no intersection
            return []
        if math.isclose(delta, 0, abs_tol=1e-7):  # One real solution, tangent intersection
            t_param = -quadratic_equation_b / (2 * quadratic_equation_a)
            return [point1 + t_param * vec]
        sqrt_delta = math.sqrt(delta)
        t_param = (-quadratic_equation_b + sqrt_delta) / (2 * quadratic_equation_a)
        s_param = (-quadratic_equation_b - sqrt_delta) / (2 * quadratic_equation_a)
        return [point1 + t_param * vec, point1 + s_param * vec]
    z_constant = circle_3d.frame.origin.z
    constant = (z_constant - line.point1.z) / direction_vector.z
    x_coordinate = constant * direction_vector.x + line.point1.x
    y_coordinate = constant * direction_vector.y + line.point1.y
    if math.isclose((x_coordinate - circle_3d.frame.origin.x) ** 2 + (y_coordinate - circle_3d.frame.origin.y) ** 2,
                    circle_3d.radius ** 2, abs_tol=1e-6):
        intersections = [volmdlr.Point3D(x_coordinate, y_coordinate, z_constant)]
    return intersections


def conic3d_line_intersections(conic3d, line3d, abs_tol: float = 1e-6):
    """
    Calculates the intersections between an Ellipse3D and a Line3D.

    :param conic3d: The Hyperbola 3D.
    :param line3d: The Line 3D.
    :param abs_tol: Tolerance.
    :return: list of points intersecting the Hyperbola 3D.
    """
    intersections = []
    if not math.isclose(abs(conic3d.frame.w.dot(volmdlr.Z3D)), 1, abs_tol=abs_tol):
        frame_mapped_conic3d = conic3d.frame_mapping(conic3d.frame, 'new')
        frame_mapped_line = line3d.frame_mapping(conic3d.frame, 'new')
        circle_linseg_intersections = conic3d_line_intersections(frame_mapped_conic3d, frame_mapped_line)
        for inter in circle_linseg_intersections:
            intersections.append(conic3d.frame.local_to_global_coordinates(inter))
        return intersections

    if line3d.point1.z == line3d.point2.z == conic3d.frame.origin.z:
        conic2d = conic3d.self_2d
        line2d = line3d.to_2d(conic3d.frame.origin, conic3d.frame.u, conic3d.frame.v)
        intersections_2d = conic2d.line_intersections(line2d)
        for intersection in intersections_2d:
            intersections.append(volmdlr.Point3D(intersection[0], intersection[1], conic3d.frame.origin.z))
        return intersections
    plane_lineseg_intersections = get_plane_line_intersections(conic3d.frame, line3d)
    if conic3d.point_belongs(plane_lineseg_intersections[0]):
        return plane_lineseg_intersections
    return []


def ellipse2d_line_intersections(ellipse2d, line2d):
    """
    Calculates the intersections between a line and an ellipse.

    :param ellipse2d: Ellipse to calculate intersections
    :param line2d: line to calculate intersections
    :return: list of points intersections, if there are any
    """
    if line2d.point_distance(ellipse2d.center) > ellipse2d.major_axis + 1e-6:
        return []
    theta = volmdlr.geometry.clockwise_angle(ellipse2d.major_dir, volmdlr.X2D)
    if not math.isclose(theta, 0.0, abs_tol=1e-6) and not math.isclose(theta, 2 * math.pi, abs_tol=1e-6):
        frame = volmdlr.Frame2D(ellipse2d.center, ellipse2d.major_dir, ellipse2d.minor_dir)
        frame_mapped_ellipse = ellipse2d.frame_mapping(frame, 'new')
        line_inters = frame_mapped_ellipse.line_intersections(line2d.frame_mapping(frame, 'new'))
        line_intersections = [frame.local_to_global_coordinates(point) for point in line_inters]
        return line_intersections

    if line2d.point2.x == line2d.point1.x:
        x1 = line2d.point1.x
        x2 = x1
        y1 = ellipse2d.minor_axis * math.sqrt((1 - x1 ** 2 / ellipse2d.major_axis ** 2))
        y2 = -y1
    else:
        m = (line2d.point2.y - line2d.point1.y) / (line2d.point2.x - line2d.point1.x)
        c = - m * (line2d.point1.x + ellipse2d.center.x) + line2d.point1.y + ellipse2d.center.y
        if ellipse2d.major_axis ** 2 * m ** 2 + ellipse2d.minor_axis ** 2 > c ** 2:
            x1 = - (2 * (ellipse2d.major_axis ** 2) * m * c + math.sqrt(
                (2 * (ellipse2d.major_axis ** 2) * m * c) ** 2 - 4 * (
                        ellipse2d.major_axis ** 2 * m ** 2 + ellipse2d.minor_axis ** 2) *
                ellipse2d.major_axis ** 2 * (c ** 2 - ellipse2d.minor_axis ** 2))) / (
                         2 * (ellipse2d.major_axis ** 2 * (m ** 2) +
                              ellipse2d.minor_axis ** 2))

            x2 = - (2 * (ellipse2d.major_axis ** 2) * m * c - math.sqrt(
                (2 * (ellipse2d.major_axis ** 2) * m * c) ** 2 - 4 * (
                        ellipse2d.major_axis ** 2 * m ** 2 + ellipse2d.minor_axis ** 2) *
                ellipse2d.major_axis ** 2 * (c ** 2 - ellipse2d.minor_axis ** 2))) / (
                         2 * (ellipse2d.major_axis ** 2 * (m ** 2) +
                              ellipse2d.minor_axis ** 2))
            y1 = m * x1 + c
            y2 = m * x2 + c
            point1 = volmdlr.Point2D(x1, y1)
            point2 = volmdlr.Point2D(x2, y2)
            if point1 == point2:
                return [point1]
            return [point1, point2]
    return []


def get_circle_intersections(circle1, circle2):
    """
    Calculates the intersections between two circle 2d.

    :param circle1: circle 1 verify intersection with bspline
    :param circle2: circle 2 to search for intersections.
    :return: a list with all intersections between the two circles.
    """
    if circle1.center.is_close(circle2.center):
        return []
    x0, y0 = circle1.center
    x1, y1 = circle2.center

    d_param = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

    # non-intersecting
    if d_param > circle1.radius + circle2.radius:
        return []
    # One circle within other
    if d_param < abs(circle1.radius - circle2.radius):
        return []
    # coincident circles
    if d_param == 0 and circle1.radius == circle2.radius:
        return []
    a = (circle1.radius ** 2 - circle2.radius ** 2 + d_param ** 2) / (2 * d_param)
    h_param = math.sqrt(circle1.radius ** 2 - a ** 2)
    x2 = x0 + a * (x1 - x0) / d_param
    y2 = y0 + a * (y1 - y0) / d_param
    x3 = x2 + h_param * (y1 - y0) / d_param
    y3 = y2 - h_param * (x1 - x0) / d_param

    x4 = x2 - h_param * (y1 - y0) / d_param
    y4 = y2 + h_param * (x1 - x0) / d_param

    return [volmdlr.Point2D(x3, y3), volmdlr.Point2D(x4, y4)]


def bspline_intersections_initial_conditions(primitive, bsplinecurve, resolution: float = 10):
    """
    Gets the initial conditions to calculate intersections between a bspline curve 2d and another edge 2d.

    :param primitive: primitive to verify intersection with bspline
    :param bsplinecurve: bsplinecurve to search for intersections.
    :param resolution: bspline discretization resolution, to search for initial intersection conditions.
    :return: a list with all initial sections where there may exist an intersection.
    """
    line_seg_class_ = getattr(volmdlr.edges, 'LineSegment'+bsplinecurve.__class__.__name__[-2:])
    abscissa1 = 0
    abscissa2 = bsplinecurve.length()
    bspline_discretized_points, points_abscissas = get_abscissa_discretization(bsplinecurve, abscissa1, abscissa2,
                                                                               max_number_points=resolution)
    if bsplinecurve.periodic:
        bspline_discretized_points += [bspline_discretized_points[0]]
        if points_abscissas[0] == 0.0:
            points_abscissas += [bsplinecurve.length()]
        else:
            points_abscissas += [points_abscissas[0]]
    param_intersections = []
    for point1, point2, abscissa1, abscissa2 in zip(bspline_discretized_points[:-1], bspline_discretized_points[1:],
                                                    points_abscissas[:-1], points_abscissas[1:]):
        line_seg = line_seg_class_(point1, point2)
        intersection = primitive.linesegment_intersections(line_seg)
        if intersection:
            param_intersections.append((abscissa1, abscissa2))
    if not param_intersections:
        param_intersections.append((0.0, bsplinecurve.length()))
    return param_intersections


def get_bsplinecurve_intersections(primitive, bsplinecurve, abs_tol: float = 1e-6):
    """
    Calculates the intersections between a primitive and a BSpline Curve.

    This method calculates the intersections between a given primitive and a BSpline Curve.
    The primitive can be any edge type or a plane3D. It first determines the initial intersection
    conditions using the `bspline_intersections_initial_conditions` function. Then, it iteratively
    checks for intersections by discretizing the BSpline Curve and checking for intersections with
    line segments. It maintains a list of intersections and validates each intersection against
    the specified tolerance. The method continues until there are no more intersection conditions
    to process. The resulting intersections are returned as a list.

    :param primitive: The primitive object to verify intersection with the BSpline Curve.
    It can be any edge type or a plane3D.
    :param bsplinecurve: The BSpline Curve object to search for intersections.
    :param abs_tol: The tolerance to be considered while validating an intersection.

    :return: A list with all intersections between the edge and BSpline Curve.
    :rtype: [volmdlr.Point3D].
    """
    param_intersections = bspline_intersections_initial_conditions(primitive, bsplinecurve)
    line_seg_class_ = getattr(volmdlr.edges, 'LineSegment' + bsplinecurve.__class__.__name__[-2:])
    intersections = []
    while True:
        if not param_intersections:
            break
        abscissa1, abscissa2 = param_intersections[0]
        discretized_points_between_1_2, points_abscissas = get_abscissa_discretization(bsplinecurve, abscissa1,
                                                                                       abscissa2, max_number_points=10)
        for point1, point2, abscissa_point1, abscissa_point2 in zip(
                discretized_points_between_1_2[:-1], discretized_points_between_1_2[1:],
                points_abscissas[:-1], points_abscissas[1:]):
            line_seg = line_seg_class_(point1, point2)
            intersection = primitive.linesegment_intersections(line_seg, 1e-6)
            if not intersection:
                continue
            if bsplinecurve.point_distance(intersection[0]) > 1e-6:
                param_intersections.insert(0, (abscissa_point1, abscissa_point2))
            elif not volmdlr.core.point_in_list(intersection[0], intersections):
                intersections.append(intersection[0])
        param_intersections.remove((abscissa1, abscissa2))
    return intersections


def conic_intersections(conic1, conic2, abs_tol: float = 1e-6):
    """
    Gets intersections between a two conic curves 3D.

    :param conic1: First conic curve 3D.
    :param conic2: Other conic curve 3D.
    :param abs_tol: tolerance.
    :return: A list of points, containing all intersections between the Line 3D and the Parabola3D.
    """
    intersections = []
    if conic1.frame.w.is_colinear_to(conic2.frame.w) and \
            math.isclose(conic1.frame.w.dot(conic2.frame.origin - conic1.frame.origin), 0, abs_tol=abs_tol):
        frame_mapped_conic1 = conic1.frame_mapping(conic1.frame, 'new')
        frame_mapped_conic2 = conic2.frame_mapping(conic1.frame, 'new')
        conic1_2d = frame_mapped_conic1.self_2d
        conic2_2d = frame_mapped_conic2.to_2d(frame_mapped_conic1.frame.origin,
                                              frame_mapped_conic1.frame.u,frame_mapped_conic1.frame.v)
        intersections_2d = conic1_2d.curve_intersections(conic2_2d, abs_tol)
        local_intersections = []
        for intersection in intersections_2d:
            local_intersections.append(volmdlr.Point3D(intersection[0], intersection[1], conic1.frame.origin.z))
        # circle_linseg_intersections = conic3d_line_intersections(frame_mapped_conic1, frame_mapped_conic2, abs_tol)
        for inter in local_intersections:
            intersections.append(conic1.frame.local_to_global_coordinates(inter))
        return intersections

    plane_intersections = get_two_planes_intersections(conic1.frame, conic2.frame)
    if not plane_intersections:
        return []
    plane_intersections = volmdlr.curves.Line3D(plane_intersections[0], plane_intersections[1])
    self_ellipse3d_line_intersections = conic3d_line_intersections(conic1, plane_intersections)
    ellipse3d_line_intersections = conic3d_line_intersections(conic2, plane_intersections)
    for intersection in self_ellipse3d_line_intersections + ellipse3d_line_intersections:
        if volmdlr.core.point_in_list(intersection, intersections):
            continue
        if conic1.point_belongs(intersection, abs_tol) and conic2.point_belongs(intersection, abs_tol):
            intersections.append(intersection)
    return intersections


def get_plane_linesegment_intersections(plane_frame, linesegment, abs_tol: float = 1e-6):
    """
    Gets the intersections of a plane a line segment 3d.

    :param plane_frame: the plane's frame.
    :param linesegment: other line segment.
    :param abs_tol: tolerance allowed.
    :return: a list with the intersecting point.
    """
    u_vector = linesegment.end - linesegment.start
    w_vector = linesegment.start - plane_frame.origin
    normaldotu = plane_frame.w.dot(u_vector)
    if normaldotu == 0.0 or math.isclose(plane_frame.w.unit_vector().dot(u_vector.unit_vector()),
                                         0.0, abs_tol=abs_tol):
        return []
    intersection_abscissea = - plane_frame.w.dot(w_vector) / normaldotu
    if intersection_abscissea < 0 or intersection_abscissea > 1:
        if math.isclose(abs(intersection_abscissea), 0, abs_tol=abs_tol):
            return [linesegment.start]
        if math.isclose(intersection_abscissea, 1, abs_tol=abs_tol):
            return [linesegment.end]
        return []
    return [linesegment.start + intersection_abscissea * u_vector]


def get_plane_line_intersections(plane_frame, line, abs_tol: float = 1e-8):
    """
    Find the intersection with a line.

    :param plane_frame: the plane's frame.
    :param line: Line to evaluate the intersection
    :type line: :class:`edges.Line`.
    :param abs_tol: tolerance.
    :type abs_tol: float.
    :return: ADD DESCRIPTION
    :rtype: List[volmdlr.Point3D]
    """
    u_vector = line.point2 - line.point1
    w_vector = line.point1 - plane_frame.origin
    if math.isclose(plane_frame.w.dot(u_vector), 0, abs_tol=abs_tol):
        return []
    intersection_abscissea = - plane_frame.w.dot(w_vector) / plane_frame.w.dot(u_vector)
    return [line.point1 + intersection_abscissea * u_vector]


def get_two_planes_intersections(plane1_frame, plane2_frame):
    """
    Calculates the intersections between two planes, given their frames.

    :param plane1_frame: Plane's 1 frame.
    :param plane2_frame: Plane's 2 frame.
    :return: A list containing two points that define an infinite line if there is any intersections,
    or an empty list if the planes are parallel.
    """
    if plane1_frame.w.is_colinear_to(plane2_frame.w):
        return []
    line_direction = plane1_frame.w.cross(plane2_frame.w)

    if line_direction.norm() < 1e-6:
        return None

    a1, b1, c1, d1 = get_plane_equation_coefficients(plane1_frame)
    a2, b2, c2, d2 = get_plane_equation_coefficients(plane2_frame)
    if not math.isclose(a1 * b2 - a2 * b1, 0.0, abs_tol=1e-10):
        x0 = (b1 * d2 - b2 * d1) / (a1 * b2 - a2 * b1)
        y0 = (a2 * d1 - a1 * d2) / (a1 * b2 - a2 * b1)
        point1 = volmdlr.Point3D(x0, y0, 0)
    elif a2 * c1 != a1 * c2:
        x0 = (c2 * d1 - c1 * d2) / (a2 * c1 - a1 * c2)
        z0 = (a1 * d2 - a2 * d1) / (a2 * c1 - a1 * c2)
        point1 = volmdlr.Point3D(x0, 0, z0)
    elif c1 * b2 != b1 * c2:
        y0 = (- c2 * d1 + c1 * d2) / (b1 * c2 - c1 * b2)
        z0 = (- b1 * d2 + b2 * d1) / (b1 * c2 - c1 * b2)
        point1 = volmdlr.Point3D(0, y0, z0)
    else:
        raise NotImplementedError
    return [point1, point1 + line_direction]
