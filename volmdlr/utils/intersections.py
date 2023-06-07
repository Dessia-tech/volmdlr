"""
volmdlr utils for calculating curves intersections.

"""
import math

import numpy as npy

import volmdlr


def circle_3d_line_intersections(circle_3d, line):
    """
    Calculates the intersections between a Circle3D and a Line3D.

    :param circle_3d: Circle3D or Arc3D
    :param line: Line3D to verify intersections
    :return: list of points intersecting Circle
    """
    intersections = []
    if not math.isclose(abs(circle_3d.frame.w.dot(volmdlr.Z3D)), 1, abs_tol=1e-6):
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
    if line.point1.z == line.point2.z == circle_3d.frame.origin.z:
        quadratic_equation_a = 1 + (direction_vector.y ** 2 / direction_vector.x ** 2)
        quadratic_equation_b = -2 * (direction_vector.y ** 2 / direction_vector.x ** 2) * line.point1.x + \
                               2 * (direction_vector.y / direction_vector.x) * line.point1.y
        quadratic_equation_c = (line.point1.y - (direction_vector.y / direction_vector.x) *
                                line.point1.x) ** 2 - circle_3d.radius ** 2
        delta = quadratic_equation_b ** 2 - 4 * quadratic_equation_a * quadratic_equation_c
        if delta < 0:  # No real solutions, no intersection
            return []
        if delta == 0:  # One real solution, tangent intersection
            x = -quadratic_equation_b / (2 * quadratic_equation_a)
            y = (direction_vector.y / direction_vector.x) * (x - line.point1.x) + line.point1.y
            return [volmdlr.Point3D(x, y, circle_3d.frame.origin.z)]
        x1 = (- quadratic_equation_b + math.sqrt(delta)) / (2 * quadratic_equation_a)
        x2 = (- quadratic_equation_b - math.sqrt(delta)) / (2 * quadratic_equation_a)
        y1 = (direction_vector.y / direction_vector.x) * (x1 - line.point1.x) + line.point1.y
        y2 = (direction_vector.y / direction_vector.x) * (x2 - line.point1.x) + line.point1.y
        return [volmdlr.Point3D(x1, y1, circle_3d.frame.origin.z), volmdlr.Point3D(x2, y2, circle_3d.frame.origin.z)]
    z_constant = circle_3d.frame.origin.z
    constant = (z_constant - line.point1.z) / direction_vector.z
    x_coordinate = constant * direction_vector.x + line.point1.x
    y_coordinate = constant * direction_vector.y + line.point1.y
    if math.isclose((x_coordinate - circle_3d.frame.origin.x) ** 2 + (y_coordinate - circle_3d.frame.origin.y) ** 2,
                    circle_3d.radius ** 2, abs_tol=1e-6):
        intersections = [volmdlr.Point3D(x_coordinate, y_coordinate, z_constant)]
    return intersections


def ellipse2d_line_intersections(ellipse2d, line2d):
    """
    Calculates the intersections between a line and an ellipse.

    :param ellipse2d: Ellipse to calculate intersections
    :param line2d: line to calculate intersections
    :return: list of points intersections, if there are any
    """
    theta = volmdlr.geometry.clockwise_angle(ellipse2d.major_dir, volmdlr.X2D)
    if not math.isclose(theta, 0.0, abs_tol=1e-6) and not math.isclose(theta, 2 * math.pi, abs_tol=1e-6):
        frame = volmdlr.Frame2D(ellipse2d.center, ellipse2d.major_dir, ellipse2d.minor_dir)
        frame_mapped_ellipse = ellipse2d.frame_mapping(frame, 'new')
        frame_mapped_line = line2d.frame_mapping(frame, 'new')
        line_inters = frame_mapped_ellipse.line_intersections(frame_mapped_line)
        line_intersections = [frame.local_to_global_coordinates(point) for point in line_inters]
        return line_intersections

    if line2d.point2.x == line2d.point1.x:
        x1 = line2d.point1.x
        x2 = x1
        y1 = ellipse2d.minor_axis * math.sqrt((1 - x1 ** 2 / ellipse2d.major_axis ** 2))
        y2 = -y1
        c = ellipse2d.center.y + line2d.point1.y
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

    d = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

    # non-intersecting
    if d > circle1.radius + circle2.radius:
        return []
    # One circle within other
    if d < abs(circle1.radius - circle2.radius):
        return []
    # coincident circles
    if d == 0 and circle1.radius == circle2.radius:
        return []
    a = (circle1.radius ** 2 - circle2.radius ** 2 + d ** 2) / (2 * d)
    h = math.sqrt(circle1.radius ** 2 - a ** 2)
    x2 = x0 + a * (x1 - x0) / d
    y2 = y0 + a * (y1 - y0) / d
    x3 = x2 + h * (y1 - y0) / d
    y3 = y2 - h * (x1 - x0) / d

    x4 = x2 - h * (y1 - y0) / d
    y4 = y2 + h * (x1 - x0) / d

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
    bspline_discretized_points = bsplinecurve.discretization_points(number_points=resolution)
    param_intersections = []
    for point1, point2 in zip(bspline_discretized_points[:-1], bspline_discretized_points[1:]):
        line_seg = line_seg_class_(point1, point2)
        intersection = primitive.linesegment_intersections(line_seg)
        if intersection:
            abscissa1 = bsplinecurve.abscissa(point1)
            abscissa2 = bsplinecurve.abscissa(point2)
            param_intersections.append((abscissa1, abscissa2))
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
        discretized_points_between_1_2 = []
        for abscissa in npy.linspace(abscissa1, abscissa2, num=10):
            abscissa_point = bsplinecurve.point_at_abscissa(abscissa)
            if not volmdlr.core.point_in_list(abscissa_point, discretized_points_between_1_2):
                discretized_points_between_1_2.append(abscissa_point)
        for point1, point2 in zip(discretized_points_between_1_2[:-1], discretized_points_between_1_2[1:]):
            line_seg = line_seg_class_(point1, point2)
            intersection = primitive.linesegment_intersections(line_seg, 1e-6)
            if not intersection:
                continue
            if bsplinecurve.point_distance(intersection[0]) > abs_tol:
                param_intersections.insert(0, (bsplinecurve.abscissa(point1),
                                               bsplinecurve.abscissa(point2)))
            elif not volmdlr.core.point_in_list(intersection[0], intersections):
                intersections.append(intersection[0])
        param_intersections.remove((abscissa1, abscissa2))
    return intersections
