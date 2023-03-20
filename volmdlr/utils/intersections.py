"""
volmdlr utils for calculating curves intersections.

"""
import math

import volmdlr


def circle_3d_line_intersections(circle_3d, line):
    """
    Calculates the intersections between a Circle3D and a Line3D.

    :param circle_3d: Circle3D or Arc3D
    :param line: Line3D to verify intersections
    :return: list of points intersecting Circle
    """
    if not math.isclose(abs(circle_3d.frame.w.dot(volmdlr.Z3D)), 1, abs_tol=1e-6):
        frame_mapped_circle = circle_3d.frame_mapping(circle_3d.frame, 'new')
        frame_mapped_lineseg = line.frame_mapping(circle_3d.frame, 'new')
        circle_linseg_intersections = frame_mapped_circle.line_intersections(frame_mapped_lineseg)
        intersections = []
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
        return [volmdlr.Point3D(x_coordinate, y_coordinate, z_constant)]
    return []


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
