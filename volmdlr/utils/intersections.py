"""
volmdlr utils for calculating curves intersetions
"""
import math

import volmdlr


def circle_3d_linesegment_intersections(circle_3d, linesegment):
    """
    Calculates the intersections between a Circle3D and a LineSegment3D.

    :param circle_3d: Circle3D or Arc3D
    :param linesegment: LineSegment3D to verify intersections
    :return: list of points intersecting Circle
    """
    if not math.isclose(abs(circle_3d.frame.w.dot(volmdlr.Z3D)), 1, abs_tol=1e-6):
        frame_mapped_circle = circle_3d.frame_mapping(circle_3d.frame, 'new')
        frame_mapped_lineseg = linesegment.frame_mapping(circle_3d.frame, 'new')
        circle_linseg_intersections = frame_mapped_circle.linesegment_intersections(frame_mapped_lineseg)
        intersections = []
        for inter in circle_linseg_intersections:
            intersections.append(circle_3d.frame.old_coordinates(inter))
        return intersections
    distance_center_lineseg = linesegment.point_distance(circle_3d.frame.origin)
    if distance_center_lineseg > circle_3d.radius:
        return []
    direction_vector = linesegment.direction_vector()
    if linesegment.start.z == linesegment.end.z == circle_3d.frame.origin.z:
        quadratic_equation_a = (1 + (direction_vector.y ** 2 / direction_vector.x ** 2))
        quadratic_equation_b = (-2 * (direction_vector.y ** 2 / direction_vector.x ** 2) * linesegment.start.x +
                                2 * (direction_vector.y / direction_vector.x) * linesegment.start.y)
        quadratic_equation_c = ((linesegment.start.y - (direction_vector.y / direction_vector.x) *
                                 linesegment.start.x) ** 2 - circle_3d.radius ** 2)
        delta = (quadratic_equation_b ** 2 - 4 * quadratic_equation_a * quadratic_equation_c)
        x1 = (- quadratic_equation_b + math.sqrt(delta)) / (2 * quadratic_equation_a)
        x2 = (- quadratic_equation_b - math.sqrt(delta)) / (2 * quadratic_equation_a)
        y1 = (direction_vector.y / direction_vector.x) * (x1 - linesegment.start.x) + linesegment.start.y
        y2 = (direction_vector.y / direction_vector.x) * (x2 - linesegment.start.x) + linesegment.start.y
        return [volmdlr.Point3D(x1, y1, circle_3d.frame.origin.z), volmdlr.Point3D(x2, y2, circle_3d.frame.origin.z)]
    z_constant = circle_3d.frame.origin.z
    constant = (z_constant - linesegment.start.z) / direction_vector.z
    x_coordinate = constant * direction_vector.x + linesegment.start.x
    y_coordinate = constant * direction_vector.y + linesegment.start.y
    if math.isclose((x_coordinate - circle_3d.frame.origin.x) ** 2 + (y_coordinate - circle_3d.frame.origin.y) ** 2,
                    circle_3d.radius ** 2, abs_tol=1e-6):
        return [volmdlr.Point3D(x_coordinate, y_coordinate, z_constant)]
    return []
