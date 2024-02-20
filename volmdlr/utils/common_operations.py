"""
Concatenate common operation for two or more objects.

"""
import math
import random

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
import scipy.integrate as scipy_integrate
from sklearn.cluster import DBSCAN

import volmdlr.core
from volmdlr.core import EdgeStyle


def plot_circle(circle, ax=None, edge_style: EdgeStyle = EdgeStyle()):
    """
    Create a Matplotlib plot for a circle 2d or fullarc 2d.

    :param circle: circle to plot.
    :param ax: Matplotlib plot axis.
    :param edge_style: Edge Style to implement.
    :return: Matplotlib plot axis.
    """
    if ax is None:
        _, ax = plt.subplots()
        x_min, x_max = circle.center[0] - circle.radius, circle.center[0] + circle.radius
        y_min, y_max = circle.center[1] - circle.radius, circle.center[1] + circle.radius
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
    if circle.radius > 0:
        ax.add_patch(matplotlib.patches.Arc((circle.center.x, circle.center.y),
                                            2 * circle.radius,
                                            2 * circle.radius,
                                            angle=0,
                                            theta1=0,
                                            theta2=360,
                                            color=edge_style.color,
                                            alpha=edge_style.alpha,
                                            linestyle=edge_style.linestyle,
                                            linewidth=edge_style.linewidth))
    if edge_style.plot_points:
        ax.plot([circle.start.x], [circle.start.y], 'o',
                color=edge_style.color, alpha=edge_style.alpha)
    if edge_style.equal_aspect:
        ax.set_aspect('equal')

    return ax


def random_color():
    """Random color generator."""
    return random.random(), random.random(), random.random()


def split_wire_by_plane(wire, plane3d):
    """
    Splits a wire into two parts using a plane.

    This method splits a wire into two parts based on the intersection points between the wire's primitives
    (edges) and a given 3D plane. It first finds the intersection points between each primitive and the plane,
    excluding duplicate points. Then, it checks if the number of intersection points is greater than one. If so,
    it raises a NotImplementedError, as the split is ambiguous. Otherwise, it performs the split using the
    `split_with_sorted_points` method of the wire object. The resulting wire objects are returned as a tuple.
    Note: The method assumes that the wire and the plane are in the same coordinate system.

    :param wire: The wire object to be split.
    :param plane3d: The 3D plane object used for splitting the wire.
    :return: A tuple containing two wire objects resulting from the split.
    :raises: NotImplementedError: If the wire intersects the plane at more than one point.

    :Example:
    >>> from volmdlr import Point3D, Vector3D
    >>> from volmdlr.surfaces import Plane3D
    >>> from volmdlr.core import EdgeStyle
    >>> from volmdlr.utils.common_operations import random_color
    >>> from volmdlr.models.open_rounded_line_segments import open_rounded_line_segements
    >>> plane = Plane3D.from_plane_vectors(Point3D(0.4, 0.4, 0.2), Vector3D(1, 0, 0), Vector3D(0, 1, 0))
    >>> split_wire1,split_wire2 = split_wire_by_plane(open_rounded_line_segements, plane)
    >>> ax = open_rounded_line_segements.plot()
    >>> plane.plot(ax)
    >>> split_wire1.plot(ax, EdgeStyle(random_color()))
    >>> split_wire2.plot(ax, EdgeStyle(random_color()))
    """
    wire_plane_intersections = []
    for primitive in wire.primitives:
        intersections = plane3d.edge_intersections(primitive)
        for intersection in intersections:
            if not intersection.in_list(wire_plane_intersections):
                wire_plane_intersections.append(intersection)
    if len(wire_plane_intersections) > 1:
        raise NotImplementedError
    wire1, wire2 = wire.split_with_sorted_points([wire_plane_intersections[0], wire.primitives[-1].end])
    return wire1, wire2


def plot_components_from_points(points, close_plot: bool = False):
    """
    Gets Matplotlib components from points.

    :param points: given points.
    :param close_plot: Weather to close the plot or not.
    :return:
    """
    components = [[], [], []]
    for point in points:
        for i, component in enumerate(point):
            components[i].append(component)
    valid_components = []
    for list_components in components:
        if list_components:
            if close_plot:
                list_components.append(list_components[0])
            valid_components.append(list_components)
    return valid_components


def plot_from_discretization_points(ax, edge_style, element, number_points: int = None, close_plot: bool = False):
    """
    General plot method using discretization_points method to generate points.

    :param ax: Matplotlib plot.
    :param edge_style: edge_style to be applied to plot.
    :param element: volmdlr element to be plotted (either 2D or 3D).
    :param number_points: number of points to be used in the plot.
    :param close_plot: specifies if plot is to be closed or not.
    :return: Matplotlib plot axis.
    """
    points = element.discretization_points(number_points=number_points)
    valid_components = plot_components_from_points(points, close_plot)
    ax.plot(*valid_components, color=edge_style.color, alpha=edge_style.alpha)
    return ax


def minimum_distance_points_circle3d_linesegment3d(circle3d,  linesegment3d):
    """
    Gets the points from the arc and the line that gives the minimal distance between them.

    :param circle3d: Circle 3d or Arc 3d.
    :param linesegment3d: Other line segment 3d.
    :return: Minimum distance points.
    """
    def distance_squared(x, u_param, v_param, k_param, w_param):
        """Calculates the squared distance."""
        return (u_param.dot(u_param) * x[0] ** 2 + w_param.dot(w_param) + v_param.dot(v_param) * (
                (math.sin(x[1])) ** 2) * radius ** 2 + k_param.dot(k_param) * ((math.cos(x[1])) ** 2) * radius ** 2
                - 2 * x[0] * w_param.dot(u_param) - 2 * x[0] * radius * math.sin(x[1]) * u_param.dot(v_param) - 2 * x[
                    0] * radius * math.cos(x[1]) * u_param.dot(k_param)
                + 2 * radius * math.sin(x[1]) * w_param.dot(v_param) +
                2 * radius * math.cos(x[1]) * w_param.dot(k_param)
                + math.sin(2 * x[1]) * v_param.dot(k_param) * radius ** 2)
    radius = circle3d.radius
    linseg_direction_vector = linesegment3d.direction_vector()
    vector_point_origin = circle3d.point_at_abscissa(0.0) - circle3d.frame.origin
    vector_point_origin = vector_point_origin.unit_vector()
    w = circle3d.frame.origin - linesegment3d.start
    v = circle3d.frame.w.cross(vector_point_origin)

    results = []
    for initial_value in [np.array([0.5, circle3d.angle / 2]), np.array([0.5, 0]), np.array([0.5, circle3d.angle])]:
        results.append(least_squares(distance_squared, initial_value,
                                     bounds=[(0, 0), (1, circle3d.angle)],
                                     args=(linseg_direction_vector, v, vector_point_origin, w)))

    point1 = linesegment3d.point_at_abscissa(results[0].x[0] * linesegment3d.length())
    point2 = circle3d.point_at_abscissa(results[1].x[1] * circle3d.radius)

    for couple in results[1:]:
        ptest1 = linesegment3d.point_at_abscissa(couple.x[0] * linesegment3d.length())
        ptest2 = circle3d.point_at_abscissa(couple.x[1] * circle3d.radius)
        if ptest1.point_distance(ptest2) < v.dot(v):
            point1, point2 = ptest1, ptest2

    return point1, point2


def get_abscissa_discretization(primitive, abscissa1, abscissa2, max_number_points: int = 10,
                                return_abscissas: bool = True):
    """
    Gets n discretization points between two given points of the edge.

    :param primitive: Primitive to discretize locally.
    :param abscissa1: Initial abscissa.
    :param abscissa2: Final abscissa.
    :param max_number_points: Expected number of points to discretize locally.
    :param return_abscissas: By default, returns also a list of abscissas corresponding to the
        discretization points
    :return: list of locally discretized point and a list containing the abscissas' values.
    """
    discretized_points_between_1_2 = []
    points_abscissas = []
    for abscissa in np.linspace(abscissa1, abscissa2, num=max_number_points):
        if abscissa > primitive.length() + 1e-6:
            continue
        abscissa_point = primitive.point_at_abscissa(abscissa)
        if not abscissa_point.in_list(discretized_points_between_1_2):
            discretized_points_between_1_2.append(abscissa_point)
            points_abscissas.append(abscissa)
    if return_abscissas:
        return discretized_points_between_1_2, points_abscissas
    return discretized_points_between_1_2


def get_point_distance_to_edge(edge, point, start, end):
    """
    Calculates the distance from a given point to an edge.

    :param edge: Edge to calculate distance to point.
    :param point: Point to calculate the distance to edge.
    :param start: Edge's start point.
    :param end: Edge's end point.
    :return: distance to edge.
    """
    best_distance = math.inf
    if start != end:
        if start.is_close(end):
            if not edge.periodic:
                return point.point_distance(start)
            number_points = 10 if abs(0 - edge.length()) > 5e-6 else 2
        else:
            number_points = 10 if abs(edge.abscissa(start) - edge.abscissa(end)) > 5e-6 else 2
        # number_points = 10 if abs(0 - edge.length()) > 5e-6 else 2
    elif edge.periodic:
        number_points = 10 if abs(0 - edge.length()) > 5e-6 else 2
    distance = best_distance
    point1_ = start
    point2_ = end
    linesegment_class_ = getattr(volmdlr.edges, 'LineSegment' + edge.__class__.__name__[-2:])
    while True:
        discretized_points_between_1_2 = edge.local_discretization(point1_, point2_, number_points)
        if not discretized_points_between_1_2:
            break
        distance = point.point_distance(discretized_points_between_1_2[0])
        for point1, point2 in zip(discretized_points_between_1_2[:-1], discretized_points_between_1_2[1:]):
            if point1.is_close(point2):
                continue
            line = linesegment_class_(point1, point2)
            dist = line.point_distance(point)
            if dist < distance:
                point1_ = point1
                point2_ = point2
                distance = dist
        if not point1_ or math.isclose(distance, best_distance, abs_tol=1e-7):
            break
        best_distance = distance
        # if math.isclose(abscissa1, abscissa2, abs_tol=1e-6):
        #     break
    return distance


def generic_minimum_distance(self, element, point1_edge1_, point2_edge1_, point1_edge2_,
                              point2_edge2_, return_points=False):
    """
    Gets the minimum distance between two elements.

    This is a generalized method in a case an analytical method has not yet been defined.

    :param element: other element.
    :param return_points: Weather to return the corresponding points or not.
    :return: distance to edge.
    """
    n = max(1, int(self.length() / element.length()))
    best_distance = math.inf
    distance_points = None
    distance = best_distance

    linesegment_class_ = getattr(volmdlr.edges, 'LineSegment' + self.__class__.__name__[-2:])
    while True:
        edge1_discretized_points_between_1_2 = self.local_discretization(point1_edge1_, point2_edge1_,
                                                                         number_points=10 * n)
        edge2_discretized_points_between_1_2 = element.local_discretization(point1_edge2_, point2_edge2_)
        if not edge1_discretized_points_between_1_2:
            break
        distance = edge2_discretized_points_between_1_2[0].point_distance(edge1_discretized_points_between_1_2[0])
        distance_points = [edge2_discretized_points_between_1_2[0], edge1_discretized_points_between_1_2[0]]
        for point1_edge1, point2_edge1 in zip(edge1_discretized_points_between_1_2[:-1],
                                              edge1_discretized_points_between_1_2[1:]):
            lineseg1 = linesegment_class_(point1_edge1, point2_edge1)
            for point1_edge2, point2_edge2 in zip(edge2_discretized_points_between_1_2[:-1],
                                                  edge2_discretized_points_between_1_2[1:]):
                lineseg2 = linesegment_class_(point1_edge2, point2_edge2)
                dist, min_dist_point1_, min_dist_point2_ = lineseg1.minimum_distance(lineseg2, True)
                if dist < distance:
                    point1_edge1_, point2_edge1_ = point1_edge1, point2_edge1
                    point1_edge2_, point2_edge2_ = point1_edge2, point2_edge2
                    distance = dist
                    distance_points = [min_dist_point1_, min_dist_point2_]
        if math.isclose(distance, best_distance, abs_tol=1e-6):
            break
        best_distance = distance
        n = 1
    if return_points:
        return distance, distance_points[0], distance_points[1]
    return distance


def ellipse_abscissa_angle_integration(ellipse3d, point_abscissa, angle_start, initial_angle):
    """
    Calculates the angle for a given abscissa point by integrating the ellipse.

    :param ellipse3d: the Ellipse3D.
    :param point_abscissa: the given abscissa for given point.
    :param angle_start: Ellipse3D / ArcEllipse3D start angle. (0 for Ellipse3D).
    :param initial_angle: angle abscissa's initial value.
    :return: final angle abscissa's value.
    """
    def ellipse_arc_length(theta):
        return math.sqrt((ellipse3d.major_axis ** 2) * math.sin(theta) ** 2 +
                         (ellipse3d.minor_axis ** 2) * math.cos(theta) ** 2)

    iter_counter = 0
    while True:
        res, _ = scipy_integrate.quad(ellipse_arc_length, angle_start, initial_angle)
        if math.isclose(res, point_abscissa, abs_tol=1e-8):
            abscissa_angle = initial_angle
            break
        if res > point_abscissa:
            increment_factor = (abs(initial_angle - angle_start) * (point_abscissa - res)) / (2 * abs(res))
        elif res == 0.0:
            increment_factor = 1e-5
        else:
            increment_factor = (abs(initial_angle - angle_start) * (point_abscissa - res)) / abs(res)
        initial_angle += increment_factor
        iter_counter += 1
    return abscissa_angle


def get_plane_equation_coefficients(plane_frame):
    """
    Returns the a,b,c,d coefficient from equation ax+by+cz+d = 0.

    """
    a, b, c = plane_frame.w
    d = -plane_frame.origin.dot(plane_frame.w)
    return round(a, 12), round(b, 12), round(c, 12), round(d, 12)


def get_plane_point_distance(plane_frame, point3d):
    """
    Gets distance between a point and plane, using its frame.

    :param plane_frame: plane's frame.
    :param point3d: other point.
    :return: point plane distance.
    """
    coefficient_a, coefficient_b, coefficient_c, coefficient_d = get_plane_equation_coefficients(plane_frame)
    return abs(plane_frame.w.dot(point3d) + coefficient_d) / math.sqrt(coefficient_a ** 2 +
                                                                       coefficient_b ** 2 + coefficient_c ** 2)


def order_points_list_for_nearest_neighbor(points):
    """
    Given a list of unordered points defining a path, it will order these points considering the nearest neighbor.

    """
    ordered_points = []
    remaining_points = np.array([[*point] for point in points])
    current_point = remaining_points[0, :]
    remaining_points = np.delete(remaining_points, 0, axis=0)

    while remaining_points.any():
        nearest_point_idx = np.argmin(np.linalg.norm(remaining_points - current_point, axis=1))
        nearest_point = remaining_points[nearest_point_idx, :]
        remaining_points = np.delete(remaining_points, nearest_point_idx, axis=0)
        ordered_points.append(volmdlr.Point3D(*current_point))
        current_point = nearest_point

    # Add the last point to complete the loop
    ordered_points.append(volmdlr.Point3D(*current_point))

    return ordered_points


def separate_points_by_closeness(points):
    """
    Separates a list of 3D Cartesian points into two groups based on their spatial closeness using DBSCAN.

    This function applies the DBSCAN (Density-Based Spatial Clustering of Applications with Noise) algorithm to
    the given list of 3D Cartesian points. DBSCAN clusters the points based on their spatial proximity.
    The points are separated into two groups, 'group1' and 'group2', depending on their spatial closeness
    as determined by the DBSCAN clustering.

    Please note that the 'eps' parameter inside the function can be adjusted to control the closeness threshold.

    :param points: A list of 3D Cartesian points, where each point is represented as a list of three coordinates.

    :return:
    - group1 (list of lists): The first group of points based on their closeness.
    - group2 (list of lists): The second group of points based on their closeness.
    """
    if not points:
        return []
    points_ = np.array([[*point] for point in points])

    # Apply DBSCAN clustering with a small epsilon to separate close points
    distances = sorted(np.linalg.norm(points_[1:] - points_[0], axis=1))
    eps = max(min(np.mean(distances[:max(int(len(points)*0.1), 30)]) / 2, 0.35), 0.02)
    # eps = np.mean(distances[:max(int(len(points)*0.1), 30)]) / 2

    dbscan = DBSCAN(eps=eps, min_samples=1)
    labels = dbscan.fit_predict(points_)

    # Initialize two empty lists for the two groups
    groups = {}
    # Assign points to group1 or group2 based on DBSCAN labels
    for i, label in enumerate(labels):
        if label not in groups:
            groups[label] = [points[i]]
            continue
        groups[label].append(points[i])
    keys = list(groups.keys())
    for key in keys:
        groups[key] = order_points_list_for_nearest_neighbor(groups[key])
        groups[key].append(groups[key][0])
    return list(groups.values())


def get_center_of_mass(list_points):
    """
    Gets the center of mass of a given list of points.

    :param list_points: list of points to get the center of mass from.
    :return: center of mass point.
    """
    center_mass = volmdlr.O3D
    for point in list_points:
        center_mass += point
    center_mass /= len(list_points)
    return center_mass
