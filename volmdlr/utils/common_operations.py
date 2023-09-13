"""
Concatenate common operation for two or more objects.

"""
import math
import random

import numpy as np
from scipy.optimize import least_squares


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
            if not point_in_list(intersection, wire_plane_intersections):
                wire_plane_intersections.append(intersection)
    if len(wire_plane_intersections) > 1:
        raise NotImplementedError
    wire1, wire2 = wire.split_with_sorted_points([wire_plane_intersections[0], wire.primitives[-1].end])
    return wire1, wire2


def plot_from_discretization_points(ax, edge_style, element, number_points: int = None, close_plot: bool = False):
    """
    General plot method using discretization_points method to generate points.

    :param ax: Matplotlib plot.
    :param edge_style: edge_style to be applied to plot.
    :param element: volmdlr element to be plotted (either 2D or 3D).
    :param number_points: number of points to be used in the plot.
    :param close_plot: specifies if plot is to be closed or not.
    :return: Matplolib plot axis.
    """
    components = [[], [], []]
    for point in element.discretization_points(number_points=number_points):
        for i, component in enumerate(point):
            components[i].append(component)
    valid_components = []
    for list_components in components:
        if list_components:
            if close_plot:
                list_components.append(list_components[0])
            valid_components.append(list_components)
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
    circle_point = circle3d.point_at_abscissa(0.0)
    radius = circle3d.radius
    linseg_direction_vector = linesegment3d.direction_vector()
    vector_point_origin = circle_point - circle3d.frame.origin
    vector_point_origin.normalize()
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
        dtest = ptest1.point_distance(ptest2)
        if dtest < v.dot(v):
            point1, point2 = ptest1, ptest2

    return point1, point2


def element_in_list(element, list_elements, tol: float = 1e-6):
    """
    Verifies if a volmdlr element is inside a list  of elements, considering a certain tolerance.

    :param element: Element to be verified inside list.
    :param list_elements: List of elements to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: True if there is an element inside the list close to the element to given tolerance.
    """
    for element_i in list_elements:
        if element.is_close(element_i, tol):
            return True
    return False


def point_in_list(point, list_points, tol: float = 1e-6):
    """
    Verifies if a point is inside a list  of points, considering a certain tolerance.

    :param point: Point to be verified inside list.
    :param list_points: List of points to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: True if there is a point inside the list close to the point to given tolerance.
    """

    return element_in_list(point, list_points, tol)


def edge_in_list(edge, list_edges, tol: float = 1e-6):
    """
    Verifies if an edge is inside a list  of edges, considering a certain tolerance.

    :param edge: Edge to be verified inside list.
    :param list_edges: List of edges to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: True if there is an edge inside the list close to the edge to given tolerance.
    """

    return element_in_list(edge, list_edges, tol)


def get_element_index_in_list(element, list_elements, tol: float = 1e-6):
    """
    Gets the index an element inside a list of elements, considering a certain tolerance.

    :param point: Element to be verified inside list.
    :param list_elements: List of elements to be used.
    :param tol: Tolerance to consider if two elements are the same.
    :return: The element index.
    """
    for i, element_i in enumerate(list_elements):
        if element_i.is_close(element, tol):
            return i
    raise ValueError(f'{element} is not in list')


def get_point_index_in_list(point, list_points, tol: float = 1e-6):
    """
    Gets the index a point inside a list of points, considering a certain tolerance.

    :param point: Point to be verified inside list.
    :param list_points: List of points to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: The point index.
    """

    return get_element_index_in_list(point, list_points, tol)


def get_edge_index_in_list(edge, list_edges, tol: float = 1e-6):
    """
    Gets the index a edge inside a list of edges, considering a certain tolerance.

    :param edge: Edge to be verified inside list.
    :param list_edges: List of edges to be used.
    :param tol: Tolerance to consider if two edges are the same.
    :return: The edge index.
    """

    return get_element_index_in_list(edge, list_edges, tol)


def determinant(vec1, vec2, vec3):
    """
    Calculates the determinant for a three vector matrix.

    """
    # TODO: to be removed
    a = np.array((vec1.vector, vec2.vector, vec3.vector))
    return np.linalg.det(a)


def delete_double_point(list_point):
    """
    Delete duplicate points from a list of points.

    :param list_point: The initial list of points
    :type list_point: Union[List[:class:`volmdlr.Point2D`],
        List[:class:`volmdlr.Point3D`]]
    :return: The final list of points containing no duplicates
    :rtype: Union[List[:class:`volmdlr.Point2D`],
        List[:class:`volmdlr.Point3D`]]
    """
    # TODO : this method would be faster using sets
    points = []
    for point in list_point:
        if point not in points:
            points.append(point)
        else:
            continue
    return points

# def get_triangulation(primitives):
#     display_meshes = []
#     for primitive in primitives:
#         display_meshes.append(primitive.triangulation())
#     return DisplayMesh3D.merge_meshes(display_meshes)
