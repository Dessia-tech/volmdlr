"""
Concatenate common operation for two or more objects.

"""
import random

import matplotlib
import matplotlib.pyplot as plt

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
    x_min, x_max = circle.center[0] - circle.radius, circle.center[0] + circle.radius
    y_min, y_max = circle.center[1] - circle.radius, circle.center[1] + circle.radius
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

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
            if not volmdlr.core.point_in_list(intersection, wire_plane_intersections):
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
