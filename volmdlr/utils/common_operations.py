"""
Concatenate common operation for two or more objects.

"""
import math
import numpy as npy
import matplotlib
import matplotlib.pyplot as plt
import volmdlr.core
from volmdlr import edges

from volmdlr.core import EdgeStyle


def plot_circle(circle, ax=None, edge_style: EdgeStyle = EdgeStyle()):
    """
    Create a matplotlib plot for a circle 2d or fullarc 2d.

    :param circle: circle to plot.
    :param ax: matplotlib plot axis.
    :param edge_style: Edge Style to implement.
    :return: matplotlib plot axis.
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
    return ax


def point_distance_to_edge(edge, point):
    """
    Calculates the distance from a given point to an edge.

    :param point: point.
    :param edge: edge.
    :return: distance to edge.
    """
    best_distance = math.inf
    abscissa1 = 0
    abscissa2 = edge.abscissa(edge.end)
    distance = best_distance
    point1_ = None
    point2_ = None
    linesegment_class_ = getattr(edges, 'LineSegment'+edge.__class__.__name__[-2:])
    while True:
        discretized_points_between_1_2 = []
        for abscissa in npy.linspace(abscissa1, abscissa2, num=8):
            abscissa_point = edge.point_at_abscissa(abscissa)
            if not volmdlr.core.point_in_list(abscissa_point, discretized_points_between_1_2):
                discretized_points_between_1_2.append(abscissa_point)
        if not discretized_points_between_1_2:
            break
        distance = point.point_distance(discretized_points_between_1_2[0])
        for point1, point2 in zip(discretized_points_between_1_2[:-1], discretized_points_between_1_2[1:]):
            line = linesegment_class_(point1, point2)
            dist = line.point_distance(point)
            if dist < distance:
                point1_ = point1
                point2_ = point2
                distance = dist
        if not point1_ or math.isclose(distance, best_distance, abs_tol=1e-6):
            break
        abscissa1 = edge.abscissa(point1_)
        abscissa2 = edge.abscissa(point2_)
        best_distance = distance
        if math.isclose(abscissa1, abscissa2, abs_tol=1e-6):
            break
    return distance
