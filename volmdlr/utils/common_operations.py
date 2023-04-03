"""
Concatenate common operation for two or more objects.

"""
import matplotlib
import matplotlib.pyplot as plt

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
