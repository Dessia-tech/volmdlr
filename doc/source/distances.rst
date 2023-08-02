=================
D I S T A N C E S
=================

Distance between two points
***************************

.. plot::
    :include-source:
    :align: center

    import volmdlr
    import matplotlib.pyplot as plt
    import mplcyberpunk
    plt.style.use("cyberpunk")

    point1 = volmdlr.Point2D(0.0, 0.0)
    point2 = volmdlr.Point2D(1.0, 1.0)

    distance_point1_point2 = point1.point_distance(point2)
    print(distance_point1_point2)

    ax = point1.plot(color='white')
    point2.plot(ax, color='c')


Distance between point and edge
*******************************

Distance LineSegment-Point
==========================

Distance Arc-Point
==================

Distance ArcEllipse-Point
=========================

Distance BSplineCurve-Point
==========================

Distance between point and face
*******************************

Distance between point and surface
**********************************

Distance between point and shell
********************************

