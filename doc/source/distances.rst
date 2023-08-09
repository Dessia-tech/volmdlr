=================
D I S T A N C E S
=================

Distance between two points
***************************

.. plot::
    :include-source:
    :align: center

    import volmdlr


    point1 = volmdlr.Point2D(0.0, 0.0)
    point2 = volmdlr.Point2D(1.0, 1.0)

    distance_point1_point2 = point1.point_distance(point2)
    print(distance_point1_point2)

    ax = point1.plot(color='k')
    point2.plot(ax, color='c')


Distance between point and another object
*****************************************

To calculate the distance between a point and **any** other volmdlr object is enough to use just the `point_distance` method.
Here are some examples:

Distance LineSegment2D to Point2D:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import edges
    from volmdlr.core import EdgeStyle

    point1 = volmdlr.Point2D(0.0, 0.0)
    linesegment2d = edges.LineSegment2D(volmdlr.Point2D(1.0, 1.0), volmdlr.Point2D(2.0, -1))

    distance_linesegment2d_point1, other_point = linesegment2d.point_distance(point1, return_other_point=True)

    ax = point1.plot()
    linesegment2d.plot(ax, EdgeStyle(color='c'))
    other_point.plot(ax, 'g')

.. code-block:: python

   print('distance_linesegment2d_point1: ', distance_linesegment2d_point1)
   >>> distance_linesegment2d_point1:  1.3416407864998738



2-Dimensional
=============


Distance LineSegment2D-Point2D
------------------------------



Distance Arc2D-Point2D
----------------------

Distance ArcEllipse2D-Point2D
-----------------------------

Distance BSplineCurve2D-Point2D
-------------------------------


3-Dimensional
=============


Distance LineSegment3D-Point3D
------------------------------

Distance Arc3D-Point3D
----------------------

Distance ArcEllipse3D-Point3D
-----------------------------

BSplineCurve3D-Point3D
----------------------

Distance between edge to edge
*****************************

2-Dimensional
=============

Distance LineSegment2D-LineSegment2D
------------------------------------

Distance LineSegment2D-LineSegment2D
------------------------------------

Distance LineSegment2D-Arc2D
----------------------------

Distance LineSegment2D-ArcEllipse2D
-----------------------------------

Distance LineSegment2D-BSplineCurve2D
-------------------------------------

Distance Arc2D-Arc2D
--------------------

Distance Arc2D-ArcEllipse2D
---------------------------

Distance Arc3D-BSplineCurve2D
-----------------------------

Distance ArcEllipse2D-ArcEllipse2D
----------------------------------

Distance ArcEllipse2D-ArcEllipse2D
----------------------------------



Distance



Distance between point and face
*******************************

Distance between point and surface
**********************************

Distance between point and shell
********************************

