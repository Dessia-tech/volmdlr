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

    ax = point1.plot(color='white')
    point2.plot(ax, color='c')


Distance between point and edge
*******************************

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

