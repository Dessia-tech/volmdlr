==============================================
Section on how to create Basic volmdlr objects
==============================================

Basic Objects
*************


How to create a Vector
======================

To instantiate the Vector, you need to provide the required arguments. In this case, the constructor
expects two mandatory arguments x and y for the 2-dimensional vector, and an additional argument z for the
3-dimensional vector, which are floats representing the vector's coordinates.
Additionally, there is an optional argument name, which is a string representing the vector's name.

.. grid:: 2

    .. grid-item-card::  Vector2D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           import matplotlib.pyplot as plt
           import mplcyberpunk
           plt.style.use("cyberpunk")

           vector = volmdlr.Vector2D(1.0, 1.0)
           vector.plot(color='orange')

    .. grid-item-card::  Vector3D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           import matplotlib.pyplot as plt
           import mplcyberpunk
           plt.style.use("cyberpunk")

           vector = volmdlr.Vector3D(1.0, 1.0, 1.0)
           vector.plot(color='orange')

How to create a Point
======================

Point 2D
--------
To instantiate the Point2D class, you need to provide the two mandatory arguments required x and y,
which are both floats representing the vector's coordinates. Additionally, there is an optional argument name,
which is a string representing the vector's name.

.. code-block:: python

   import volmdlr

   point2d = volmdlr.Point2D(x=1.0, y=2.0, name='vector name')
   print('point2d:', point2d)
   >>> point2d: Point2D: [1.0, 2.0]

Point 3D
--------
To instantiate the Point3D class, you need to provide three mandatory arguments required x, y, z
which are floats representing the vector's coordinates. Additionally, there is an optional argument name,
which is a string representing the vector's name.

.. code-block:: python

   import volmdlr

   point2d = volmdlr.Point3D(x=1.0, y=2.0, name='vector name')
   print('point2d:', point2d)
   >>> point2d: Point2D: [1.0, 2.0]

How to create a Frame3D
=======================
The given class Frame3D represents a 3D frame. It defines a frame by specifying its origin point and
three basis vectors (u, v, and w) that determine the orientation of the frame.

Exmaple:

.. plot::
   :include-source:
   :align: center

   import volmdlr
   import matplotlib.pyplot as plt
   import mplcyberpunk
   plt.style.use("cyberpunk")

   origin = volmdlr.Point3D(0, 0, 0)
   u = volmdlr.Vector3D(1, 0, 0)
   v = volmdlr.Vector3D(0, 1, 0)
   w = volmdlr.Vector3D(0, 0, 1)
   frame = volmdlr.Frame3D(origin, u, v, w)
   frame.plot()

Curves
******

How to create Line in 2D and 3D
===============================

Line2D and Line3D represents an infinite lines in both 2 and 3D that passes through two points.
They are a subclass of Line, which handles line-related operations. The class takes two Point objects
as inputs to define the line and an optional name for identification.

To instantiate then, you need to create an object of of the corresponding class by calling its constructor (__init__)
and providing the required arguments. Here's how you can do it:

.. grid:: 2

    .. grid-item-card::  Line2D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import curves
           from volmdlr.core import EdgeStyle
           import matplotlib.pyplot as plt
           import mplcyberpunk
           plt.style.use("cyberpunk")

           point1 = volmdlr.Point2D(1.0, 1.0)
           point2 = volmdlr.Point2D(-2.0, -3.0)
           line2d = curves.Line2D(point1, point2, name='line2d_name_is_optional')
           line2d.plot(edge_style=EdgeStyle('r'))

    .. grid-item-card::  Line3D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import curves
           from volmdlr.core import EdgeStyle
           import matplotlib.pyplot as plt
           import mplcyberpunk
           plt.style.use("cyberpunk")

           point1 = volmdlr.Point3D(1.0, 1.0, 1.0)
           point2 = volmdlr.Point3D(-2.0, -3.0, -1.0)
           line3d = curves.Line3D(point1, point2, name='line3d_name_is_optional')
           line3d.plot(color='r')



How to create a Circle
=====================

Todo: describe how to create circle

.. grid:: 2

    .. grid-item-card::  Circle2D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import curves
           from volmdlr.core import EdgeStyle
           import matplotlib.pyplot as plt
           import mplcyberpunk
           plt.style.use("cyberpunk")

           center2d = volmdlr.Point2D(0.0, 0.0)
           circle2d = curves.Circle2D(center=center2d, radius=1, name='optional_circle_name')
           circle2d.plot(edge_style=EdgeStyle('orange'))


    .. grid-item-card::  Circle3D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import curves
           from volmdlr.core import EdgeStyle
           import matplotlib.pyplot as plt
           import mplcyberpunk
           plt.style.use("cyberpunk")

           center3D = volmdlr.Point3D(0.0, 0.0, 0.0)
           u_vector = volmdlr.Vector3D(1.0, 0.0, 0.0)
           v_vector = volmdlr.Vector3D(0.0, 1.0, 0.0)
           w_vector = volmdlr.Vector3D(0.0, 0.0, 1.0)
           frame3d = volmdlr.Frame3D(center3D, u_vector, v_vector, w_vector)
           circle3d = curves.Circle3D(frame=frame3d, radius=1, name='optional_circle_name')
           circle3d.plot(edge_style=EdgeStyle('orange'))

How to create an Ellipse
========================


Edges
*****

How to create a LineSegment
===========================

LineSegment2D
-------------

LineSegment3D
-------------

How to create an Arc
====================
Arc2D
-----
Arc3D
-----

How to create an ArcEllipse
===========================
ArcEllipse2D
------------
ArcEllipse3D
------------

How to create a BSplineCurve
============================
BSplineCurve2D
--------------
BSplineCurve3D
--------------

Wires
*****

Surfaces
********

Faces
*****

Shells
*****
