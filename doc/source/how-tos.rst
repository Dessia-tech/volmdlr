==============================================
Section on how to create Basic volmdlr objects
==============================================

Basic Objects
*************


How to create a Vector
======================

Vector 2D
---------
To instantiate the Vector2D class, you need to provide the required arguments. In this case, the constructor
expects two mandatory arguments x and y, which are both floats representing the vector's coordinates.
Additionally, there is an optional argument name, which is a string representing the vector's name.

Here's how you can instantiate the Vector2D class:

.. plot::
   :include-source:
   :align: center

   import volmdlr
   import matplotlib.pyplot as plt
   import mplcyberpunk
   plt.style.use("cyberpunk")

   vector = volmdlr.Vector2D(1.0, 1.0)
   vector.plot(color='orange')

Vector 3D
---------

Similar to a Vector2D, Vector3D will also have an x and y component, plus the z component , which are both floats
representing the vector's coordinates. It also has an optional name argument.

Here's how you can instantiate the Vector3D class:

.. plot::
   :include-source:
   :align: center

   import volmdlr
   import matplotlib.pyplot as plt
   import mplcyberpunk
   plt.style.use("cyberpunk")

   vector = volmdlr.Vector3D(1.0, 1.0, 1.0)
   vector.plot(color='orange')

How to create a Vector
======================

Point 2D
--------

Point 3D
--------


How to create a Frame3D
=======================

Curves
******

How to create a Line
====================

How to create a Cicle
=====================

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
