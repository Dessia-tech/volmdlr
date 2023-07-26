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

.. code-block:: python

   import volmdlr

   vector = volmdlr.Vector2D(1.0, 1.0)
   vector.plot(color='r')

Vector 3D
---------

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
