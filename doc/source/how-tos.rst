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

           vector = volmdlr.Vector2D(1.0, 1.0)
           vector.plot(color='orange')

    .. grid-item-card::  Vector3D

        .. plot::
           :include-source:
           :align: center

           import volmdlr

           vector = volmdlr.Vector3D(1.0, 1.0, 1.0)
           vector.plot(color='orange')

How to create a Point
=====================

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

   point3d = volmdlr.Point3D(x=1.0, y=2.0, z=3.0, name='vector name')
   print('point3d:', point3d)
   >>> point2d: Point3D: [1.0, 2.0, 3.0]

How to create a Frame3D
=======================
The given class Frame3D represents a 3D frame. It defines a frame by specifying its origin point and
three basis vectors (u, v, and w) that determine the orientation of the frame.

Exmaple:

.. plot::
   :include-source:
   :align: center

   import volmdlr

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

           point1 = volmdlr.Point2D(1.0, 1.0)
           point2 = volmdlr.Point2D(-2.0, -3.0)
           line2d = curves.Line2D(point1, point2, name='line2d_name_is_optional')
           line2d.plot(edge_style=EdgeStyle('orange'))

    .. grid-item-card::  Line3D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import curves
           from volmdlr.core import EdgeStyle

           point1 = volmdlr.Point3D(1.0, 1.0, 1.0)
           point2 = volmdlr.Point3D(-2.0, -3.0, -1.0)
           line3d = curves.Line3D(point1, point2, name='line3d_name_is_optional')
           line3d.plot(edge_style=EdgeStyle('orange'))

How to create a Circle
======================

Circle2D
--------

The circle2d  is defined by its center point (Point2D) and a radius (float),
and it also accepts an optional name for identification.

.. grid-item-card::

    .. plot::
       :include-source:
       :align: center

       import volmdlr
       from volmdlr import curves
       from volmdlr.core import EdgeStyle

       center2d = volmdlr.Point2D(0.0, 0.0)
       circle2d = curves.Circle2D(center=center2d, radius=1, name='optional_circle_name')
       circle2d.plot(edge_style=EdgeStyle('orange'))

Circle3D
--------

The circle is defined by a Frame3D object that includes information about the center and orientation of the
circle in 3D space, along with a radius (float). The frame's u and v vectors define the plane in which the
circle lies, and w represents the normal vector to the plane.

.. grid-item-card::

    .. plot::
       :include-source:
       :align: center

       import volmdlr
       from volmdlr import curves
       from volmdlr.core import EdgeStyle

       center3D = volmdlr.Point3D(0.0, 0.0, 0.0)
       u_vector = volmdlr.Vector3D(1.0, 0.0, 0.0)
       v_vector = volmdlr.Vector3D(0.0, 1.0, 0.0)
       w_vector = volmdlr.Vector3D(0.0, 0.0, 1.0)
       frame3d = volmdlr.Frame3D(center3D, u_vector, v_vector, w_vector)
       circle3d = curves.Circle3D(frame=frame3d, radius=1, name='optional_circle_name')
       circle3d.plot(edge_style=EdgeStyle('orange'))

How to create an Ellipse
========================

An ellipse in defined by three arguments: a major axis (A), e minor axis (B) and a Frame (2D or 3D).

.. grid:: 1

    .. grid-item-card::  Ellipse2D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import curves
           from volmdlr.core import EdgeStyle

           u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
           v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
           ellipse2d = curves.Ellipse2D(major_axis=2, minor_axis=1, frame=volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector))
           ellipse2d.plot(edge_style=EdgeStyle(color='orange'))

    .. grid-item-card::  Ellipse3D

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import curves
           from volmdlr.core import EdgeStyle

           vector1 = volmdlr.Vector3D(1, 1, 1)
           vector1 = vector1.unit_vector()
           vector2 = vector1.deterministic_unit_normal_vector()
           vector3 = vector1.cross(vector2)
           frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
           ellipse3d = curves.Ellipse3D(major_axis=2, minor_axis=1, frame=frame)
           ellipse3d.plot(edge_style=EdgeStyle('orange'))

Edges
*****

How to create a LineSegment
===========================

To instanciate a LineSegment in 2D or 3D, you will need to provide two mandatory arguments, the start and end points.
Additionally you will have two optional arguments: an infinite line which the line segment lies on and a name argument.

LineSegment2D
-------------
.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges
           from volmdlr.core import EdgeStyle

           start_point = volmdlr.Point2D(1.0, 1.0)
           end_point = volmdlr.Point2D(3.0, 4.0)
           linesegment2d = edges.LineSegment2D(start=start_point, end=end_point, line=None, name='linesegment\'s name')
           linesegment2d.plot(edge_style=EdgeStyle(color='orange'))

LineSegment3D
-------------

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges
           from volmdlr.core import EdgeStyle

           start_point = volmdlr.Point3D(1.0, 1.0, 1.0)
           end_point = volmdlr.Point3D(3.0, 4.0, 6.0)
           linesegment3d = edges.LineSegment3D(start=start_point, end=end_point, line=None, name='linesegment\'s name')
           linesegment3d.plot(edge_style=EdgeStyle(color='orange'))


How to create an Arc
====================

Arc2D
-----

An Arc2D is defined by the base circle curve along with a start and end points.
There is also a boolean is_trigo argument that defines if the arc is in the trigo-wise direction or not and a last and optional name argument.

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges, curves
           from volmdlr.core import EdgeStyle

           circle2d = curves.Circle2D(volmdlr.O2D, 1)
           arc2d = edges.Arc2D(circle2d, volmdlr.Point2D(-1, 0), volmdlr.Point2D(1, 0), True)
           ax = arc2d.plot(edge_style=EdgeStyle('orange'))
           ax.set_aspect('equal')


Arc3D
-----

Just as the Arc2D, Arc3D is defined by the base circle curve along with a start and end points. There is also an optional name argument.

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges, curves
           from volmdlr.core import EdgeStyle

           vector1 = volmdlr.Vector3D(1, 1, 1)
           vector1 = vector1.unit_vector()
           vector2 = vector1.deterministic_unit_normal_vector()
           vector3 = vector1.cross(vector2)
           frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
           circle3d = curves.Circle3D(frame, 1)
           arc3d = edges.Arc3D(circle3d, start=volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                       end=volmdlr.Point3D(-0.9855985596534886, -0.11957315586905026, -0.11957315586905026))
           ax = arc3d.plot(edge_style=EdgeStyle('orange'))


How to create an ArcEllipse
===========================

Both ArcEllipse2D and ArcEllipse3D require a base Ellipse curve along with a start end end points.
There also an optional name argument.

ArcEllipse2D
------------

Object's descrition

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges, curves
           from volmdlr.core import EdgeStyle

           u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
           v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
           ellipse2d = curves.Ellipse2D(2, 1, volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector))
           u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
           v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
           ellipse2d = curves.Ellipse2D(2, 1, volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector))
           arc_ellipse2d = edges.ArcEllipse2D(ellipse2d, start=volmdlr.Point2D(0.5, 1.5), end=volmdlr.Point2D(1.5, 0.5))
           arc_ellipse2d.plot(edge_style=EdgeStyle('orange'))



ArcEllipse3D
------------

Object's descrition

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges, curves
           from volmdlr.core import EdgeStyle

           vector1 = volmdlr.Vector3D(1, 1, 1)
           vector1 = vector1.unit_vector()
           vector2 = vector1.deterministic_unit_normal_vector()
           vector3 = vector1.cross(vector2)
           frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
           start_point = volmdlr.Point3D(0.2391463117381003, 1.1051717155225391, 1.1051717155225391)
           end_point = volmdlr.Point3D(-1.393846850117352, -0.5278214463329132, -0.5278214463329132)
           ellipse3d = curves.Ellipse3D(2, 1, frame)
           arc_ellipse3d = edges.ArcEllipse3D(ellipse3d, start=start_point, end=end_point)
           arc_ellipse3d.plot(edge_style=EdgeStyle('orange'))



How to create a BSplineCurve
============================

To instanciate a BsplineCurve 2D or 3D, we must provide the necessary parameters, such as the degree, control points,
knot multiplicities, knot vector, and optional weights and name.

BSplineCurve2D
--------------

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges
           from volmdlr.core import EdgeStyle
           from geomdl import utilities

           DEGREE = 3
           points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
           knotvector = utilities.generate_knot_vector(DEGREE, len(points))
           knot_multiplicity = [1] * len(knotvector)
           bspline1 = edges.BSplineCurve2D(DEGREE, points, knot_multiplicity, knotvector, None, False)
           bspline1.plot(edge_style=EdgeStyle('orange'))


BSplineCurve3D
--------------

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges
           from volmdlr.core import EdgeStyle

           degree = 5
           control_points = [volmdlr.Point3D(0, 3, 0),
                             volmdlr.Point3D(3, 2, 1),
                             volmdlr.Point3D(5, -1, 4),
                             volmdlr.Point3D(5, -4, 0),
                             volmdlr.Point3D(-1, -2, -3),
                             volmdlr.Point3D(-3, 4, 1)]
           knots = [0.0, 1.0]
           knot_multiplicities = [6, 6]
           weights = None  # [1, 2, 1, 2, 1, 2]
           bspline_curve3d = edges.BSplineCurve3D(degree=degree, control_points=control_points,
                                           knot_multiplicities=knot_multiplicities,
                                           knots=knots,
                                           weights=weights,
                                           periodic=False,
                                           name='B Spline Curve 3D 1')
           bspline_curve3d.plot(edge_style=EdgeStyle('orange'))


Wires
*****

How to create a Wire
====================

A wire is an object composed of a list of primitives that does not form a closed and an optional name. This primitives list can contain any set of edges following each other.

Wire2D
------

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import wires, edges
           from volmdlr.core import EdgeStyle

           line_segment1 = edges.LineSegment2D(volmdlr.Point2D(1, -1), volmdlr.Point2D(1.5, 1))
           arc = edges.Arc2D.from_3_points(volmdlr.Point2D(1.5, 1), volmdlr.Point2D(1.3, 1.5), volmdlr.Point2D(0.5, 1.5))
           points2d = [volmdlr.Point2D(-1, 1), volmdlr.Point2D(2, 2), volmdlr.Point2D(-2, -2), volmdlr.Point2D(1, -1)]
           bspline = edges.BSplineCurve2D(3, points2d, knot_multiplicities=[4, 4], knots=[0.0, 1.0])
           wire2d = wires.Wire2D([bspline, line_segment1, arc])
           wire2d.plot(edge_style=EdgeStyle('orange'))


Wire3D
------

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges, wires
           from volmdlr.core import EdgeStyle

           degree = 5
           control_points = [volmdlr.Point3D(0, 3, 0),
                            volmdlr.Point3D(3, 2, 1),
                            volmdlr.Point3D(5, -1, 4),
                            volmdlr.Point3D(5, -4, 0),
                            volmdlr.Point3D(-1, -2, -3),
                            volmdlr.Point3D(-3, 4, 1)]
           knots = [0.0, 1.0]
           knot_multiplicities = [6, 6]
           weights = None  # [1, 2, 1, 2, 1, 2]
           bspline_curve3d = edges.BSplineCurve3D(degree=degree, control_points=control_points,
                                          knot_multiplicities=knot_multiplicities,
                                          knots=knots,
                                          weights=weights,
                                          periodic=False,
                                          name='B Spline Curve 3D 1')
           lineseg1 = edges.LineSegment3D(volmdlr.Point3D(3, 3, 2), bspline_curve3d.start)
           lineseg2 = edges.LineSegment3D(bspline_curve3d.end, volmdlr.Point3D(-3, -3, 0))
           wire3d = wires.Wire3D([lineseg1, bspline_curve3d, lineseg2])
           wire3d.plot(edge_style=EdgeStyle('orange'))

How to create a Contour
=======================

As a  wire, A contour is also an object composed of a list of primitives but now it forms a closed loop.
The optional name argument is always present. The primitives list can contain any set of edges following each other.

Contour2D
---------

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges, wires
           from volmdlr.core import EdgeStyle

           line_segment1 = edges.LineSegment2D(volmdlr.Point2D(1, -1), volmdlr.Point2D(1.5, 1))
           line_segment2 = edges.LineSegment2D(volmdlr.Point2D(0.5, 1.5), volmdlr.Point2D(-2, 1))
           line_segment3 = edges.LineSegment2D(volmdlr.Point2D(-2, 1), volmdlr.Point2D(-2, 0.7))
           line_segment4 = edges.LineSegment2D(volmdlr.Point2D(-2, 0.7), volmdlr.Point2D(-1, 1))
           arc = edges.Arc2D.from_3_points(volmdlr.Point2D(1.5, 1), volmdlr.Point2D(1.3, 1.5), volmdlr.Point2D(0.5, 1.5))
           points2d = [volmdlr.Point2D(-1, 1), volmdlr.Point2D(2, 2), volmdlr.Point2D(-2, -2), volmdlr.Point2D(1, -1)]
           bspline = edges.BSplineCurve2D(3, points2d, knot_multiplicities=[4, 4], knots=[0.0, 1.0])
           wire2d = wires.Wire2D([bspline, line_segment1, arc, line_segment2, line_segment3, line_segment4])
           wire2d.plot(edge_style=EdgeStyle('orange'))

Contour3D
---------

.. grid:: 1

    .. grid-item-card::

        .. plot::
           :include-source:
           :align: center

           import volmdlr
           from volmdlr import edges, wires
           from volmdlr.core import EdgeStyle

           degree = 5
           control_points = [volmdlr.Point3D(0, 3, 0),
                            volmdlr.Point3D(3, 2, 1),
                            volmdlr.Point3D(5, -1, 4),
                            volmdlr.Point3D(5, -4, 0),
                            volmdlr.Point3D(-1, -2, -3),
                            volmdlr.Point3D(-3, 4, 1)]
           knots = [0.0, 1.0]
           knot_multiplicities = [6, 6]
           weights = None  # [1, 2, 1, 2, 1, 2]
           bspline_curve3d = edges.BSplineCurve3D(degree=degree, control_points=control_points,
                                          knot_multiplicities=knot_multiplicities,
                                          knots=knots,
                                          weights=weights,
                                          periodic=False,
                                          name='B Spline Curve 3D 1')
           lineseg1 = edges.LineSegment3D(volmdlr.Point3D(3, 3, 2), bspline_curve3d.start)
           lineseg2 = edges.LineSegment3D(bspline_curve3d.end, volmdlr.Point3D(-3, -3, 0))
           arc = edges.Arc3D.from_3_points(volmdlr.Point3D(-3, -3, 0), volmdlr.Point3D(6.324555320336761, -5.692099788303083, -0.8973665961010275), volmdlr.Point3D(3, 3, 2))
           wire3d = wires.Wire3D([lineseg1, bspline_curve3d, lineseg2, arc])
           wire3d.plot(edge_style=EdgeStyle('orange'))

Faces
*****

PlaneFace3D
===========

To create a `PlaneFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `Plane3D`, which is constructed using a `Frame3D` and an optional `name` parameter.

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `PlaneFace3D`.

.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

            import volmdlr
            from volmdlr import edges, curves, surfaces, wires, faces
            from volmdlr.core import EdgeStyle

            surface3d = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(0.0, 0.0, 0.0), volmdlr.Vector3D(1.0, 0.0, 0.0),
                                                        volmdlr.Vector3D(0.0, 1.0, 0.0), volmdlr.Vector3D(0.0, 0.0, 1.0)))

            outer_contour2d = wires.Contour2D.from_points(points=[volmdlr.Point2D(0., 0.), volmdlr.Point2D(2, 0),
                                                                 volmdlr.Point2D(2, 2), volmdlr.Point2D(1, 2),
                                                                 volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)])
            inner_contours2d = []
            surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=inner_contours2d)

            plane_face = faces.PlaneFace3D(surface3d=surface3d, surface2d=surface2d)

            plane_face.babylonjs()

        .. figure:: ../source/_static/index-images/planeface3d.png


Triangle3D
==========

A Triangle3D receives three mandatory arguments: The three vertices points of the triaangle, along with a last optional name argument.

.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle

           triangle3d = faces.Triangle3D(volmdlr.Point3D(0., 0., 1.0), volmdlr.Point3D(2, 0, 0.2), volmdlr.Point3D(2, 2, 3.0))
           triangle3d.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/triangle3d.png

CylindricalFace3D
=================

To create a `CylindricalFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `CylindricalSurface3D`, which is constructed using a `Frame3D`, a float value for the cylinder radius and an optional `name` parameter.

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `CylindricalFace3D`.

.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle

           vector1 = volmdlr.Vector3D(1, 1, 1)
           vector1 = vector1.unit_vector()
           vector2 = vector1.deterministic_unit_normal_vector()
           vector3 = vector1.cross(vector2)
           frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)

           surface3d = surfaces.CylindricalSurface3D(frame, 1)

           outer_contour2d = wires.Contour2D.from_points(points=[volmdlr.Point2D(0., 0.), volmdlr.Point2D(4, 0),
                                                                            volmdlr.Point2D(4, 4), volmdlr.Point2D(2, 4),
                                                                            volmdlr.Point2D(2, 2), volmdlr.Point2D(0, 2)])
           surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=[])

           face3d = faces.CylindricalFace3D(surface3d, surface2d)

           face3d.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/cylindricalface3d.png

ToroidalFace3D
==============

To create a `ToroidalFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `ToroidalSurface3D`, which is constructed using three main arguments:

    - `Frame3D`: the three dimensional frame where the toroidal face is at.
    - tore_radius: The distance from the center of the torus to the center of the tube (the larger radius).
    - small_radius: The radius of the tube (the smaller radius).

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `ToroidalFace3D`.

.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle


           surface3d = surfaces.ToroidalSurface3D(volmdlr.OXYZ, tore_radius=0.2, small_radius=0.03, name='optional_toroidalsurface3d\'s_name')

           points = [volmdlr.Point2D(-1.0, 0), volmdlr.Point2D(1, 0), volmdlr.Point2D(1, 3.5), volmdlr.Point2D(-1, 3.5)]
           outer_contour2d = wires.Contour2D.from_points(points=points)
           surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=[])

           toroidal_face3d = faces.ToroidalFace3D(surface3d, surface2d)

           toroidal_face3d.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/toroidalface3d.png

ConicalFace3D
=============

To create a `ConicalFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `ConicalSurface3D`, which is constructed using two main arguments:

    - `Frame3D`: the three dimensional frame where the conical face is at. The frame.w is the cone's axis
    - semi_angle: The semi-angle of a cone refers to the angle between the central axis of the cone and a line connecting the apex (top) of the cone to a point on the base.

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `ConicalFace3D`.

.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle


           surface3d = surfaces.ConicalSurface3D(volmdlr.OXYZ, semi_angle=0.2, name='optional_conicalsurface3d\'s_name')

           points = [volmdlr.Point2D(-1.0, 0.0), volmdlr.Point2D(3.0, 0.0), volmdlr.Point2D(3.0, 4.0), volmdlr.Point2D(-1.0, 4.0)]
           outer_contour2d = wires.Contour2D.from_points(points=points)
           surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=[])

           toroidal_face3d = faces.ConicalFace3D(surface3d, surface2d)

           toroidal_face3d.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/conicalface3d.png

SphericalFace3D
===============

To create a `SphericalFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `SphericalSurface3D`, which is constructed using two main arguments:

    * `Frame3D`: the three dimensional frame where the spherical face is at. The frame.origin is the spheres' center.
    * radius: the radius of the sphere.

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `ConicalFace3D`.


.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle


           surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, radius=0.2, name='optional_sphericalsurface3d\'s_name')

           points = [volmdlr.Point2D(0.0, 0.0), volmdlr.Point2D(2.5, 0.0), volmdlr.Point2D(2.5, 1.5), volmdlr.Point2D(0.0, 1.5)]
           outer_contour2d = wires.Contour2D.from_points(points=points)
           surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=[])

           spherical_face3d = faces.SphericalFace3D(surface3d, surface2d)

           spherical_face3d.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/sphericalface3d.png

RuledFace3D
===========

ExtrusionFace3D
===============

To create a `ExtrusionFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `ExtrusionSurface3D`, which is constructed using two main arguments:

    * `edge`: the edge to be estruded.
    * direction: The extrusion direction vector.

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `ExtrusionFace3D`.


.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle


           arc2 = volmdlr.edges.Arc3D(curves.Circle3D(volmdlr.OXYZ, 1), volmdlr.Point3D(1, 0, 0), volmdlr.Point3D(0, 1, 0))
           surface3d = surfaces.ExtrusionSurface3D(edge=arc2, direction=volmdlr.Z3D)

           outer_contour2d = wires.Contour2D.from_points(points=[volmdlr.Point2D(0., 0.), volmdlr.Point2D(1, 0),
                                                                            volmdlr.Point2D(1, 1), volmdlr.Point2D(0.5, 1),
                                                                            volmdlr.Point2D(0.5, 0.5), volmdlr.Point2D(0, 0.5)])
           inner_contours2d = []
           surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=inner_contours2d)

           face = faces.ExtrusionFace3D(surface3d, surface2d)

           face.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/extrusionface3d.png

RevolutionFace3D
================

To create a `RevolutionFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `RevolutionSurface3D`, which is constructed using three main arguments:

    * `edge`: the revolution edge.
    * `axis_point`: revolution's axis point.
    * `axis`: The axis of revolution.

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `RevolutionFace3D`.


.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle

           fullarc = edges.FullArc3D(circle=curves.Circle3D(
                        volmdlr.Frame3D(
                            volmdlr.Point3D(0.003516498393599, -0.01267818173491, 0.0), volmdlr.Vector3D(1.0, 0.0, 0.0),
                            volmdlr.Vector3D(0.0, 1.0, 0.0), volmdlr.Vector3D(0.0, 0.0, 1.0)), radius=0.024102542625267),
                            start_end=volmdlr.Point3D(0.027619041018866, -0.01267818173491, 0.0))

           surface3d = surfaces.RevolutionSurface3D(
                edge=fullarc, axis_point=volmdlr.Point3D(0, 0, 0), axis=volmdlr.Vector3D(0, 1, 0))


           outer_contour2d = wires.Contour2D(primitives=[edges.LineSegment2D(volmdlr.Point2D(0.0, 0.023550776716126855),
                                                                  volmdlr.Point2D(6.283185307179586, 0.023550776716126855)),
                                              edges.LineSegment2D(volmdlr.Point2D(6.283185307179586, 0.023550776716126855),
                                                                  volmdlr.Point2D(6.283185307179586, 0.016162537035284696)),
                                              edges.LineSegment2D(volmdlr.Point2D(6.283185307179586, 0.016162537035284696),
                                                                  volmdlr.Point2D(0.0, 0.016162537035284696)),
                                              edges.LineSegment2D(volmdlr.Point2D(0.0, 0.016162537035284696),
                                                                  volmdlr.Point2D(0.0, 0.023550776716126855))])
           inner_contours2d = []
           surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=inner_contours2d)
           face = faces.RevolutionFace3D(surface3d, surface2d)

           face.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/revolutionface3d.png

BSplineFace3D
=============

To create a `RevolutionFace3D`, you need to provide two arguments: a `surface3d` and a `surface2d`.

1. For the `surface3d`, you must create a `BSplineSurface3D`, for which we have to provide the necessary parameters,
such as the degrees (degree_u and degree_v), control points (instances of Point3D), number of control points
in u and v directions (nb_u and nb_v), knot multiplicities, knot vectors (u_knots and v_knots), optional weights, and name.

2. For the `surface2d`, you instantiate it by providing an outer contour in 2D, which will serve as the outer border of the face. Additionally, you need to provide a list of inner contours in 2D, representing any holes within the face, if applicable. The `surface2d` can also have an optional `name` argument.

Ensure to provide the necessary information for both `surface3d` and `surface2d` to successfully create the `RevolutionFace3D`.


.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces
           from volmdlr.core import EdgeStyle

           control_points = [volmdlr.Point3D(0, 0, 0), volmdlr.Point3D(0.1, 0.02, 0), volmdlr.Point3D(0.2, 0.02, 0),
                             volmdlr.Point3D(0, 0, 0.15), volmdlr.Point3D(0.1, 0.02, 0.15), volmdlr.Point3D(0.2, 0.02, 0.15),
                             volmdlr.Point3D(0, 0, 0.3), volmdlr.Point3D(0.1, 0.021, 0.3), volmdlr.Point3D(0.2, 0.022, 0.3)
                  ]

           surface3d = surfaces.BSplineSurface3D(degree_u=2, degree_v=2, control_points=control_points, nb_u=3, nb_v=3,
                                               u_multiplicities=[1, 2, 2, 1], v_multiplicities=[1, 2, 2, 1],
                                               u_knots=[0.1, 0.3, 0.5, 0.7], v_knots=[0.1, 0.3, 0.5, 0.7])

           outer_contour2d = wires.Contour2D.from_points(points=[volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0),
                                                                 volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)])
           inner_contours2d = []
           surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=inner_contours2d)

           face = faces.RevolutionFace3D(surface3d, surface2d)

           face.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/bsplineface3d.png




Shells
******

A shell is defined as a collection of connected faces. A Shell can a `ClosedShell3D` or an `OpenShell3D`.
it receives as parameters a list of faces (instances of Face3D), optional color, alpha (transparency), name, and a bounding box.

In the example bellow, it is shown the definition of the shell's lateral faces.

.. grid:: 1

    .. grid-item-card::

        .. code-block:: python

           import volmdlr
           from volmdlr import edges, curves, surfaces, wires, faces, shells
           from volmdlr.core import EdgeStyle
           import math

           polygon1_vol1 = wires.ClosedPolygon3D([volmdlr.Point3D(-0.1, -0.05, 0), volmdlr.Point3D(-0.15, 0.1, 0),
                               volmdlr.Point3D(0.05, 0.2, 0), volmdlr.Point3D(0.12, 0.15, 0), volmdlr.Point3D(0.1, -0.02, 0)])

           polygon2_vol1 = polygon1_vol1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi).translation(0.2*volmdlr.Z3D)
           polygon3_vol1 = polygon2_vol1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi/8).translation(0.1*(volmdlr.Z3D+volmdlr.X3D+volmdlr.Y3D))
           faces_ = [faces.Triangle3D(*points)
                   for points in polygon1_vol1.sewing(polygon2_vol1, volmdlr.X3D, volmdlr.Y3D)] + \
                   [faces.Triangle3D(*points)
                   for points in polygon2_vol1.sewing(polygon3_vol1, volmdlr.X3D, volmdlr.Y3D)]


OpenShell3D
===========

.. grid:: 1

    .. grid-item-card::

        With these faces we can instantiate an OpenShell3D:

        .. code-block:: python

           shell1 = shells.OpenShell3D(faces_)
           shell1.babylonjs(dark_mode=True)

        .. figure:: ../source/_static/index-images/openshell3d.png


ClosedShell3D
=============


.. grid:: 1

    .. grid-item-card::

        Then the bottom and top faces can be created so a closedshell3d can be instantiated:

        .. code-block:: python

           bottom_surface3d = surfaces.Plane3D.from_plane_vectors(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
           bottom_surface2d = surfaces.Surface2D(polygon1_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D),[])

           top_surface3d = surfaces.Plane3D.from_plane_vectors(0.3*volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D)
           top_surface2d = surfaces.Surface2D(polygon3_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D),[])

           bottom_face = faces.PlaneFace3D(bottom_surface3d, bottom_surface2d)
           top_face = faces.PlaneFace3D(top_surface3d, top_surface2d)
           faces_ += [bottom_face, top_face]

           shell1 = shells.ClosedShell3D(faces_)
           shell1.babylonjs(dark_mode=True)

    .. figure:: ../source/_static/index-images/closedshell3d.png