=========================
I N T E R S E C T I O N S
=========================

Here You will find various explanations and examples for many intersection operations,
and between many kinds of objects, like edges with edges, edges with surfaces, surface with a shell,
shell with another shell and much more.


Intersection between curves:
****************************


To calculate the intersections between any two curves, you can call the `curve_intersections` methods.
It will work for any edges in curves.py module.

For example, to calculate the intersections between a Ellipse2D and Circle2D you can do as follows:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import edges, curves
    from volmdlr.core import EdgeStyle

    circle = curves.Circle2D(frame=volmdlr.Frame2D(origin=volmdlr.Point2D(0.0, 0.0), u=volmdlr.Vector2D(1.0, -2.2204460492503126e-16), v=volmdlr.Vector2D(-2.4980018054066027e-16, 1.0)), radius=1)
    ellipse2d = curves.Ellipse2D(frame=volmdlr.Frame2D(origin=volmdlr.Point2D(1.6329931618544102, -1.1547005383785152), u=volmdlr.Vector2D(-0.8164965809277259, 0.5773502691896258), v=volmdlr.Vector2D(0.5773502691896258, 0.8164965809277259)),
                                 major_axis=1.7320508075677725, minor_axis=0.9999999999993615)

    intersections = circle.curve_intersections(ellipse2d)

    ax = ellipse2d.plot()
    circle.plot(ax, EdgeStyle('r'))
    for p in intersections:
        p.plot(ax, 'b')


Intersections betweeen edges
****************************

Edges are generally objects having a curve as its base tragectoty and a start and end points on that curve.
We Have edges like LineSegment, Arc, FullArc, ArcEllipse, FullArcEllipse, BSplineCurve, both in 2-D and in 3-D.

To calculate the intersections between any two edges, you can call the `intersections` methods.
It will work for any edges in edges.py module.

For example, to calculate the intersections between a BSplineCurve2D and Arc2D you can do as follows:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import curves, edges
    from geomdl import utilities
    from volmdlr.core import EdgeStyle

    degree = 3
    points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
    knotvector = utilities.generate_knot_vector(degree, len(points))
    knot_multiplicity = [1] * len(knotvector)
    bspline1 = edges.BSplineCurve2D(degree, points, knot_multiplicity, knotvector, None, False)

    arc = edges.Arc2D.from_3_points(volmdlr.Point2D(0, 0.3), volmdlr.Point2D(1, -0.3), volmdlr.Point2D(2, 2))
    intersections = bspline1.intersections(arc)
    ax = bspline1.plot()
    arc.plot(ax, EdgeStyle('r'))
    for i in intersections:
        i.plot(ax, 'c')


Another example in 3D: BSplineCurve3D and LineSegment3D:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import edges
    from volmdlr.core import EdgeStyle
    from volmdlr.models.edges import bspline_curve3d

    #Here we are using a bspline curve previously defined as a model curve.
    bs3d = bspline_curve3d()

    #take two points of the bspline to use them to defined an intersecting line segment.
    pt1, point2 = bs3d.points[30], bs3d.points[85]
    vec = point2 - pt1
    vec = vec.to_vector()
    lineseg = edges.LineSegment3D(pt1 - vec*0.5, point2 + vec * 0.5)

    #Search for intersections
    intersections = bs3d.intersections(lineseg)

    #plot results, the point on blue are the intersection points.
    ax = bs3d.plot()
    lineseg.plot(ax, EdgeStyle('r'))
    for intersection in intersections:
        intersection.plot(ax, 'b')



Intersections betweeen an edge and a curve
******************************************



Intersections betweeen a wire/contour and an edge
*************************************************

If you ever need to calcule the intersections of a wire/contour with any edge, you can use the `edge_intersections` method.

Example:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr.models.contours import contour2_unittest
    from volmdlr.models.edges import bspline1
    from volmdlr.core import EdgeStyle

    # Here we are going to use a model contour2d and also a model bsplinecurve2d.

    contour2d = contour2_unittest.copy()
    contour2d = contour2d.translation(volmdlr.Vector2D(1, 0.0))

    #search for intersections betweeen a contour2d and an edge.
    edge_intersections = contour2d.edge_intersections(bspline1)

    #plot
    ax = contour2d.plot()
    bspline1.plot(ax, EdgeStyle('r'))
    for intersection in edge_intersections:
        intersection.plot(ax, 'g')


Intersections between surface and edges
***************************************

In the same way, it is also possible to calculate intersections between a Surface3D and an edge of
any type using the same command: `edge_intersections` . #todo: structure for method *edge_structure*.

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import edges, primitives3d, wires, surfaces
    from volmdlr.core import EdgeStyle
    from volmdlr.utils.common_operations import random_color


    spherical_surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
    ax = spherical_surface3d.plot(color='b')

    linesegment = edges.LineSegment3D(volmdlr.Point3D(-0.8, -0.8, -0.8), volmdlr.Point3D(0.8, 0.8, 0.8))
    linesegment_intersections = spherical_surface3d.linesegment_intersections(linesegment)
    linesegment.plot(ax, EdgeStyle('g'))
    for p in linesegment_intersections:
        p.plot(ax, random_color())

Intersections between two surfaces
**********************************

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import faces, surfaces, edges
    from volmdlr.core import EdgeStyle

    R = 0.15
    cylindricalsurface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, R)

    plane = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.O3D, volmdlr.Vector3D(0.0, 0.7071067811865475, -0.7071067811865476), volmdlr.Vector3D(1.0, 0.0, 0.0), volmdlr.Vector3D(0.0, -0.7071067811865476, -0.7071067811865475)))
    intersections = cylindricalsurface.plane_intersections(plane)
    ax = cylindricalsurface.plot(edge_style=EdgeStyle('k'))
    plane.plot(ax, edge_style=EdgeStyle('b'))
    for intersection in intersections:
        intersection.plot(ax, EdgeStyle('r'))


Intersections betweeen surfaces and faces
*****************************************

The intersections between a Face3D and a Surface3D are also possible, and you can do it by also calling the `surface_intersection` method.

Example:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import faces, surfaces, edges
    from volmdlr.core import EdgeStyle

Intersections between two faces
*******************************

If you have any two faces intersecting, you can search these intersections by calling `face_intersections` method.

Example:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import surfaces, faces
    from volmdlr.core import EdgeStyle

    face1 = faces.PlaneFace3D.from_surface_rectangular_cut(surfaces.Plane3D(volmdlr.OXYZ), -1, 1, -1, 1)
    face2 = faces.PlaneFace3D.from_surface_rectangular_cut(surfaces.Plane3D(volmdlr.OYZX), -1, 1, -1, 1)
    face_intersections = face1.face_intersections(face2)
    ax = face1.plot()
    face2.plot(ax, 'r')
    for intersection in face_intersections:
        intersection.plot(ax, EdgeStyle('b'))



Intersection between a face and an edge
***************************************

Edges can intersect faces too. The method `edge_intersections` will do all the work, and it will return the intersections between two objects.

Example:

.. plot::
    :include-source:
    :align: center

    import volmdlr
    from volmdlr import faces, surfaces, edges
    from volmdlr.core import EdgeStyle

    face1 = faces.PlaneFace3D.from_surface_rectangular_cut(surfaces.Plane3D(volmdlr.OXYZ), -1, 1, -1, 1)
    lineseg = edges.LineSegment3D(volmdlr.Point3D(-.5, -.5, -.5), volmdlr.Point3D(.5, .5, .5))
    linseg_intersections = face1.linesegment_intersections(lineseg)

    ax = face1.plot()
    lineseg.plot(ax, EdgeStyle('g'))
    for ls in linseg_intersections:
        ls.plot(ax, 'y')


Intersections between two shells
********************************

.. toctree::
   :maxdepth: 1

   Shells Boolean Operations <boolean_operations>
