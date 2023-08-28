=========
DISTANCES
=========

This section presents how to compute distances between `volmdlr` objects.

Distance between a point and another object
-------------------------------------------------

To compute the distance between a point and **any** other `volmdlr` object, the `point_distance` method proves
sufficient.
In cases where an object lacks the `point_distance` method, attempting to use it will result in an error being raised.

Here are some examples:

Distance between two ``Point2D``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :include-source:
    :align: center

    from volmdlr import Point2D

    # Define the points
    point1 = Point2D(0.0, 0.0)
    point2 = Point2D(1.0, 1.0)

    # Compute and print the distance
    distance = point1.point_distance(point2)
    print(f"Distance between point1 and point2: {distance}")

    # Plot the points
    ax = point1.plot(color="k")
    point2.plot(ax, color="c")

Expected output:

.. code-block:: none

    Distance between point1 and point2: 1.4142135623730951


Distance between ``Point2D`` and a ``LineSegment2D``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :include-source:
    :align: center

    from volmdlr import Point2D
    from volmdlr.core import EdgeStyle
    from volmdlr.edges import LineSegment2D

    # Define a point and a line segment
    point = Point2D(0.0, 0.0)
    line_segment = LineSegment2D(Point2D(1.0, 1.0), Point2D(2.0, -1.0))

    # Compute and print the distance
    distance, closest_point = line_segment.point_distance(point, return_other_point=True)
    print(f"Distance between point and line_segment: {distance}")

    # Plot the points and the line segment
    ax = point.plot()
    line_segment.plot(ax, EdgeStyle(color="c"))
    closest_point.plot(ax, "g")

Expected output:

.. code-block:: none

    Distance between point and line_segment: 1.3416407864998738


Distance between a ``Point3D`` and a ``PlaneFace3D``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :include-source:
    :align: center

    from volmdlr import OXYZ, Point3D
    from volmdlr.faces import PlaneFace3D
    from volmdlr.surfaces import Plane3D

    # Define a point and a plane face
    point = Point3D(0.2, -0.3, 0.7)
    plane = Plane3D(OXYZ)
    plane_face = PlaneFace3D.from_surface_rectangular_cut(plane, -1, 1, -1, 1)

    # Compute and print the distance
    distance, closest_point = plane_face.point_distance(point, True)
    print(f"Distance between point and plane_face: {distance}")

    # Plot the points and the plane face
    ax = plane_face.plot()
    point.plot(ax, "r")
    closest_point.plot(ax, "b")


Expected output:

.. code-block:: none

    Distance between point and plane_face: 0.7

Distance between two edges
--------------------------

For calculating the distance between any two edges, employing the `minimum_distance` method is sufficient.

Here are some examples:

Distance beteween a ``BSplineCurve2D`` and a ``LineSegment2D``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :include-source:
    :align: center

    from volmdlr import Point2D
    from volmdlr.core import EdgeStyle
    from volmdlr.edges import BSplineCurve2D, LineSegment2D
    from volmdlr.nurbs.helpers import generate_knot_vector

    # Define LineSegment2D
    line_segment = LineSegment2D(Point2D(1, 0.5), Point2D(3, 1))

    # Define a BSplineCurve2D
    DEGREE = 3
    points = [Point2D(0, 0), Point2D(1, 1), Point2D(2, -1), Point2D(3, 0)]
    knot_vector = generate_knot_vector(DEGREE, len(points))
    knot_multiplicity = [1] * len(knot_vector)
    b_spline_curve = BSplineCurve2D(DEGREE, points, knot_multiplicity, knot_vector, None)

    # Compute and print the minimum distance
    distance, point1, point2 = b_spline_curve.minimum_distance(line_segment, return_points=True)
    print(f"Minimum distance between b_spline_curve and line_segment: {distance}")

    # Plot the line_segment, the b_spline_curve and the points
    ax = b_spline_curve.plot()
    line_segment.plot(ax, EdgeStyle("g"))
    point1.plot(ax, "r")
    point2.plot(ax, "b")


Expected output:

.. code-block:: none

    Minimum distance between b_spline_curve and line_segment: 0.2655463988751082

Distance beteween an ``Arc3D`` and a ``LineSegment2D``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. plot::
    :include-source:
    :align: center

    from volmdlr import O3D, Frame3D, Point3D, Vector3D
    from volmdlr.core import EdgeStyle
    from volmdlr.curves import Circle3D, Ellipse3D
    from volmdlr.edges import Arc3D, ArcEllipse3D

    vector1 = Vector3D(1, 1, 1)
    vector1 = vector1.unit_vector()
    vector2 = vector1.deterministic_unit_normal_vector()
    vector3 = vector1.cross(vector2)

    # Define the Arc3D
    circle3d = Circle3D(Frame3D(O3D, vector1, vector2, vector3), 1)
    arc_3d = Arc3D(
        circle=circle3d,
        start=Point3D(0.57, 0.57, 0.57),
        end=Point3D(-0.98, -0.11, -0.11),
    )

    # Define the ArcEllipse3D
    ellipse3d = Ellipse3D(2, 1, Frame3D(Point3D(1, 0, 1), vector3, vector1, vector2))
    arc_ellipse_3d = ArcEllipse3D(
        ellipse=ellipse3d,
        start=Point3D(0.42, -0.57, 0.42),
        end=Point3D(1.57, 0.57, 1.57),
    )

    # Compute and print the minimum distance
    distance, point1, point2 = arc_3d.minimum_distance(arc_ellipse_3d, return_points=True)
    print(f"Minimum distance between arc_3d and arc_ellipse_3d: {distance}")

    # Plot the arc_3d, the arc_ellipse_3d and the points
    ax = arc_3d.plot()
    arc_ellipse_3d.plot(ax, EdgeStyle("g"))
    point1.plot(ax, "r")
    point2.plot(ax, "b")


Expected output:

.. code-block:: none

    Minimum distance between arc_3d and arc_ellipse_3d: 0.8102555248814264

Distances betweeen two faces
----------------------------

Calculating the distance between two faces is also feasible. To achieve this, rely on the `face_minimum_distance` method.
This method is universally applicable for any pair of faces.

.. plot::
    :include-source:
    :align: center

    from volmdlr import OXYZ, TWO_PI, Frame3D, Point3D, Vector3D
    from volmdlr.faces import CylindricalFace3D
    from volmdlr.surfaces import CylindricalSurface3D

    # Define two cylindrical faces
    R = 0.15
    cylindrical_surface = CylindricalSurface3D(OXYZ, R)
    cylindrical_face1 = CylindricalFace3D.from_surface_rectangular_cut(
        cylindrical_surface, 0, TWO_PI, -0.25, 0.25
    )

    u_vector = Vector3D(-1, -1, -1)
    u_vector = u_vector.unit_vector()
    v_vector = u_vector.deterministic_unit_normal_vector()
    w_vector = u_vector.cross(v_vector)
    cylindrical_face2 = cylindrical_face1.frame_mapping(
        Frame3D(Point3D(-0.5, 0.5, -0.1), u_vector, v_vector, w_vector), "new"
    )

    # Compute and print the minimum distance
    minimum_distance, point1, point2 = cylindrical_face1.face_minimum_distance(cylindrical_face2, True)
    print(f"Minimum distance between cylindrical_face1 and cylindrical_face2: {minimum_distance}")

    # Plot the cylindrical_face1, the cylindrical_face2 and the points
    ax = cylindrical_face1.plot()
    cylindrical_face2.plot(ax, "r")
    point1.plot(ax, "y")
    point2.plot(ax, "b")

Expected output:

.. code-block:: none

    Minimum distance between cylindrical_face1 and cylindrical_face2: 0.30994820373017934

Distances betweeen two shells
-----------------------------

Similarly, you can determine the distance between two shells, as illustrated in the following example:

.. plot::
    :include-source:
    :align: center

    import math

    import volmdlr
    from volmdlr.faces import PlaneFace3D, Triangle3D
    from volmdlr.shells import ClosedShell3D
    from volmdlr.surfaces import Plane3D, Surface2D
    from volmdlr.wires import ClosedPolygon3D

    # Create the contours of the faces of a random shape
    polygon1_vol1 = ClosedPolygon3D(
        [
            volmdlr.Point3D(-0.1, -0.05, 0),
            volmdlr.Point3D(-0.15, 0.1, 0),
            volmdlr.Point3D(0.05, 0.2, 0),
            volmdlr.Point3D(0.12, 0.15, 0),
            volmdlr.Point3D(0.1, -0.02, 0),
        ]
    )
    polygon2_vol1 = polygon1_vol1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi).translation(
        0.2 * volmdlr.Z3D
    )
    polygon3_vol1 = polygon2_vol1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 8).translation(
        0.1 * (volmdlr.Z3D + volmdlr.X3D + volmdlr.Y3D)
    )

    # Create the faces of the shape
    faces = [
        Triangle3D(*points)
        for points in polygon1_vol1.sewing(polygon2_vol1, volmdlr.X3D, volmdlr.Y3D)
    ] + [
        Triangle3D(*points)
        for points in polygon2_vol1.sewing(polygon3_vol1, volmdlr.X3D, volmdlr.Y3D)
    ]

    # Create the top and bottom faces of the shape
    bottom_surface3d = Plane3D.from_plane_vectors(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
    bottom_surface2d = Surface2D(polygon1_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D), [])
    top_surface3d = Plane3D.from_plane_vectors(0.3 * volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D)
    top_surface2d = Surface2D(polygon3_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D), [])
    bottom_face = PlaneFace3D(bottom_surface3d, bottom_surface2d)
    top_face = PlaneFace3D(top_surface3d, top_surface2d)
    faces += [bottom_face, top_face]

    # Create a closed shell from these faces
    shell1 = ClosedShell3D(faces)

    # Create a second shell from the first one, by rotating and translating it.
    shell2 = shell1.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 5)
    shell2 = shell2.translation(volmdlr.Vector3D(0.5, 0.5, 0.5))

    # Compute and print the distance
    minimum_distance, point1, point2 = shell1.minimum_distance(shell2, True)
    print(f"Minimum distance between shell1 and shell2: {minimum_distance}")

    # Plot the shells
    ax = shell1.plot()
    shell2.plot(ax, "r")
    point1.plot(ax, "b")
    point2.plot(ax, "g")

Expected output:

.. code-block:: none

    Minimum distance between shell1 and shell2: 0.3374259086917476
