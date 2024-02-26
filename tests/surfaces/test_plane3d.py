import math
import unittest
import os
import numpy as np
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import edges, surfaces, curves, wires
from volmdlr.surfaces import Plane3D
from volmdlr.models.edges import bspline_curve3d


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "objects_plane_test")


class TestPlane3D(unittest.TestCase):

    plane1 = surfaces.Plane3D(volmdlr.OXYZ)
    plane2 = surfaces.Plane3D(
        volmdlr.OYZX.translation(volmdlr.Vector3D(1, 2, 1)).rotation(
            volmdlr.Point3D(1, 2, 1), volmdlr.Y3D, math.pi / 4
        )
    )
    plane3 = surfaces.Plane3D(
        volmdlr.OXYZ.translation(volmdlr.Vector3D(0, 0, 1)).rotation(
            volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4
        )
    )
    plane4 = surfaces.Plane3D(volmdlr.OXYZ.rotation(volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4))

    def setUp(self):
        self.point1 = volmdlr.Point3D(0, 0, 0)
        self.point2 = volmdlr.Point3D(1, 0, 0)
        self.point3 = volmdlr.Point3D(0, 1, 0)
        self.point4 = volmdlr.Point3D(0, 0, 1)
        self.vector1 = volmdlr.Vector3D(1, 0, 0)
        self.vector2 = volmdlr.Vector3D(0, 1, 0)
        self.vector3 = volmdlr.Vector3D(0, 0, 1)
        self.plane5 = Plane3D.from_plane_vectors(
            volmdlr.Point3D(1, 2, 3), volmdlr.Vector3D(1, 0, 0), volmdlr.Vector3D(0, 1, 0)
        )

    def test_parametric_points_to_3d(self):
        parametric_points = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0],
                                      [-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]])
        points3d = self.plane2.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[0.292893218813, 1.0, 0.292893218813], [0.292893218813, 3.0, 0.292893218813],
                                    [1.707106781187, 3.0, 1.707106781187], [1.707106781187, 1.0, 1.707106781187],
                                    [0.646446609407, 1.5, 0.646446609407], [0.646446609407, 2.5, 0.646446609407],
                                    [1.353553390593, 2.5, 1.353553390593], [1.353553390593, 1.5, 1.353553390593]])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0)

    def test_from_normal(self):
        plane = Plane3D.from_normal(self.point1, self.vector3)
        self.assertEqual(plane.frame, volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D, -volmdlr.Y3D, volmdlr.Z3D))

    def test_from_plane_vectors(self):
        plane = Plane3D.from_plane_vectors(self.point1, self.vector1, self.vector2)
        self.assertEqual(plane.frame, volmdlr.OXYZ)

    def test_from_points(self):
        # Test with only two points
        points = [self.point1, self.point2]
        with self.assertRaises(ValueError):
            Plane3D.from_points(points)

        # Test with three points
        plane = Plane3D.from_points([self.point1, self.point2, self.point3])
        self.assertEqual(plane.frame, volmdlr.OXYZ)

        # Test with more than three points
        points = [self.point1, self.point2, self.point3, self.point4]
        plane = Plane3D.from_points(points)
        self.assertEqual(plane.frame, volmdlr.OXYZ)

    def test_angle_between_planes(self):
        # Test with two orthogonal planes
        plane1 = Plane3D.from_normal(self.point1, self.vector1)
        plane2 = Plane3D.from_normal(self.point1, self.vector2)
        angle = plane1.angle_between_planes(plane2)
        self.assertAlmostEqual(angle, math.pi / 2)

        # Test with two parallel planes
        plane1 = Plane3D.from_normal(self.point1, self.vector1)
        plane2 = Plane3D.from_normal(self.point2, self.vector1)
        angle = plane1.angle_between_planes(plane2)
        self.assertAlmostEqual(angle, 0)

    def test_point_belongs(self):
        # Test with point on the plane
        plane = Plane3D.from_3_points(self.point1, self.point2, self.point3)
        point = self.point1
        self.assertTrue(plane.point_belongs(point))

        # Test with point above the plane
        point = volmdlr.Point3D(0, 0, 1)
        self.assertFalse(plane.point_belongs(point))

        # Test with point below the plane
        point = volmdlr.Point3D(0, 0, -1)
        self.assertFalse(plane.point_belongs(point))

    def test_point_distance(self):
        # test point above the plane
        point = volmdlr.Point3D(1, 2, 4)
        self.assertAlmostEqual(self.plane5.point_distance(point), 1.0)

        # test point below the plane
        point = volmdlr.Point3D(1, 2, 2)
        self.assertAlmostEqual(self.plane5.point_distance(point), 1.0)

        # test point on the plane
        point = volmdlr.Point3D(1, 2, 3)
        self.assertAlmostEqual(self.plane5.point_distance(point), 0.0)

    def test_is_parallel(self):
        self.assertTrue(self.plane1.is_parallel(self.plane3))
        self.assertFalse(self.plane1.is_parallel(self.plane2))

    def test_is_coincident(self):
        self.assertTrue(self.plane1.is_coincident(self.plane4))
        self.assertFalse(self.plane1.is_coincident(self.plane2))

    def test_plane_intersections(self):
        plane_intersections = self.plane1.plane_intersections(self.plane2)
        self.assertEqual(len(plane_intersections), 1)
        self.assertTrue(plane_intersections[0].is_close(
            curves.Line3D(volmdlr.O3D, volmdlr.Point3D(0, 0.7071067811865476, 0))))
        no_plane_intersections = self.plane1.plane_intersections(self.plane3)
        self.assertFalse(no_plane_intersections)
        plane1 = surfaces.Plane3D(
            volmdlr.Frame3D(
                volmdlr.Point3D(2.47172762684, 0.709056119825, 0.533657243895),
                volmdlr.Vector3D(0.08730196938518492, 0.9961818941044193, 0.0),
                volmdlr.Vector3D(-0.36467438001762453, 0.031958813694844, -0.9305864982826579),
                volmdlr.Vector3D(-0.9270334204872172, 0.08124203398333905, 0.3660720819920861),
            )
        )
        plane2 = surfaces.Plane3D(
            volmdlr.Frame3D(
                volmdlr.Point3D(2.535691031746372, 0.7426189496471666, 0.6712946669810791),
                volmdlr.Vector3D(0.08730196938518722, 0.9961818941044189, 0.0),
                volmdlr.Vector3D(0.9270334204872168, -0.08124203398334119, -0.36607208199208596),
                volmdlr.Vector3D(-0.3646743800176244, 0.03195881369484484, -0.9305864982826579),
            )
        )
        expected_line = curves.Line3D(
            volmdlr.Point3D(2.4648333822539743, 0.0, 0.6735585604963772),
            volmdlr.Point3D(2.377531412868789, -0.9961818941044192, 0.6735585604963764),
        )
        plane_intersections2 = plane1.plane_intersections(plane2)
        self.assertEqual(expected_line, plane_intersections2[0])
        plane1, plane2 = DessiaObject.from_json(
            os.path.join(folder, "test_plane_plane_intersections301123.json")
        ).primitives
        plane_intersections = plane1.plane_intersections(plane2)
        self.assertTrue(plane_intersections[0].point1.is_close(volmdlr.Point3D(0.0, 0.303510045726, 0.0)))
        self.assertTrue(plane_intersections[0].point2.is_close(
            volmdlr.Point3D(1.0000000000000002, 0.3035100457150452, 0.0)))

    def test_line_intersections(self):
        # test line intersects the plane
        line = curves.Line3D(volmdlr.Point3D(1, 2, 1), volmdlr.Point3D(1, 2, 5))
        expected_intersection = [volmdlr.Point3D(1.0, 2.0, 3.0)]
        self.assertEqual(self.plane5.line_intersections(line), expected_intersection)

        # test line is parallel to the plane
        line = curves.Line3D(volmdlr.Point3D(1, 1, 3), volmdlr.Point3D(2, 2, 3))
        expected_intersection = []
        self.assertEqual(self.plane5.line_intersections(line), expected_intersection)

    def test_linesegment_intersections(self):
        # test linesegment intersects the plane
        linesegment = edges.LineSegment3D(volmdlr.Point3D(1, 2, 1), volmdlr.Point3D(1, 2, 5))
        expected_intersection = [volmdlr.Point3D(1.0, 2.0, 3.0)]
        self.assertEqual(self.plane5.linesegment_intersections(linesegment), expected_intersection)

        # test linesegment does not intersect the plane
        linesegment = edges.LineSegment3D(volmdlr.Point3D(1, 1, 3), volmdlr.Point3D(2, 2, 3))
        expected_intersection = []
        self.assertEqual(self.plane5.linesegment_intersections(linesegment), expected_intersection)

    def test_fullarc_intersections(self):
        start_end1 = self.plane2.frame.origin + self.plane2.frame.u * 3
        circle1 = curves.Circle3D(self.plane2.frame, self.plane2.frame.origin.point_distance(start_end1))
        fullarc1 = edges.FullArc3D(circle1, start_end1)
        start_end2 = self.plane3.frame.origin + self.plane3.frame.u * 3
        circle2 = curves.Circle3D(self.plane3.frame, self.plane3.frame.origin.point_distance(start_end2))
        fullarc2 = edges.FullArc3D(circle2, start_end2)
        fullarc_intersections = self.plane1.fullarc_intersections(fullarc1)
        self.assertEqual(len(fullarc_intersections), 2)
        self.assertTrue(fullarc_intersections[0].is_close(volmdlr.Point3D(0, 4.645751311, 0)))
        self.assertTrue(fullarc_intersections[1].is_close(volmdlr.Point3D(0, -0.645751311, 0)))
        fullarc_intersections2 = self.plane1.fullarc_intersections(fullarc2)
        self.assertFalse(fullarc_intersections2)

    def test_arc_intersections(self):
        plane = Plane3D.from_plane_vectors(
            volmdlr.Point3D(1, 2, 3), volmdlr.Vector3D(1, 0, 0), volmdlr.Vector3D(0, 1, 0)
        )
        arc = edges.Arc3D.from_3_points(volmdlr.Point3D(2, 2, 4), volmdlr.Point3D(2, 3, 3), volmdlr.Point3D(2, 2, 2))
        arc_intersections = plane.arc_intersections(arc)
        self.assertEqual(len(arc_intersections), 1)
        self.assertTrue(arc_intersections[0].is_close(volmdlr.Point3D(2.0, 3.0, 3.0)))
        arc2 = edges.Arc3D.from_3_points(
            volmdlr.Point3D(2, 1, 2.5), volmdlr.Point3D(2, 2, 3.5), volmdlr.Point3D(2, 3, 2.5)
        )
        arc_intersections2 = plane.arc_intersections(arc2)
        self.assertEqual(len(arc_intersections2), 2)
        self.assertTrue(arc_intersections2[0].is_close(volmdlr.Point3D(2.0, 1.1339745962155614, 3.0)))
        self.assertTrue(arc_intersections2[1].is_close(volmdlr.Point3D(2.0, 2.8660254037844384, 3.0)))

    def test_bspline_intersections(self):
        plane = Plane3D(volmdlr.OZXY)
        intersections = plane.bsplinecurve_intersections(bspline_curve=bspline_curve3d())
        self.assertTrue(len(intersections), 2)
        self.assertTrue(intersections[0].is_close(volmdlr.Point3D(3.4032921805712286, 0.0, 1.5267489712255913)))
        self.assertTrue(intersections[1].is_close(volmdlr.Point3D(-0.49607254405814866, 0.0, -0.6842127097914265)))

    def test_from_3_points(self):
        """
        Test we can create from 3 points, that should be on the surface & a that the frame.w is orthogonal.
        """
        for p1, p2, p3 in [
            (
                volmdlr.Point3D(-0.18141727653547446, -0.3129888988150751, -0.7869168525937733),
                volmdlr.Point3D(0.1756754261634219, -0.22088840984091496, -0.2482699866635696),
                volmdlr.Point3D(0.7373327193311012, -0.8849251090904118, -0.9464563060452031),
            ),
            (
                volmdlr.Point3D(9.80579460566209, 6.917108655469294, 0.0),
                volmdlr.Point3D(12.0, 0.0, 0.6956521739130439),
                volmdlr.Point3D(12.0, 0.0, 0.5217391304347831),
            ),
            (
                volmdlr.Point3D(8.773150355532506, 8.18729704110092, 0.0),
                volmdlr.Point3D(12.0, 0.0, 0.8695652173913047),
                volmdlr.Point3D(12.0, 0.0, 0.6956521739130439),
            ),
        ]:

            surface = surfaces.Plane3D.from_3_points(p1, p2, p3)
            self.assertTrue(surface.point_belongs(p1))
            self.assertTrue(surface.point_belongs(p2))
            self.assertTrue(surface.point_belongs(p3))
            self.assertAlmostEqual(surface.frame.w.dot(p1 - p2), 0.0)
            self.assertAlmostEqual(surface.frame.w.dot(p3 - p2), 0.0)
            self.assertAlmostEqual(surface.frame.w.dot(p3 - p1), 0.0)

    def test_plane_between_two_planes(self):
        plane1 = Plane3D(volmdlr.OXYZ)
        plane2 = Plane3D(volmdlr.OYZX)
        plane_b = Plane3D.plane_between_two_planes(plane1, plane2)
        expected_frame = volmdlr.Frame3D(
            volmdlr.O3D,
            volmdlr.Vector3D(0, 1, 0),
            volmdlr.Vector3D(0.7071067811865475, 0.0, -0.7071067811865475),
            volmdlr.Vector3D(0.7071067811865475, 0.0, 0.7071067811865475),
        )
        self.assertEqual(plane_b.frame, expected_frame)

    def test_rotation(self):
        plane1 = Plane3D(volmdlr.OXYZ)
        rotated_plane1 = plane1.rotation(volmdlr.O3D, volmdlr.Vector3D(0, 1, 0), math.pi / 4)
        expected_frame = volmdlr.Frame3D(
            volmdlr.O3D,
            volmdlr.Vector3D(0.7071067811865476, 0.0, -0.7071067811865475),
            volmdlr.Vector3D(0, 1, 0),
            volmdlr.Vector3D(0.7071067811865475, 0.0, 0.7071067811865476),
        )
        self.assertEqual(rotated_plane1.frame, expected_frame)

    def test_contour3d_to_2d(self):
        plane = surfaces.Plane3D.from_json(os.path.join(folder, "plane_parametric_operation_bug_surface.json"))
        contour3d = wires.Contour3D.from_json(os.path.join(folder, "plane_parametric_operation_bug_contour.json"))
        contour = plane.contour3d_to_2d(contour3d)
        self.assertTrue(contour.is_ordered())
        self.assertAlmostEqual(contour.area(), 8.120300532917004e-06)

        plane = surfaces.Plane3D.from_json(os.path.join(folder, "planesurface_arc3d_to_2d.json"))
        contour3d = wires.Contour3D.from_json(os.path.join(folder, "planesurface_arc3d_to_2d_contour.json"))
        contour = plane.contour3d_to_2d(contour3d)
        self.assertTrue(contour.is_ordered())
        self.assertAlmostEqual(contour.area(), math.pi * 0.202 ** 2, 4)


if __name__ == "__main__":
    unittest.main()
