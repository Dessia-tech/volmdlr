import unittest

import volmdlr
from volmdlr import edges, wires


class TestCircle3D(unittest.TestCase):

    def test_frame_mapping(self):
        for side in ['old', 'new']:
            circle = wires.Circle3D(volmdlr.Frame3D(volmdlr.Point3D(1.456, -0.12456, -0.457),
                                                    volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D), 0.32)
            frame = volmdlr.Frame3D(volmdlr.Point3D(0.12, 0.1, -0.07), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
            framed_mapped_circle = circle.frame_mapping(frame, side)
            framed_center = circle.frame.origin.frame_mapping(frame, side)
            framed_normal = circle.frame.w.frame_mapping(frame, side)
            self.assertAlmostEqual(framed_center.point_distance(framed_mapped_circle.frame.origin), 0.)
            self.assertAlmostEqual(framed_normal.point_distance(framed_mapped_circle.frame.w), 0.)
            self.assertEqual(circle.radius, framed_mapped_circle.radius)

    def test_linsegment_intersections(self):
        circle = wires.Circle3D(volmdlr.OXYZ, 1)
        lineseg = edges.LineSegment3D(volmdlr.Point3D(0, 0, -1), volmdlr.Point3D(2, 0, 1))
        line_seg2 = edges.LineSegment3D(volmdlr.Point3D(1, 2, 0), volmdlr.Point3D(-1, -2, 0))
        circle_linseg_intersections = circle.linesegment_intersections(lineseg) + circle.linesegment_intersections(
            line_seg2)
        self.assertEqual(len(circle_linseg_intersections), 3)
        expected_intersections = [volmdlr.Point3D(1.0, 0.0, 0.0),
                                  volmdlr.Point3D(0.447213595499958, 0.894427190999916, 0.0),
                                  volmdlr.Point3D(-0.447213595499958, -0.8944271909999157, 0.0)]
        for expected_point, point in zip(expected_intersections, circle_linseg_intersections):
            self.assertTrue(expected_point.is_close(point))
        circle2 = wires.Circle3D(volmdlr.OYZX, 1)
        lineseg2 = edges.LineSegment3D(volmdlr.Point3D(-1, 0, -2), volmdlr.Point3D(2, 0, 1))
        circle_linseg_intersections = circle2.linesegment_intersections(lineseg2)
        self.assertEqual(len(circle_linseg_intersections), 1)
        self.assertTrue(circle_linseg_intersections[0].is_close(volmdlr.Point3D(0.0, 0.0, -1.0)))
        circle3 = wires.Circle3D(
            volmdlr.Frame3D(volmdlr.O3D, volmdlr.Vector3D(0.2672612419124244, 0.5345224838248488, 0.8017837257372732),
                            volmdlr.Vector3D(0.9636241116594316, -0.14824986333222026, -0.22237479499833038),
                            volmdlr.Vector3D(0, 0.8320502943378438, -0.5547001962252291)), 1)
        circle_linseg_intersections1 = circle3.linesegment_intersections(lineseg2)
        self.assertEqual(len(circle_linseg_intersections1), 1)
        self.assertTrue(circle_linseg_intersections1[0].is_close(volmdlr.Point3D(1, 0.0, 0.0)))


if __name__ == '__main__':
    unittest.main()
