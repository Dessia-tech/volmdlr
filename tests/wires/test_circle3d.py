import unittest

import volmdlr
from volmdlr import edges, wires


class TestCircle3D(unittest.TestCase):
    def test_linsegment_intersections(self):
        circle = wires.Circle3D(volmdlr.OXYZ, 1)
        lineseg = edges.LineSegment3D(volmdlr.Point3D(0, 0, -1), volmdlr.Point3D(2, 0, 1))
        line_seg2 = edges.LineSegment3D(volmdlr.Point3D(1, 2, 0), volmdlr.Point3D(-1, -2, 0))
        circle_linseg_intersections = circle.linesegment_intersections(lineseg) + circle.linesegment_intersections(
            line_seg2
        )
        self.assertEqual(len(circle_linseg_intersections), 3)
        expected_intersections = [
            volmdlr.Point3D(1.0, 0.0, 0.0),
            volmdlr.Point3D(0.447213595499958, 0.894427190999916, 0.0),
            volmdlr.Point3D(-0.447213595499958, -0.8944271909999157, 0.0),
        ]
        for expected_point, point in zip(expected_intersections, circle_linseg_intersections):
            self.assertEqual(expected_point, point)
        circle2 = wires.Circle3D(volmdlr.OYZX, 1)
        lineseg2 = edges.LineSegment3D(volmdlr.Point3D(-1, 0, -2), volmdlr.Point3D(2, 0, 1))
        circle_linseg_intersections = circle2.linesegment_intersections(lineseg2)
        self.assertEqual(len(circle_linseg_intersections), 1)
        self.assertEqual(circle_linseg_intersections[0], volmdlr.Point3D(0.0, 0.0, -1.0))
        circle3 = wires.Circle3D(
            volmdlr.Frame3D(
                volmdlr.O3D,
                volmdlr.Vector3D(0.2672612419124244, 0.5345224838248488, 0.8017837257372732),
                volmdlr.Vector3D(0.9636241116594316, -0.14824986333222026, -0.22237479499833038),
                volmdlr.Vector3D(0, 0.8320502943378438, -0.5547001962252291),
            ),
            1,
        )
        circle_linseg_intersections1 = circle3.linesegment_intersections(lineseg2)
        self.assertEqual(len(circle_linseg_intersections1), 1)
        self.assertEqual(circle_linseg_intersections1[0], volmdlr.Point3D(1, 0.0, 0.0))


if __name__ == "__main__":
    unittest.main()
