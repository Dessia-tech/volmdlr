import math
import unittest
import volmdlr
from volmdlr import wires, edges


class TestArc3D(unittest.TestCase):
    arc3d = edges.Arc3D(volmdlr.Point3D(-1, 0, 0), volmdlr.Point3D(0, -1, 0), volmdlr.Point3D(1, 0, 0))

    def test_point_belongs(self):
        point1 = volmdlr.Point3D(0, 1, 0)
        point2 = volmdlr.Point3D(-1, -1, 0)
        point3 = volmdlr.Point3D(-math.sqrt(2) / 2, -math.sqrt(2) / 2, 0)
        point4 = volmdlr.Point3D(-.5, -.5, 0)
        point5 = volmdlr.Point3D(0, 0, 1)
        self.assertFalse(self.arc3d.point_belongs(point1))
        self.assertFalse(self.arc3d.point_belongs(point2))
        self.assertTrue(self.arc3d.point_belongs(point3))
        self.assertFalse(self.arc3d.point_belongs(point4))
        self.assertFalse(self.arc3d.point_belongs(point5))

    def test_linesegment_intersections(self):
        lineseg1 = edges.LineSegment3D(volmdlr.Point3D(2, 2, 0), volmdlr.Point3D(-2, -2, 0))
        lineseg2 = edges.LineSegment3D(volmdlr.Point3D(0, -1, -1), volmdlr.Point3D(0, 1, 0))
        lineseg3 = edges.LineSegment3D(volmdlr.Point3D(-1, -1, -1),
                                       volmdlr.Point3D(-math.sqrt(2) / 2, -math.sqrt(2) / 2, 0))
        inters1 = self.arc3d.linesegment_intersections(lineseg1)
        inters2 = self.arc3d.linesegment_intersections(lineseg2)
        inters3 = self.arc3d.linesegment_intersections(lineseg3)
        expected_point = volmdlr.Point3D(-0.7071067811865476, -0.7071067811865475, 0.0)
        self.assertEqual(len(inters1), 1)
        self.assertEqual(inters1[0], expected_point)
        self.assertFalse(inters2)
        self.assertEqual(len(inters3), 1)
        self.assertEqual(inters3[0], expected_point)


if __name__ == '__main__':
    unittest.main()
