import math
import unittest

import volmdlr
from volmdlr import edges, wires


class TestArc3D(unittest.TestCase):
    arc3d = edges.Arc3D(volmdlr.Point3D(-1, 0, 0), volmdlr.Point3D(0, -1, 0), volmdlr.Point3D(1, 0, 0))
    arc3d_2 = edges.Arc3D(volmdlr.Point3D(-0.2672612419124244, -0.5345224838248488, -0.8017837257372732),
                            volmdlr.Point3D(0.9636241116594316, -0.14824986333222026, -0.22237479499833038),
                            volmdlr.Point3D(0.2672612419124244, 0.5345224838248488, 0.8017837257372732))

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
        self.assertTrue(inters3[0].is_close(expected_point))
        lineseg = edges.LineSegment3D(volmdlr.Point3D(-1, 0, -2), volmdlr.Point3D(2, 0, 1))
        arc3d_lineseg_inters = self.arc3d_2.linesegment_intersections(lineseg)
        self.assertEqual(arc3d_lineseg_inters[0], volmdlr.Point3D(1.0, 0.0, 0.0))
        lineseg_ = edges.LineSegment3D(volmdlr.Point3D(-1, 0, -2), volmdlr.Point3D(-1, 0, -2) +
                                       lineseg.length() * 0.6 * lineseg.unit_direction_vector())
        no_arc3d_lineseg_inters = self.arc3d_2.linesegment_intersections(lineseg_)
        self.assertFalse(no_arc3d_lineseg_inters)

    def test_line_intersections(self):
        lineseg = edges.LineSegment3D(volmdlr.Point3D(-1, 0, -2), volmdlr.Point3D(2, 0, 1))
        line = edges.Line3D(volmdlr.Point3D(-1, 0, -2), volmdlr.Point3D(-1, 0, -2) +
                            lineseg.length() * 0.6 * lineseg.unit_direction_vector())
        arc3d_line_inters = self.arc3d_2.line_intersections(line)
        self.assertEqual(arc3d_line_inters[0], volmdlr.Point3D(1.0, 0.0, 0.0))

    def test_split(self):
        split1 = self.arc3d.split(self.arc3d.start)
        self.assertIsNone(split1[0])
        self.assertEqual(split1[1], self.arc3d)
        split2 = self.arc3d.split(self.arc3d.end)
        self.assertEqual(split2[0], self.arc3d)
        self.assertIsNone(split2[1])
        split3 = self.arc3d.split(self.arc3d.interior)
        self.assertTrue(split3[0].start.is_close(self.arc3d.start))
        self.assertTrue(split3[0].end.is_close(self.arc3d.interior))
        self.assertTrue(split3[1].start.is_close(self.arc3d.interior))
        self.assertTrue(split3[1].end.is_close(self.arc3d.end))


if __name__ == '__main__':
    unittest.main()
