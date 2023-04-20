import math
import unittest

import volmdlr
from volmdlr import edges, faces


class TestCylindricalFace3D(unittest.TestCase):
    surface = faces.CylindricalSurface3D(volmdlr.OXYZ, 0.32)
    cylindrical_face = surface.rectangular_cut(-0.01, 1.3, -0.1, 0.3)

    def test_linesegment_intersections(self):
        lineseg3d = edges.LineSegment3D(volmdlr.O3D, volmdlr.Point3D(0.3, 0.3, .3))
        line_inters = self.cylindrical_face.linesegment_intersections(lineseg3d)
        self.assertEqual(len(line_inters), 1)
        self.assertTrue(line_inters[0].is_close(volmdlr.Point3D(0.22627416, 0.22627416, 0.22627416)))
        cylindrical_face1 = self.surface.rectangular_cut(-math.pi/4, 1.5 * math.pi, -0.1, 0.3)
        lineseg3d_2 = edges.LineSegment3D(volmdlr.Point3D(-0.3, -0.3, -.1), volmdlr.Point3D(0.3, 0.3, .3))
        line_inters_2 = cylindrical_face1.linesegment_intersections(lineseg3d_2)
        self.assertEqual(len(line_inters_2), 2)
        self.assertTrue(line_inters_2[0].is_close(volmdlr.Point3D(0.2262741, 0.2262741, 0.2508494)))
        self.assertTrue(line_inters_2[1].is_close(volmdlr.Point3D(-0.2262741, -0.2262741, -0.0508494)))


if __name__ == '__main__':
    unittest.main()
