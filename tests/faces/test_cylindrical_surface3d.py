import unittest

import volmdlr
from volmdlr import faces, edges


class TestCylindricalSurface3D(unittest.TestCase):
    cylindrical_surface = faces.CylindricalSurface3D(volmdlr.OXYZ, 0.32)

    def test_linesegment_intersections(self):
        line3d = edges.LineSegment3D(volmdlr.O3D, volmdlr.Point3D(0.3, 0.3, .3))
        lineseg_inters = self.cylindrical_surface.linesegment_intersections(line3d)
        self.assertEqual(len(lineseg_inters), 2)
        self.assertEqual(lineseg_inters[0], volmdlr.Point3D(0.22627416997969518, 0.22627416997969518,
                                                            0.22627416997969518))
        self.assertEqual(lineseg_inters[1], volmdlr.Point3D(-0.22627416997969518, -0.22627416997969518,
                                                            -0.22627416997969518))


if __name__ == '__main__':
    unittest.main()
