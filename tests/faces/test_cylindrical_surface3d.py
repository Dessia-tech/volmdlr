import unittest

import volmdlr
from volmdlr import faces, edges


class TestCylindricalSurface3D(unittest.TestCase):
    cylindrical_surface = faces.CylindricalSurface3D(volmdlr.OXYZ, 0.32)

    def test_linesegment_intersections(self):
        line3d = edges.Line3D(volmdlr.O3D, volmdlr.Point3D(0.3, 0.3, .3))
        line_inters = self.cylindrical_surface.line_intersections(line3d)
        self.assertEqual(len(line_inters), 2)
        self.assertEqual(line_inters[0], volmdlr.Point3D(0.22627416997969518, 0.22627416997969518,
                                                         0.22627416997969518))
        self.assertEqual(line_inters[1], volmdlr.Point3D(-0.22627416997969518, -0.22627416997969518,
                                                         -0.22627416997969518))


if __name__ == '__main__':
    unittest.main()
