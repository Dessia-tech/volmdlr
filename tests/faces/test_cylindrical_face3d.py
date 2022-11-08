import unittest
import volmdlr
from volmdlr import faces, edges


class TestCylindricalFace3D(unittest.TestCase):
    surface = faces.CylindricalSurface3D(volmdlr.OXYZ, 0.32)
    cylindrical_face = surface.rectangular_cut(-0.01, 1.3, -0.1, 0.3)

    def test_linesegment_intersections(self):
        lineseg3d = edges.LineSegment3D(volmdlr.O3D, volmdlr.Point3D(0.3, 0.3, .3))
        line_inters = self.cylindrical_face.linesegment_intersections(lineseg3d)
        self.assertEqual(len(line_inters), 1)
        self.assertEqual(line_inters[0], volmdlr.Point3D(0.22627416997969518, 0.22627416997969518, 0.22627416997969518))


if __name__ == '__main__':
    unittest.main()
