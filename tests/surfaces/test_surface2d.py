import unittest

import volmdlr

from volmdlr import surfaces, wires


class TestSurface2D(unittest.TestCase):
    contour = wires.Contour2D.from_points([volmdlr.O2D, volmdlr.Point2D(1, 0), volmdlr.Point2D(1, 1),
                                           volmdlr.Point2D(0, 1)])
    surface2d = surfaces.Surface2D(contour, [])

    def test_triangulation(self):
        tri = self.surface2d.triangulation()
        self.assertAlmostEqual(tri.area(), 1)

        surface2d = surfaces.Surface2D.load_from_file("surfaces/objects_surface2d_test/self_intersections.json")
        tri = surface2d.triangulation()
        self.assertTrue(tri)


if __name__ == '__main__':
    unittest.main()
