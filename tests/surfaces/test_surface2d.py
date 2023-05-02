import unittest

import volmdlr

from volmdlr import surfaces, wires


class TestSurface2D(unittest.TestCase):
    contour = wires.Contour2D.from_points([volmdlr.O2D, volmdlr.Point2D(1, 0), volmdlr.Point2D(1, 1),
                                           volmdlr.Point2D(0, 1)])
    surface2d = surfaces.Surface2D(contour, [])

    def test_triangulation(self):
        tri = self.surface2d.triangulation()
        tri.plot()
        self.assertAlmostEqual(self.surface2d.triangulation().area(), 1)


if __name__ == '__main__':
    unittest.main()
