"""
Tests for faces 3D
"""
import unittest
import volmdlr
from volmdlr import faces, surfaces


class TestFace3D(unittest.TestCase):
    def test_point_distance(self):
        radius = 0.15
        cylindricalsurface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, radius)
        cylindricalface = faces.CylindricalFace3D.from_surface_rectangular_cut(cylindricalsurface, 0,
                                                                               volmdlr.TWO_PI / 3, -.25,
                                                                               .25)
        point3d = volmdlr.Point3D(.05, .05, -0.05)
        distance, point1 = cylindricalface.point_distance(point3d, True)
        self.assertAlmostEqual(distance, 0.07884386917468453)
        self.assertTrue(point1.is_close(volmdlr.Point3D(0.101631724322, 0.109586186278, -0.049915413859)))


if __name__ == '__main__':
    unittest.main()
