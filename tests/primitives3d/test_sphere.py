import unittest
from math import pi

import volmdlr
from volmdlr.primitives3d import Sphere


class TestCylinder(unittest.TestCase):
    def setUp(self):
        self.sphere = Sphere(
            center=volmdlr.O3D,
            radius=0.02,
        )

    def test_point_belongs(self):
        self.assertTrue(self.sphere.point_belongs(volmdlr.O3D))
        self.assertFalse(self.sphere.point_belongs(volmdlr.Point3D(1, 0, 0)))

    def test_volume(self):
        volume = 4 / 3 * pi * self.sphere.radius**3
        self.assertEqual(self.sphere.volume(), volume)


if __name__ == "__main__":
    unittest.main()
