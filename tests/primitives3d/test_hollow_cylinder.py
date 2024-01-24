import math
import unittest

import volmdlr
from volmdlr.primitives3d import HollowCylinder


class TestHollowCylinder(unittest.TestCase):
    def setUp(self) -> None:
        self.inner_radius = 0.1
        self.outer_radius = 0.2
        self.length = 0.3
        self.hollow_cylinder = HollowCylinder(
            frame=volmdlr.OXYZ,
            inner_radius=self.inner_radius,
            outer_radius=self.outer_radius,
            length=self.length,
        )

    def test_point_belongs(self):
        self.assertFalse(self.hollow_cylinder.point_inside(volmdlr.O3D))
        self.assertTrue(
            self.hollow_cylinder.point_inside(volmdlr.Point3D((self.inner_radius + self.outer_radius) / 2, 0, 0))
        )

    def test_from_end_points(self):
        point1 = volmdlr.Point3D(-1, 0, 0)
        point2 = volmdlr.Point3D(1, 0, 0)

        cylinder = HollowCylinder.from_end_points(point1, point2, self.inner_radius, self.outer_radius)

        self.assertEqual(cylinder.position, volmdlr.O3D)
        self.assertEqual(cylinder.axis, volmdlr.X3D)
        self.assertEqual(cylinder.length, 2.0)
        self.assertEqual(cylinder.inner_radius, self.inner_radius)
        self.assertEqual(cylinder.outer_radius, self.outer_radius)

    def test_from_center_point_and_axis(self):
        cylinder = HollowCylinder.from_center_point_and_axis(
            volmdlr.O3D, volmdlr.Z3D, self.inner_radius, self.outer_radius, self.length
        )

        self.assertEqual(cylinder.position, volmdlr.O3D)
        self.assertEqual(cylinder.axis, volmdlr.Z3D)
        self.assertEqual(cylinder.inner_radius, self.inner_radius)
        self.assertEqual(cylinder.outer_radius, self.outer_radius)

    def test_rotation(self):
        rotated_hollow_cylinder = self.hollow_cylinder.rotation(
            center=volmdlr.O3D, axis=volmdlr.X3D, angle=-math.pi / 2
        )
        self.assertTrue(rotated_hollow_cylinder.axis.is_close(volmdlr.Y3D))

    def test_translation(self):
        translated_hollow_cylinder = self.hollow_cylinder.translation(offset=volmdlr.X3D)
        self.assertTrue(translated_hollow_cylinder.position.is_close(volmdlr.Point3D(1.0, 0.0, 0.0)))


if __name__ == "__main__":
    unittest.main()
