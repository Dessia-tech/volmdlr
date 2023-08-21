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
        self.assertFalse(self.hollow_cylinder.point_belongs(volmdlr.O3D))
        self.assertTrue(
            self.hollow_cylinder.point_belongs(volmdlr.Point3D((self.inner_radius + self.outer_radius) / 2, 0, 0))
        )

    def test_from_extremal_points(self):
        point1 = volmdlr.Point3D(-1, 0, 0)
        point2 = volmdlr.Point3D(1, 0, 0)

        cylinder = HollowCylinder.from_extremal_points(point1, point2, self.inner_radius, self.outer_radius)

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


if __name__ == "__main__":
    unittest.main()
