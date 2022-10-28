import unittest
import volmdlr

from volmdlr.primitives3d import Cylinder


class TestCylinder(unittest.TestCase):
    cylinder1 = Cylinder(
        position=volmdlr.Point3D(0, 0, 0),
        axis=volmdlr.Vector3D(1, 0, 0),
        radius=0.02,
        length=0.1,
        color=(0, 0, 1),
    )

    cylinder2 = Cylinder(
        position=volmdlr.Point3D(0.05, 0, 0),
        axis=volmdlr.Vector3D(0, 1, 0),
        radius=0.005,
        length=0.01,
        color=(1, 0, 1),
    )

    cylinder3 = Cylinder(
        position=volmdlr.Point3D(0.15, 0, 0),
        axis=volmdlr.Vector3D(0, 1, 0),
        radius=0.005,
        length=0.01,
        color=(1, 0, 1),
    )

    def test_point_belongs(self):
        self.assertTrue(self.cylinder1.point_belongs(volmdlr.O3D))
        self.assertFalse(self.cylinder1.point_belongs(volmdlr.Point3D(1, 0, 0)))

    def test_random_point_inside(self):
        point = self.cylinder1.random_point_inside()
        self.assertIsInstance(point, volmdlr.Point3D)
        self.assertTrue(self.cylinder1.point_belongs(point))

    def test_interference_volume_with_other_cylinder(self):
        interference_volume = self.cylinder1.interference_volume_with_other_cylinder(
            other_cylinder=self.cylinder1
        )
        self.assertEqual(interference_volume, self.cylinder1.volume())

    def test_min_distance_to_other_cylinder(self):
        dist_1_2 = self.cylinder1.min_distance_to_other_cylinder(self.cylinder2)
        dist_1_3 = self.cylinder1.min_distance_to_other_cylinder(self.cylinder3)

        self.assertAlmostEqual(dist_1_2, 0.0, delta=1e-5)
        self.assertAlmostEqual(dist_1_3, 0.095, delta=1e-5)

    def test_is_intersecting_other_cylinder(self):
        intersecting_1_2 = self.cylinder1.is_intersecting_other_cylinder(self.cylinder2)
        intersecting_1_3 = self.cylinder1.is_intersecting_other_cylinder(self.cylinder3)

        self.assertTrue(intersecting_1_2)
        self.assertFalse(intersecting_1_3)


if __name__ == "__main__":
    unittest.main()
