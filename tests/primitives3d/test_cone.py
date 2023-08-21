import unittest

import volmdlr
from volmdlr.primitives3d import Cone


class TestCone(unittest.TestCase):
    def setUp(self) -> None:
        self.radius = 0.1
        self.length = 0.1
        self.cone = Cone(frame=volmdlr.OXYZ, radius=self.radius, length=self.length)

    def test_point_belongs(self):
        self.assertTrue(self.cone.point_belongs(volmdlr.O3D))
        self.assertTrue(self.cone.point_belongs(volmdlr.Point3D(0.0, 0.0, self.length / 2)))
        self.assertTrue(self.cone.point_belongs(volmdlr.Point3D(-self.radius, 0.0, -self.length / 2)))
        self.assertTrue(self.cone.point_belongs(volmdlr.Point3D(0.0, -self.radius, -self.length / 2)))
        self.assertFalse(self.cone.point_belongs(volmdlr.Point3D(0.0, 0.0, self.length)))
        self.assertFalse(self.cone.point_belongs(volmdlr.Point3D(0.01, 0.01, self.length / 2)))
