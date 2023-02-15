import math
import unittest

import volmdlr
from volmdlr import edges, faces


class TestPlane3D(unittest.TestCase):
    plane1 = faces.Plane3D(volmdlr.OXYZ)
    plane2 = faces.Plane3D(
        volmdlr.OYZX.translation(volmdlr.Vector3D(1, 2, 1)).rotation(volmdlr.Point3D(1, 2, 1), volmdlr.Y3D, math.pi / 4)
    )
    plane3 = faces.Plane3D(
        volmdlr.OXYZ.translation(volmdlr.Vector3D(0, 0, 1)).rotation(
            volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4
        )
    )
    plane4 = faces.Plane3D(volmdlr.OXYZ.rotation(volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4))

    def test_plane_intersections(self):
        plane_intersections = self.plane1.plane_intersection(self.plane2)
        self.assertEqual(len(plane_intersections), 1)
        self.assertEqual(plane_intersections[0], edges.Line3D(volmdlr.O3D, volmdlr.Point3D(0, 0.707106781, 0)))
        no_plane_intersections = self.plane1.plane_intersection(self.plane3)
        self.assertFalse(no_plane_intersections)

    def test_is_parallel(self):
        self.assertTrue(self.plane1.is_parallel(self.plane3))
        self.assertFalse(self.plane1.is_parallel(self.plane2))

    def test_is_coincident(self):
        self.assertTrue(self.plane1.is_coincident(self.plane4))
        self.assertFalse(self.plane1.is_coincident(self.plane2))

    def test_fullarc_intersections(self):
        fullarc1 = edges.FullArc3D(
            self.plane2.frame.origin, self.plane2.frame.origin + self.plane2.frame.u * 3, self.plane2.frame.w
        )

        fullarc2 = edges.FullArc3D(
            self.plane3.frame.origin, self.plane3.frame.origin + self.plane3.frame.u * 3, self.plane3.frame.w
        )
        fullarc_intersections = self.plane1.fullarc_intersections(fullarc1)
        self.assertEqual(len(fullarc_intersections), 2)
        self.assertEqual(fullarc_intersections[0], volmdlr.Point3D(0, 4.645751311, 0))
        self.assertEqual(fullarc_intersections[1], volmdlr.Point3D(0, -0.645751311, 0))
        fullarc_intersections2 = self.plane1.fullarc_intersections(fullarc2)
        self.assertFalse(fullarc_intersections2)


if __name__ == "__main__":
    unittest.main()
