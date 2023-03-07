import math
import unittest

import volmdlr
from volmdlr import edges, faces


class TestPlane3D(unittest.TestCase):
    plane1 = faces.Plane3D(volmdlr.OXYZ)
    plane2 = faces.Plane3D(volmdlr.OYZX.translation(volmdlr.Vector3D(1, 2, 1)).rotation(
        volmdlr.Point3D(1, 2, 1), volmdlr.Y3D, math.pi / 4))
    plane3 = faces.Plane3D(volmdlr.OXYZ.translation(volmdlr.Vector3D(0, 0, 1)).rotation(
        volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4))
    plane4 = faces.Plane3D(volmdlr.OXYZ.rotation(volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4))

    def test_plane_intersections(self):
        plane_intersections = self.plane1.plane_intersection(self.plane2)
        self.assertEqual(len(plane_intersections), 1)
        self.assertEqual(plane_intersections[0], edges.Line3D(volmdlr.O3D, volmdlr.Point3D(0, 0.7071067811865476, 0)))
        no_plane_intersections = self.plane1.plane_intersection(self.plane3)
        self.assertFalse(no_plane_intersections)
        plane1 = faces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(2.47172762684, 0.709056119825, 0.533657243895),
                                               volmdlr.Vector3D(0.08730196938518492, 0.9961818941044193, 0.0),
                                               volmdlr.Vector3D(-0.36467438001762453, 0.031958813694844,
                                                                -0.9305864982826579),
                                               volmdlr.Vector3D(-0.9270334204872172, 0.08124203398333905,
                                                                0.3660720819920861)))
        plane2 = faces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(2.535691031746372, 0.7426189496471666,
                                                               0.6712946669810791),
                                               volmdlr.Vector3D(0.08730196938518722, 0.9961818941044189, 0.0),
                                               volmdlr.Vector3D(0.9270334204872168, -0.08124203398334119,
                                                                -0.36607208199208596),
                                               volmdlr.Vector3D(-0.3646743800176244, 0.03195881369484484,
                                                                -0.9305864982826579)))
        expected_line = edges.Line3D(volmdlr.Point3D(2.4648333822539743, 0.0, 0.6735585604963772),
                                     volmdlr.Point3D(2.377531412868789, -0.9961818941044192, 0.6735585604963764))
        plane_intersections2 = plane1.plane_intersection(plane2)
        self.assertEqual(expected_line, plane_intersections2[0])

    def test_is_parallel(self):
        self.assertTrue(self.plane1.is_parallel(self.plane3))
        self.assertFalse(self.plane1.is_parallel(self.plane2))

    def test_is_coincident(self):
        self.assertTrue(self.plane1.is_coincident(self.plane4))
        self.assertFalse(self.plane1.is_coincident(self.plane2))

    def test_fullarc_intersections(self):
        fullarc1 = edges.FullArc3D(self.plane2.frame.origin, self.plane2.frame.origin +
                                   self.plane2.frame.u * 3, self.plane2.frame.w)

        fullarc2 = edges.FullArc3D(self.plane3.frame.origin, self.plane3.frame.origin +
                                   self.plane3.frame.u * 3, self.plane3.frame.w)
        fullarc_intersections = self.plane1.fullarc_intersections(fullarc1)
        self.assertEqual(len(fullarc_intersections), 2)
        self.assertTrue(fullarc_intersections[0].is_close(volmdlr.Point3D(0, 4.645751311, 0)))
        self.assertTrue(fullarc_intersections[1].is_close(volmdlr.Point3D(0, -0.645751311, 0)))
        fullarc_intersections2 = self.plane1.fullarc_intersections(fullarc2)
        self.assertFalse(fullarc_intersections2)


if __name__ == '__main__':
    unittest.main()
