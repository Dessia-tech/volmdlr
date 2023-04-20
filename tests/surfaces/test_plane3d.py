import math
import unittest

import volmdlr
from volmdlr import edges, surfaces


class TestPlane3D(unittest.TestCase):
    plane1 = surfaces.Plane3D(volmdlr.OXYZ)
    plane2 = surfaces.Plane3D(volmdlr.OYZX.translation(volmdlr.Vector3D(1, 2, 1)).rotation(
        volmdlr.Point3D(1, 2, 1), volmdlr.Y3D, math.pi / 4))
    plane3 = surfaces.Plane3D(volmdlr.OXYZ.translation(volmdlr.Vector3D(0, 0, 1)).rotation(
        volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4))
    plane4 = surfaces.Plane3D(volmdlr.OXYZ.rotation(volmdlr.O3D, volmdlr.Vector3D(0, 0, 1), math.pi / 4))

    def test_plane_intersections(self):
        plane_intersections = self.plane1.plane_intersection(self.plane2)
        self.assertEqual(len(plane_intersections), 1)
        self.assertEqual(plane_intersections[0], edges.Line3D(volmdlr.O3D, volmdlr.Point3D(0, 0.7071067811865476, 0)))
        no_plane_intersections = self.plane1.plane_intersection(self.plane3)
        self.assertFalse(no_plane_intersections)
        plane1 = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(2.47172762684, 0.709056119825, 0.533657243895),
                                               volmdlr.Vector3D(0.08730196938518492, 0.9961818941044193, 0.0),
                                               volmdlr.Vector3D(-0.36467438001762453, 0.031958813694844,
                                                                -0.9305864982826579),
                                               volmdlr.Vector3D(-0.9270334204872172, 0.08124203398333905,
                                                                0.3660720819920861)))
        plane2 = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(2.535691031746372, 0.7426189496471666,
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

    def test_from_3_points(self):
        """
        Test we can create from 3 points, that should be on the surface & a that the frame.w is orthogonal.
        """
        for p1, p2, p3 in [(volmdlr.Point3D(-0.18141727653547446, -0.3129888988150751, -0.7869168525937733),
                            volmdlr.Point3D(0.1756754261634219, -0.22088840984091496, -0.2482699866635696),
                            volmdlr.Point3D(0.7373327193311012, -0.8849251090904118, -0.9464563060452031)),
                           (volmdlr.Point3D(9.80579460566209, 6.917108655469294, 0.0),
                            volmdlr.Point3D(12.0, 0.0, 0.6956521739130439),
                            volmdlr.Point3D(12.0, 0.0, 0.5217391304347831)),
                           (volmdlr.Point3D(8.773150355532506, 8.18729704110092, 0.0),
                            volmdlr.Point3D(12.0, 0.0, 0.8695652173913047),
                            volmdlr.Point3D(12.0, 0.0, 0.6956521739130439))
                           ]:

            surface = surfaces.Plane3D.from_3_points(p1, p2, p3)
            self.assertTrue(surface.point_on_surface(p1))
            self.assertTrue(surface.point_on_surface(p2))
            self.assertTrue(surface.point_on_surface(p3))
            self.assertAlmostEqual(surface.frame.w.dot(p1 - p2), 0.)
            self.assertAlmostEqual(surface.frame.w.dot(p3 - p2), 0.)
            self.assertAlmostEqual(surface.frame.w.dot(p3 - p1), 0.)

if __name__ == '__main__':
    unittest.main()
