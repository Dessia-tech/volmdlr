import math
import unittest

import volmdlr
from volmdlr import edges, faces, surfaces


class TestCylindricalFace3D(unittest.TestCase):
    cylindrical_surface1 = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 0.32)
    cylindrical_face1 = faces.CylindricalFace3D.from_surface_rectangular_cut(
        cylindrical_surface1, -0.01, 1.3, -0.1, 0.3)

    cylindrical_surface2 = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 12.0)
    cylindrical_face2 = faces.CylindricalFace3D.from_surface_rectangular_cut(cylindrical_surface2, 0, 3.14, 0., 8.)

    def test_linesegment_intersections(self):
        lineseg3d = edges.LineSegment3D(volmdlr.O3D, volmdlr.Point3D(0.3, 0.3, .3))
        line_inters = self.cylindrical_face1.linesegment_intersections(lineseg3d)
        self.assertEqual(len(line_inters), 1)
        self.assertTrue(line_inters[0].is_close(volmdlr.Point3D(0.22627416, 0.22627416, 0.22627416)))
        cylindrical_face1 = faces.CylindricalFace3D.from_surface_rectangular_cut(
            self.cylindrical_surface1, -math.pi / 4, 1.5 * math.pi, -0.1, 0.3)
        lineseg3d_2 = edges.LineSegment3D(volmdlr.Point3D(-0.3, -0.3, -.1), volmdlr.Point3D(0.3, 0.3, .3))
        line_inters_2 = cylindrical_face1.linesegment_intersections(lineseg3d_2)
        self.assertEqual(len(line_inters_2), 2)
        self.assertTrue(line_inters_2[0].is_close(volmdlr.Point3D(0.2262741, 0.2262741, 0.2508494)))
        self.assertTrue(line_inters_2[1].is_close(volmdlr.Point3D(-0.2262741, -0.2262741, -0.0508494)))

    def test_triangulation_quality(self):
        """
        The triangle middle of triangulation should be at least at radius/20 of the surface
        """
        triangulation = self.cylindrical_face2.triangulation()
        # ax = self.cylindrical_surface2.plot()
        for i1, i2, i3 in triangulation.triangles:
            point1 = triangulation.points[i1]
            point2 = triangulation.points[i2]
            point3 = triangulation.points[i3]

            triangle = faces.Triangle3D(point1, point2, point3)
            # triangle.plot(ax=ax)
            # Test orthogonality
            self.assertAlmostEqual(triangle.surface3d.frame.w.dot(point1 - point2), 0.)
            # Test distance from middle to surface

            # if self.cylindrical_surface2.point_distance(triangle.middle()) < self.cylindrical_surface2.radius*0.05:
            #     triangle.middle().plot(ax=ax, color='g')
            # else:
            #     triangle.middle().plot(ax=ax, color='r')

            self.assertLess(self.cylindrical_surface2.point_distance(triangle.middle()),
                            self.cylindrical_surface2.radius * 0.05)


if __name__ == '__main__':
    unittest.main()
