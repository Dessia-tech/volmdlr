import math
import unittest

import volmdlr
from volmdlr import edges, faces, surfaces, wires
from dessia_common.core import DessiaObject


class TestToroidalFace3D(unittest.TestCase):
    surface1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 0.32, 0.08)
    face1 = faces.ToroidalFace3D.from_surface_rectangular_cut(surface1, -0.1, 1.3, 2, 0.3)

    def test_triangulation_quality(self):
        """
        The triangle middle of triangulation should be at least at radius/20 of the surface
        """
        triangulation = self.face1.triangulation()
        for i1, i2, i3 in triangulation.triangles:
            point1 = triangulation.points[i1]
            point2 = triangulation.points[i2]
            point3 = triangulation.points[i3]

            triangle = faces.Triangle3D(point1, point2, point3)
            # Test orthogonality
            self.assertAlmostEqual(triangle.surface3d.frame.w.dot(point1 - point2), 0.)
            # Test distance from middle to surface

            self.assertLess(self.surface1.point_distance(triangle.middle()),
                            self.surface1.small_radius * 0.05)

    def test_number_triangles(self):
        triangulation = self.face1.triangulation()
        triangulation.plot()
        n_triangles = len(triangulation.triangles)
        n_triangles_max = 220  # Could be 208 (13*8 tiles on this ex, 2 triangles per tile)
        self.assertLess(n_triangles, n_triangles_max,
                        f'Too much triangles in cylindrical face triangulation: {n_triangles}/{n_triangles_max}')


if __name__ == '__main__':
    unittest.main()
