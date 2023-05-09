import unittest
import volmdlr
from volmdlr import surfaces, wires, edges


class TestSphericalSurface3D(unittest.TestCase):
    surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)

    def test_contour3d_to_2d(self):
        surface = surfaces.SphericalSurface3D.load_from_file("surfaces/objects_spherical_tests/sphericalsurface1.json")
        contour = wires.Contour3D.load_from_file("surfaces/objects_spherical_tests/spericalsurface1_contour0.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 6)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 4.107527949001648, 2)

    def arc3d_to_2d(self):
        arc_with_two_singularities = edges.Arc3D.load_from_file(
            "surfaces/objects_spherical_tests/arc_with_two_singularities.json")
        test_arc = self.surface3d.arc3d_to_2d(arc_with_two_singularities)
        self.assertEqual(len(test_arc), 5)
        self.assertTrue(test_arc[0].start.is_close(volmdlr.Point2D(0, 0)))
        self.assertTrue(test_arc[-1].end.is_close(volmdlr.Point2D(0.0, 0.8975979010256554)))


if __name__ == '__main__':
    unittest.main()
