import math
import unittest
import volmdlr
import volmdlr.edges as vme


class TestFullArcEllipse3D(unittest.TestCase):
    start_end = volmdlr.Point3D(0.0225, 0.0, 0.0)
    major_axis = 0.0225
    minor_axis = 0.0075
    center = volmdlr.O3D
    major_dir = volmdlr.X3D
    normal = volmdlr.Z3D
    ellipse = vme.FullArcEllipse3D(start_end, major_axis, minor_axis, center, normal, major_dir)

    def test_init(self):
        self.assertAlmostEqual(self.ellipse.major_axis, 0.0225, places=4)
        self.assertEqual(self.ellipse.minor_axis, 0.0075)
        self.assertAlmostEqual(self.ellipse.theta, 0, places=4)

    def test_length(self):
        self.assertAlmostEqual(self.ellipse.length(), 0.10023669584870037)

    def test_to_2d(self):
        plane_origin = volmdlr.Point3D(1, 1, 0)
        x = volmdlr.Vector3D(0.5*math.sqrt(2), 0.5*math.sqrt(2), 0)
        y = volmdlr.Vector3D(-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0)
        ellipse2d = self.ellipse.to_2d(plane_origin, x, y)
        self.assertTrue(ellipse2d.major_dir.is_close(volmdlr.Vector2D(0.5*math.sqrt(2), -0.5*math.sqrt(2))))
        self.assertTrue(ellipse2d.minor_dir.is_close(volmdlr.Vector2D(0.5*math.sqrt(2), 0.5*math.sqrt(2))))
        self.assertAlmostEqual(ellipse2d.major_axis, 0.0225, places=4)
        self.assertAlmostEqual(ellipse2d.minor_axis, 0.0075, places=4)

    def test_reverse(self):
        self.assertEqual(self.ellipse, self.ellipse.reverse())

    def test_frame_mapping(self):
        new_frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.Z3D, volmdlr.X3D, -volmdlr.Y3D)
        new_ellipse = self.ellipse.frame_mapping(new_frame, 'new')
        self.assertEqual(new_ellipse.major_dir, volmdlr.Vector3D(0.0, 1.0, 0.0))
        self.assertEqual(new_ellipse.minor_dir, volmdlr.Vector3D(0.0, 0.0, 1.0))
        self.assertEqual(new_ellipse.normal, volmdlr.Vector3D(1.0, 0.0, 0.0))

    def test_abscissa(self):
        point1 = volmdlr.Point3D(0, -0.0075, 0)
        point2 = volmdlr.Point3D(0.0225, 0, 0)
        self.assertAlmostEqual(self.ellipse.abscissa(point1), 0.75*self.ellipse.length())
        self.assertAlmostEqual(self.ellipse.abscissa(point2), self.ellipse.length())

    def test_translation(self):
        translated_ellipse = self.ellipse.translation(volmdlr.X3D)
        self.assertEqual(translated_ellipse.center, volmdlr.Point3D(1, 0, 0))
        self.assertEqual(translated_ellipse.start_end, volmdlr.Point3D(1.0225, 0, 0))

    def test_point_belongs(self):
        ellipse = vme.FullArcEllipse3D.load_from_file("edges/fullarcellipse3d_point_belongs.json")
        point = volmdlr.Point3D(-0.45404913959118903, -0.5072162164760068, 0.5060526668248432)
        self.assertTrue(ellipse.point_belongs(point, 1e-5))

if __name__ == '__main__':
    unittest.main()
