import unittest
import volmdlr
import volmdlr.edges as vme


class TestFullArcEllipse2D(unittest.TestCase):
    start_end = volmdlr.Point2D(0.0225, 0.0)
    major_axis = 0.0225
    minor_axis = 0.0075
    center = volmdlr.O2D
    major_dir = volmdlr.X2D
    ellipse = vme.FullArcEllipse2D(start_end, major_axis, minor_axis, center, major_dir)

    def test_init(self):
        self.assertAlmostEqual(self.ellipse.major_axis, 0.0225, places=4)
        self.assertAlmostEqual(self.ellipse.minor_axis, 0.0075, places=4)
        self.assertAlmostEqual(self.ellipse.angle, volmdlr.TWO_PI, places=4)

    def test_length(self):
        self.assertAlmostEqual(self.ellipse.length(), 0.10023669584870037)

    def test_to_3d(self):
        plane_origin = volmdlr.Point3D(1, 1, 1)
        x = volmdlr.Y3D
        y = volmdlr.Z3D
        ellipse3d = self.ellipse.to_3d(plane_origin, x, y)
        self.assertEqual(ellipse3d.major_dir, x)
        self.assertEqual(ellipse3d.minor_dir, y)
        self.assertEqual(ellipse3d.normal, volmdlr.X3D)
        self.assertAlmostEqual(ellipse3d.major_axis, 0.0225, places=4)
        self.assertAlmostEqual(ellipse3d.minor_axis, 0.0075, places=4)

    def test_reverse(self):
        self.assertEqual(self.ellipse, self.ellipse.reverse())

    def test_frame_mapping(self):
        new_frame = volmdlr.Frame2D(volmdlr.O2D, -volmdlr.Y2D, volmdlr.X2D)
        new_ellipse = self.ellipse.frame_mapping(new_frame, 'new')
        self.assertEqual(new_ellipse.major_dir, volmdlr.Vector2D(0.0, 1.0))
        self.assertEqual(new_ellipse.minor_dir, volmdlr.Vector2D(-1.0, 0.0))

    def test_abscissa(self):
        point1 = volmdlr.Point2D(0, -0.0075)
        point2 = volmdlr.Point2D(0.0225, 0)
        self.assertAlmostEqual(self.ellipse.abscissa(point1), 0.75*self.ellipse.length())
        self.assertAlmostEqual(self.ellipse.abscissa(point2), self.ellipse.length())

    def test_translation(self):
        translated_ellipse = self.ellipse.translation(volmdlr.X2D)
        self.assertEqual(translated_ellipse.center, volmdlr.Point2D(1, 0))
        self.assertEqual(translated_ellipse.start_end, volmdlr.Point2D(1.0225, 0))


if __name__ == '__main__':
    unittest.main()
