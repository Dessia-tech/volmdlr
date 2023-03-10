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

    def test_reverse(self):
        self.assertEqual(self.ellipse, self.ellipse.reverse())

    def test_frame_mapping(self):
        new_frame = volmdlr.Frame2D(volmdlr.O2D, -volmdlr.Y2D, volmdlr.X2D)
        new_ellipse = self.ellipse.frame_mapping(new_frame, 'new')
        self.assertEqual(new_ellipse.major_dir, volmdlr.Vector2D(0.0, 1.0))
        self.assertEqual(new_ellipse.minor_dir, volmdlr.Vector2D(-1.0, 0.0))


if __name__ == '__main__':
    unittest.main()
