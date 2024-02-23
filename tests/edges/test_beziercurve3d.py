"""
Unit tests for volmdlr.faces.BSplineCurve
"""
import unittest
import volmdlr
import volmdlr.edges as vme


class TestBezierCurve3D(unittest.TestCase):
    # Set up the Bézier curve
    degree = 3
    ctrlpts = [volmdlr.Point3D(0, 0, 0),
              volmdlr.Point3D(1, 1, 2),
              volmdlr.Point3D(2, 1, 1),
              volmdlr.Point3D(3, 0, 4)]
    bezier_curve3d = vme.BezierCurve3D(degree=degree,
                                       control_points=ctrlpts,
                                       name='bezier curve 1')

    def test_cut_before(self):
        new_bezier = self.bezier_curve3d.cut_before(0.5)
        self.assertTrue(new_bezier.start.is_close(self.bezier_curve3d.evaluate_single(0.5)))
        self.assertTrue(new_bezier.end.is_close(self.bezier_curve3d.end))


if __name__ == '__main__':
    unittest.main()
