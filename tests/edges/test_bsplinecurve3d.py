import unittest
import volmdlr
import volmdlr.edges as vme


class TestBSplineCurve3D(unittest.TestCase):
    def test_trim(self):
        obj = vme.BSplineCurve3D.load_from_file(r"C:\Users\gabri\Documents\dessia\tests\bspline_buggy_trim.json")
        point1 = volmdlr.Point3D(1.20555954308, -0.879118549155, 0.938030639643)
        point2 = volmdlr.Point3D(1.2150653573, -0.879118549155, 0.958332154591)
        trimmed_curve = obj.trim(point1, point2)
        self.assertTrue(trimmed_curve.start.is_close(point1))
        self.assertAlmostEqual(trimmed_curve.length(), 0.03513727259692126, 2)


if __name__ == '__main__':
    unittest.main()
