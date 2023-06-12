import unittest
import volmdlr
import volmdlr.edges as vme


class TestBSplineCurve3D(unittest.TestCase):
    def test_trim(self):
        obj = vme.BSplineCurve3D.load_from_file("edges/bsplinecurve_objects/bspline_buggy_trim.json")
        point1 = volmdlr.Point3D(1.20555954308, -0.879118549155, 0.938030639643)
        point2 = volmdlr.Point3D(1.2150653573, -0.879118549155, 0.958332154591)
        trimmed_curve = obj.trim(point1, point2)
        self.assertTrue(trimmed_curve.start.is_close(point1))
        self.assertAlmostEqual(trimmed_curve.length(), 0.03513727259692126, 2)

    def test_trim(self):
        bsplinecurve = vme.BSplineCurve3D.load_from_file("edges/bsplinecurve_objects/bsplinecurve_split_test.json")
        point1 = volmdlr.Point3D(0.024529919518700004, -0.0825, 0.0528937378742)
        point2 = volmdlr.Point3D(0.0187137575132, -0.0825, 0.0556938541471)
        new_bspline = bsplinecurve.trim(point1, point2)
        self.assertTrue(point1.is_close(new_bspline.start))
        self.assertTrue(point2.is_close(new_bspline.end))



if __name__ == '__main__':
    unittest.main()
