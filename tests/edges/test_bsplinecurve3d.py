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

    def test_from_step(self):
        obj_list = volmdlr.core.VolumeModel.load_from_file(
            "edges/bsplinecurve_objects/periodic_bsplinecurve_from_step_test_object_dict.json").primitives
        object_dict = {0: obj_list[0], 1: obj_list[1], 2: obj_list[2]}
        arguments = ["''", 1, 2, 0, '.F.']
        bsplinecurve = vme.Edge.from_step(arguments, object_dict)
        self.assertTrue(bsplinecurve.start.is_close(object_dict[2], 1e-5))
        self.assertTrue(bsplinecurve.end.is_close(object_dict[1], 1e-5))


if __name__ == '__main__':
    unittest.main()
