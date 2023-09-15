import unittest
import volmdlr
import volmdlr.edges as vme


class TestBSplineCurve3D(unittest.TestCase):

    def test_bounding_box(self):
        bspline = vme.BSplineCurve3D.from_points_interpolation([
            volmdlr.Point3D(1.0, 1.0, 0.0),
            volmdlr.Point3D(0.8090169943749475, 0.8090169943749475, 0.587785252292473),
            volmdlr.Point3D(0.30901699437494745, 0.30901699437494745, 0.9510565162951533),
            volmdlr.Point3D(0.0, 0.0, 1.0),
            volmdlr.Point3D(-0.30901699437494734, -0.30901699437494734, 0.9510565162951533),
            volmdlr.Point3D(-0.8090169943749473, -0.8090169943749473, 0.587785252292473),
            volmdlr.Point3D(-1.0, -1.0, 0.0)], 2)
        bbox = bspline.bounding_box
        self.assertAlmostEqual(bbox.volume(), 4.0, 3)

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
        self.assertTrue(bsplinecurve.start.is_close(object_dict[1], 1e-5))
        self.assertTrue(bsplinecurve.end.is_close(object_dict[2], 1e-5))
        self.assertTrue(bsplinecurve.point_at_abscissa(0.5*bsplinecurve.length()).is_close(
            volmdlr.Point3D(0.04916207191208274, -0.0426452922068, 0.14332757998779702)))

    def test_bspline_linesegment_minimum_distance(self):
        points = [volmdlr.Point3D(1.2918566581549966, 2.3839907440191492, 0.5678759590090421),
                  volmdlr.Point3D(1.2067665579541171, -1.246879774203074, -0.4359328108960321),
                  volmdlr.Point3D(-1.2905737351068276, -5.961765089244547, -0.9872550297481824),
                  volmdlr.Point3D(7.33260591629263, -4.272128323147327, -0.4240427743824422),
                  volmdlr.Point3D(7.115095014105684, 0.40888620982702983, 1.1362954032756774),
                  volmdlr.Point3D(-3.0, 1.022248896290622, 0.5746069851843745),
                  volmdlr.Point3D(2.739350840642852, -5.869347626045908, -0.7880999427201254)]
        bspline = vme.BSplineCurve3D.from_points_interpolation(points, 3)
        linesegment = vme.LineSegment3D(volmdlr.Point3D(-3.0, 4.0, 1.0), volmdlr.Point3D(-3, -3, 0))
        dist, min_dist_pt1, min_dist_pt2 = bspline.minimum_distance(linesegment, True)
        self.assertAlmostEqual(dist, 4.155003073325757e-09)
        self.assertTrue(min_dist_pt1.is_close(min_dist_pt2))

    def test_bspline_linesegment_intersections(self):
        points = [volmdlr.Point3D(1.2918566581549966, 2.3839907440191492, 0.5678759590090421),
                  volmdlr.Point3D(1.2067665579541171, -1.246879774203074, -0.4359328108960321),
                  volmdlr.Point3D(-1.2905737351068276, -5.961765089244547, -0.9872550297481824),
                  volmdlr.Point3D(7.33260591629263, -4.272128323147327, -0.4240427743824422),
                  volmdlr.Point3D(7.115095014105684, 0.40888620982702983, 1.1362954032756774),
                  volmdlr.Point3D(-3.0, 1.022248896290622, 0.5746069851843745),
                  volmdlr.Point3D(2.739350840642852, -5.869347626045908, -0.7880999427201254)]
        bspline = vme.BSplineCurve3D.from_points_interpolation(points, 3)
        linesegment = vme.LineSegment3D(volmdlr.Point3D(-3.0, 4.0, 1.0), volmdlr.Point3D(-3, -3, 0))
        intersections = bspline.linesegment_intersections(linesegment)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0].is_close(
            volmdlr.Point3D(-3.000003493713931, 1.022247107552729, 0.5746061300812078)))


if __name__ == '__main__':
    unittest.main()
