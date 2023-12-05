import unittest
import os
import volmdlr
import volmdlr.edges as vme
from dessia_common.core import DessiaObject

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'bsplinecurve_objects')


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
        obj = vme.BSplineCurve3D.load_from_file(os.path.join(folder, "bspline_buggy_trim.json"))
        point1 = volmdlr.Point3D(1.20555954308, -0.879118549155, 0.938030639643)
        point2 = volmdlr.Point3D(1.2150653573, -0.879118549155, 0.958332154591)
        trimmed_curve = obj.trim(point1, point2)
        self.assertTrue(trimmed_curve.start.is_close(point1))
        self.assertAlmostEqual(trimmed_curve.length(), 0.03513727259692126, 2)
        obj = vme.BSplineCurve3D.load_from_file(os.path.join(folder, "bsplinecurve_trim.json"))
        point1 = volmdlr.Point3D(0.342947999551, -0.440408114191, 0.0132802444727)
        point2 = volmdlr.Point3D(0.342919095763, -0.44741803835000005, 0.0132953396808)
        trimmed_curve = obj.trim(point1, point2, True)
        self.assertTrue(trimmed_curve.start.is_close(point1))
        self.assertAlmostEqual(trimmed_curve.length(), 0.011010880733091775, 2)
        bspline, point1, point2 = DessiaObject.load_from_file(
            os.path.join(folder, "test_bspline_trim271123.json")).primitives
        trim = bspline.trim(point1, point2, True)
        self.assertAlmostEqual(trim.length(), 14.606177552397396)
        trim = bspline.trim(point2, point1, True)
        self.assertAlmostEqual(trim.length(), 2.5461209947115186)


    def test_from_step(self):
        obj_list = volmdlr.core.VolumeModel.load_from_file(
            os.path.join(folder, "periodic_bsplinecurve_from_step_test_object_dict.json")).primitives
        object_dict = {0: obj_list[0], 1: obj_list[1], 2: obj_list[2]}
        arguments = ["''", 1, 2, 0, '.F.']
        bsplinecurve = vme.Edge.from_step(arguments, object_dict)
        self.assertTrue(bsplinecurve.start.is_close(object_dict[1], 1e-5))
        self.assertTrue(bsplinecurve.end.is_close(object_dict[2], 1e-5))
        self.assertTrue(bsplinecurve.point_at_abscissa(0.5 * bsplinecurve.length()).is_close(
            volmdlr.Point3D(0.04916207192770078, -0.042645292206800016, 0.14332757999206563)))

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
            volmdlr.Point3D(-3.0, 1.0222488954206392, 0.5746069850600913)))

    def test_point_at_abscissa(self):
        bspline = vme.BSplineCurve3D.load_from_file(os.path.join(folder, "bsplinecurve_periodic.json"))
        self.assertTrue(bspline.start.is_close(bspline.point_at_abscissa(0)))
        self.assertTrue(bspline.end.is_close(bspline.point_at_abscissa(bspline.length())))
        self.assertTrue(bspline.point_at_abscissa(0.5 * bspline.length()).is_close(
            volmdlr.Point3D(0.3429479995510001, -0.44040811419137504, 0.01328024447265125)))


if __name__ == '__main__':
    unittest.main()
