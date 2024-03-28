import unittest
import os
import volmdlr
import volmdlr.edges as vme
from volmdlr import curves
from dessia_common.core import DessiaObject


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'bsplinecurve_objects')


class TestBSplineCurve3D(unittest.TestCase):
    b_splinecurve3d = vme.BSplineCurve3D(degree=5, control_points=[
        volmdlr.Point3D(0.5334, 4.61e-10, -2.266), volmdlr.Point3D(0.5334, 0.236642912449, -2.26599999893),
        volmdlr.Point3D(0.5334, 0.473285829931, -2.23144925183),
        volmdlr.Point3D(0.5334, 0.70316976404, -2.16234807551),
        volmdlr.Point3D(0.5334, 1.13611540546, -1.95904362568), volmdlr.Point3D(0.5334, 1.49286052971, -1.64044168585),
        volmdlr.Point3D(0.5334, 1.64654439419, -1.45604332404), volmdlr.Point3D(0.5334, 1.77109261028, -1.25188280667),
        volmdlr.Point3D(0.5334, 1.86385510975, -1.03417888209)], knot_multiplicities=[6, 3, 6],
                                         knots=[0.0, 0.4999999725155696, 1.0])

    def test_bounding_box(self):
        bspline = vme.BSplineCurve3D.from_points_interpolation([
            volmdlr.Point3D(1.0, 1.0, 0.0),
            volmdlr.Point3D(0.8090169943749475, 0.8090169943749475, 0.587785252292473),
            volmdlr.Point3D(0.30901699437494745, 0.30901699437494745, 0.9510565162951533),
            volmdlr.Point3D(0.0, 0.0, 1.0),
            volmdlr.Point3D(-0.30901699437494734, -0.30901699437494734, 0.9510565162951533),
            volmdlr.Point3D(-0.8090169943749473, -0.8090169943749473, 0.587785252292473),
            volmdlr.Point3D(-1.0, -1.0, 0.0)], 2, centripetal=True)
        bbox = bspline.bounding_box
        self.assertAlmostEqual(bbox.volume(), 3.9998457739878495, 3)

    def test_trim(self):
        obj = vme.BSplineCurve3D.from_json(os.path.join(folder, "bspline_buggy_trim.json"))
        point1 = volmdlr.Point3D(1.20555954308, -0.879118549155, 0.938030639643)
        point2 = volmdlr.Point3D(1.2150653573, -0.879118549155, 0.958332154591)
        trimmed_curve = obj.trim(point1, point2)
        self.assertTrue(trimmed_curve.start.is_close(point1))
        self.assertTrue(trimmed_curve.end.is_close(point2))
        self.assertAlmostEqual(trimmed_curve.length(), 0.03513727259692126, 2)
        obj = vme.BSplineCurve3D.from_json(os.path.join(folder, "bsplinecurve_trim.json"))
        point1 = volmdlr.Point3D(0.342947999551, -0.440408114191, 0.0132802444727)
        point2 = volmdlr.Point3D(0.342919095763, -0.44741803835000005, 0.0132953396808)
        trimmed_curve = obj.trim(point1, point2, True)
        self.assertTrue(trimmed_curve.start.is_close(point1))
        self.assertTrue(trimmed_curve.end.is_close(point2))
        self.assertAlmostEqual(trimmed_curve.length(), 0.011010880733091775, 2)
        bspline, point1, point2 = DessiaObject.from_json(
            os.path.join(folder, "test_bspline_trim271123.json")).primitives
        trim = bspline.trim(point1, point2, True)
        self.assertAlmostEqual(trim.length(), 14.607916441075464)
        trim = bspline.trim(point2, point1, True)
        self.assertAlmostEqual(trim.length(), 2.5461209947115186)

        bspline = vme.BSplineCurve3D.from_json(os.path.join(folder, "test_periodic_bspline_trim.json"))

        point1 = bspline.point_at_abscissa(1)
        point2 = bspline.point_at_abscissa(3)

        for pt1, pt2 in [(bspline.start, point1), (point1, bspline.start), (point2, bspline.end), (point1, point2),
                         (point2, point1)]:
            trim = bspline.trim(pt1, pt2)
            self.assertTrue(trim.start.is_close(pt1))
            self.assertTrue(trim.end.is_close(pt2))
        trim = bspline.trim(bspline.start, bspline.end)
        self.assertEqual(bspline, trim)

        bspline = vme.BSplineCurve3D.from_json(os.path.join(folder, "bsplinecurve3d_split_test.json"))
        point1 = volmdlr.Point3D(0.0781678147963, -0.08091364816680001, 0.112275939295)
        point2 = volmdlr.Point3D(0.0711770282536, -0.08091364816680001, 0.11191690794)
        trimmed_curve = bspline.trim(point1, point2, True)
        self.assertTrue(trimmed_curve.start.is_close(point1))
        self.assertAlmostEqual(bspline.point_to_parameter(trimmed_curve.end), 0.49999, 4)

        bspline = vme.BSplineCurve3D.from_json(os.path.join(folder, "bsplinecurve_close_points_trim_bug.json"))
        point1 = volmdlr.Point3D(-0.544134241655, -0.609582655058, 0.271301659325)
        point2 = volmdlr.Point3D(-0.544132239261, -0.609587093354, 0.271301378758)
        trimmed_curve = bspline.trim(point1, point2, True)
        self.assertEqual(len(trimmed_curve.knots), 2)
        self.assertAlmostEqual(trimmed_curve.length(), 4.876112145226163e-06)

        bspline = vme.BSplineCurve3D.from_json(
            os.path.join(folder, "bsplinecurve_trim_with_close_to_bound_points.json"))
        point1 = volmdlr.Point3D(-0.662347482412, 0.174584944052, 0.484523514816)
        point2 = volmdlr.Point3D(-0.66214998812, 0.174699062854, 0.484497837201)
        trimmed_curve = bspline.trim(point1, point2, True)
        self.assertEqual(len(trimmed_curve.knots), 4)
        self.assertAlmostEqual(trimmed_curve.length(), 0.0002325755440461709, 6)

        bspline = vme.BSplineCurve3D.from_json(
            os.path.join(folder, "bsplinecurve_split_bug_2.json"))
        point1 = volmdlr.Point3D(3.7101740508972862, -1.631997806124392e-05, -18.5000000067868)
        point2 = volmdlr.Point3D(3.493253085347731, 1.2499999999999991, -18.49999999999689)
        trimmed_curve = bspline.trim(point1, point2, True)
        self.assertTrue(bspline.start.is_close(trimmed_curve.start))
        self.assertTrue(bspline.end.is_close(trimmed_curve.end))

    def test_from_step(self):
        obj_list = volmdlr.model.VolumeModel.from_json(
            os.path.join(folder, "periodic_bsplinecurve_from_step_test_object_dict.json")).primitives
        object_dict = {0: obj_list[0], 1: obj_list[1], 2: obj_list[2]}
        arguments = ["''", 1, 2, 0, '.F.']
        bsplinecurve = vme.Edge.from_step(arguments, object_dict)
        self.assertTrue(bsplinecurve.start.is_close(object_dict[1], 1e-5))
        self.assertTrue(bsplinecurve.end.is_close(object_dict[2], 1e-5))
        self.assertTrue(bsplinecurve.point_at_abscissa(0.5 * bsplinecurve.length()).is_close(
            volmdlr.Point3D(0.04915593260514362, -0.04264529220680001, 0.14332598788877735)))

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
        bspline = vme.BSplineCurve3D.from_points_interpolation(points, 3, centripetal=True)
        linesegment = vme.LineSegment3D(volmdlr.Point3D(-3.0, 4.0, 1.0), volmdlr.Point3D(-3, -3, 0))
        intersections = bspline.linesegment_intersections(linesegment)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0].is_close(
            volmdlr.Point3D(-3.0, 1.0222488954206392, 0.5746069850600913)))

    def test_point_at_abscissa(self):
        bspline = vme.BSplineCurve3D.from_json(os.path.join(folder, "bsplinecurve_periodic.json"))
        self.assertTrue(bspline.start.is_close(bspline.point_at_abscissa(0)))
        self.assertTrue(bspline.end.is_close(bspline.point_at_abscissa(bspline.length())))
        self.assertTrue(bspline.point_at_abscissa(0.5 * bspline.length()).is_close(
            volmdlr.Point3D(0.3429479995510001, -0.44040811419137504, 0.01328024447265125)))

    def test_decompose(self):
        bspline = vme.BSplineCurve3D.from_json(os.path.join(folder, "spiral_bsplinecurve.json"))
        decompose_results = list(bspline.decompose(return_params=True))
        self.assertEqual(len(decompose_results), 37)
        for patch, param in decompose_results:
            self.assertTrue(bspline.evaluate_single(param[0]).is_close(patch.start))
            self.assertTrue(bspline.evaluate_single(param[1]).is_close(patch.end))
        bezier_patches = list(bspline.decompose())
        self.assertEqual(len(bezier_patches), 37)

    def test_line_intersections(self):
        line = curves.Line3D(volmdlr.Point3D(0.5334, -0.44659009801843536, 0.0),
                          volmdlr.Point3D(0.5334, 0.4342689853571558, -0.47337857496375274))
        bspline_line_intersections = self.b_splinecurve3d.line_intersections(line)
        self.assertTrue(bspline_line_intersections[0].is_close(
            volmdlr.Point3D(0.5334, 1.784620497933768, -1.1990649949459866)))

    def test_linesegment_intersection(self):
        linesegment1 = vme.LineSegment3D(volmdlr.Point3D(0.5334, -0.44659009801843536, 0.0),
                                         volmdlr.Point3D(0.5334, 0.4342689853571558, -0.47337857496375274))
        linesegment2 = vme.LineSegment3D(volmdlr.Point3D(0.5334, -0.44659009801843536, 0.0),
                                         volmdlr.Point3D(0.5334, 2.1959871521083385, -1.4201357248912583))
        bspline_lineseg_intersections1 = self.b_splinecurve3d.linesegment_intersections(linesegment1)
        bspline_lineseg_intersections2 = self.b_splinecurve3d.linesegment_intersections(linesegment2)
        self.assertFalse(bspline_lineseg_intersections1)
        self.assertTrue(bspline_lineseg_intersections2[0].is_close(
            volmdlr.Point3D(0.5334, 1.7846204999239552, -1.1990649960155242)))

    def test_normal(self):
        normal = self.b_splinecurve3d.normal()
        self.assertTrue(normal.is_close(volmdlr.Z3D))

    def test_abscissa(self):
        point = volmdlr.Point3D(0.18357300891283804, 0.7465725481678318, 0.44333916797214895)
        bsplinecurve = vme.BSplineCurve3D.from_json(os.path.join(folder, "bsplinecurve3d_abscissa_test.json"))
        abscissa = bsplinecurve.abscissa(point)
        self.assertTrue(bsplinecurve.point_at_abscissa(abscissa).is_close(point))

    def test_local_discretization(self):
        edge, start, end = DessiaObject.from_json(os.path.join(
            folder, 'test_bspline_local_discretizations.json')).primitives
        expected_points = [volmdlr.Point3D(0.40000000000000013, 0.3055497472688364, -0.0577802904293785),
                           volmdlr.Point3D(0.39999999999999997, 0.3054558765193882, -0.05729389615629579),
                           volmdlr.Point3D(0.4, 0.3053650847805137, -0.056818942485284164),
                           volmdlr.Point3D(0.4000000000000001, 0.30527736613101497, -0.05635539982690565),
                           volmdlr.Point3D(0.4, 0.3051925979747842, -0.055902690345743224),
                           volmdlr.Point3D(0.4000000000000001, 0.3051105616115585, -0.055459784623149766),
                           volmdlr.Point3D(0.4, 0.3050309628076735, -0.055025298319996994),
                           volmdlr.Point3D(0.39999999999999997, 0.3049534523668191, -0.054597588839424546),
                           volmdlr.Point3D(0.39999999999999997, 0.30487764670079315, -0.05417485198958898),
                           volmdlr.Point3D(0.4000000000000001, 0.3048031484002564, -0.05375521864641269)]
        discretized_points_between_1_2 = edge.local_discretization(start, end, 10)
        self.assertEqual(len(discretized_points_between_1_2), len(expected_points))
        for result, expected_point in zip(discretized_points_between_1_2, expected_points):
            self.assertTrue(result.is_close(expected_point))

    def test_move_frame_along(self):
        degree = 5
        control_points = [
            volmdlr.Point3D(-1, 0, 0),
            volmdlr.Point3D(0.3, 0.2, 0.1),
            volmdlr.Point3D(0.5, -0.1, 0.4),
            volmdlr.Point3D(0.5, -0.4, 0.0),
            volmdlr.Point3D(-0.1, -0.2, -0.3),
            volmdlr.Point3D(-0.3, 0.4, 0.1)]
        knots = [0.0, 1.0]
        knot_multiplicities = [6, 6]
        weights = None  # [1, 2, 1, 2, 1, 2]
        bspline_curve3d = vme.BSplineCurve3D(degree=degree,
                                             control_points=control_points,
                                             knot_multiplicities=knot_multiplicities,
                                             knots=knots,
                                             weights=weights,
                                             name='B Spline Curve 3D 1')
        frame = volmdlr.Frame3D(
            origin=volmdlr.Point3D(-1.0, 0.0, 0.0),
            u=volmdlr.Vector3D(0.16926811079722504, -0.8799271690658463, -0.44393297220065076),
            v=volmdlr.Vector3D(1.3877787807814457e-17, 0.4504326973383234, -0.8928103858986645),
            w=volmdlr.Vector3D(0.9855700414821559, 0.15112432732120784, 0.076243891719756))
        self.assertEqual(bspline_curve3d.move_frame_along(frame),
                         volmdlr.Frame3D(origin=volmdlr.Point3D(-0.3, 0.4, 0.1),
                                         u=volmdlr.Vector3D(0.9632101019095523, 0.22376214189648538, 0.14885161548766362),
                                         v=volmdlr.Vector3D(2.7755575615628914e-17, 0.553867484333956, -0.8326048341185484),
                                         w=volmdlr.Vector3D(-0.26874951084493176, 0.8019733871217126, 0.5334907560296971)))


if __name__ == '__main__':
    unittest.main()
