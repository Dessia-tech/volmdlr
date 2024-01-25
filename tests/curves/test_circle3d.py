import os
import math
import unittest
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import edges, curves
from volmdlr.models.curves import circle3d


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'circle3D_objects')


class TestCircle3D(unittest.TestCase):
    list_points = circle3d.discretization_points(number_points=9)

    def test_discretization_points(self):
        expected_points = [volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                           volmdlr.Point3D(0.9855985596534886, 0.1195731558690501, 0.1195731558690501),
                           volmdlr.Point3D(0.8164965809277258, -0.4082482904638631, -0.4082482904638631),
                           volmdlr.Point3D(0.1691019787257626, -0.696923425058676, -0.696923425058676),
                           volmdlr.Point3D(-0.5773502691896257, -0.5773502691896258, -0.5773502691896258),
                           volmdlr.Point3D(-0.9855985596534886, -0.11957315586905026, -0.11957315586905026),
                           volmdlr.Point3D(-0.8164965809277259, 0.408248290463863, 0.408248290463863),
                           volmdlr.Point3D(-0.16910197872576277, 0.696923425058676, 0.696923425058676),
                           volmdlr.Point3D(0.5773502691896256, 0.577350269189626, 0.577350269189626)]
        for point, expected_point in zip(self.list_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_abscissa(self):
        expected_abscissas = [0, math.pi/4, math.pi/2, 3 * math.pi / 4, math.pi,
                              5 * math.pi / 4, 3 * math.pi / 2, 7 * math.pi / 4, 0]
        abscissas = [circle3d.abscissa(point) for point in self.list_points]
        for abscissa, expected_abscissa in zip(abscissas, expected_abscissas):
            self.assertAlmostEqual(abscissa, expected_abscissa)

    def test_point_at_abscissa(self):
        abscissas = [0, math.pi/4, math.pi/2, 3 * math.pi / 4, math.pi,
                     5 * math.pi / 4, 3 * math.pi / 2, 7 * math.pi / 4, 0]
        points_at_abcscissas = [circle3d.point_at_abscissa(abscissa) for abscissa in abscissas]
        for point, expected_point in zip(points_at_abcscissas, self.list_points):
            self.assertTrue(point.is_close(expected_point))

    def test_length(self):
        self.assertAlmostEqual(circle3d.length(), 6.283185307179586)

    def test_rotation(self):
        expected_points = [volmdlr.Point3D(0.0, -0.7071067811865475, 0.7071067811865477),
                           volmdlr.Point3D(0.7765343938240264, -0.606775209136424, -0.1697591846876029),
                           volmdlr.Point3D(0.4799246488165452, 0.33209907840941155, -0.812023727225957),
                           volmdlr.Point3D(-0.4799246488165446, 0.8120237272259567, -0.33209907840941194),
                           volmdlr.Point3D(-0.7765343938240268, 0.16975918468760318, 0.606775209136424),
                           volmdlr.Point3D(0.0, -0.7071067811865475, 0.7071067811865477)]
        rotated_circle3d = circle3d.rotation(volmdlr.O3D, circle3d.frame.v, math.pi / 2)
        rotated_circle3d_points = rotated_circle3d.discretization_points(number_points=6)
        for point, expected_point in zip(rotated_circle3d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_translation(self):
        expected_points = [volmdlr.Point3D(0.5773502691896258, 1.2844570503761734, -0.12975651199692173),
                           volmdlr.Point3D(0.9549454387105718, 0.497250629161079, -0.9169629332120162),
                           volmdlr.Point3D(0.01283846933518723, 5.8277296916986465e-05, -1.414155285076178),
                           volmdlr.Point3D(-0.9470108282979028, 0.4799829261134622, -0.934230636259633),
                           volmdlr.Point3D(-0.598123348937482, 1.2737850229851062, -0.140428539387989),
                           volmdlr.Point3D(0.5773502691896258, 1.2844570503761734, -0.12975651199692173)]
        translated_circle3d = circle3d.translation(circle3d.frame.w)
        translated_circle3d_points = translated_circle3d.discretization_points(number_points=6)
        for point, expected_point in zip(translated_circle3d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_frame_mapping(self):
        expected_points = [[volmdlr.Point3D(1.576, 0.29544000000000004, -0.527),
                            volmdlr.Point3D(1.576, 0.07432543819998319, -0.22266191478555086),
                            volmdlr.Point3D(1.576, -0.28344543819998314, -0.3389087192664086),
                            volmdlr.Point3D(1.576, -0.28344543819998314, -0.7150912807335914),
                            volmdlr.Point3D(1.576, 0.07432543819998312, -0.8313380852144492),
                            volmdlr.Point3D(1.576, 0.29544000000000004, -0.527)],
                           [volmdlr.Point3D(1.3359999999999999, 0.09544, -0.387),
                            volmdlr.Point3D(1.3359999999999999, -0.12567456180001682, -0.08266191478555085),
                            volmdlr.Point3D(1.3359999999999999, -0.4834454381999832, -0.19890871926640857),
                            volmdlr.Point3D(1.3359999999999999, -0.4834454381999832, -0.5750912807335914),
                            volmdlr.Point3D(1.3359999999999999, -0.1256745618000169, -0.6913380852144492),
                            volmdlr.Point3D(1.3359999999999999, 0.09544, -0.387)]]
        for i, side in enumerate(['old', 'new']):
            circle = curves.Circle3D(volmdlr.Frame3D(volmdlr.Point3D(1.456, -0.12456, -0.457),
                                                     volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D), 0.32)
            frame = volmdlr.Frame3D(volmdlr.Point3D(0.12, 0.1, -0.07), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
            framed_mapped_circle = circle.frame_mapping(frame, side)
            for expected_point, point in zip(expected_points[i],
                                             framed_mapped_circle.discretization_points(number_points=6)):
                self.assertTrue(expected_point.is_close(point))

    def test_to_2d(self):
        arc2d = circle3d.to_2d(circle3d.center, circle3d.frame.u, circle3d.frame.v)
        point2d_ = circle3d.point_at_abscissa(circle3d.length() * 0.5).to_2d(circle3d.center, circle3d.frame.u,
                                                                             circle3d.frame.v)
        self.assertTrue(arc2d.point_belongs(point2d_))

    def test_from_center_normal(self):
        from_center_normal = curves.Circle3D.from_center_normal(volmdlr.O3D, circle3d.normal, 1)
        self.assertTrue(from_center_normal.is_close(circle3d))

    def test_from_3_points(self):
        from_3_points = curves.Circle3D.from_3_points(*self.list_points[:3])
        self.assertTrue(from_3_points.is_close(circle3d))

    def test_extrusion(self):
        extrusion = circle3d.extrusion(circle3d.normal)
        extrusion_points = extrusion[0].outer_contour3d.discretization_points(number_points=10)
        expected_points = [volmdlr.Point3D(1.0, 0.0, 0.0),
                           volmdlr.Point3D(-0.047671222696162176, -0.7063028580314018, -0.7063028580314018),
                           volmdlr.Point3D(-0.9954549090533058, 0.06734064167230154, 0.06734064167230154),
                           volmdlr.Point3D(0.14258032800309864, 0.6998824365800763, 0.6998824365800763),
                           volmdlr.Point3D(1.0, 0.13488570125933452, -0.13488570125933444),
                           volmdlr.Point3D(0.6900466568659818, 1.2188858406194645, -0.19532772175363078),
                           volmdlr.Point3D(-0.7558373925663837, 1.1700915735925457, -0.24412198878054941),
                           volmdlr.Point3D(-0.6179832715397452, 0.15118561946618558, -1.2630279429069096),
                           volmdlr.Point3D(0.8147574288865314, 0.2971248717843101, -1.1170886905887851),
                           volmdlr.Point3D(1.0, 0.0, 0.0)]
        for point, expected_point in zip(extrusion_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_revolution(self):
        revolution_point = circle3d.center + circle3d.frame.u * 4 * circle3d.radius
        revolution = circle3d.revolution(revolution_point, circle3d.normal, math.pi / 4)
        rev_points = revolution[0].outer_contour3d.discretization_points(number_points=10)
        expected_points = [volmdlr.Point3D(-0.5773502691896257, -0.5773502691896257, -0.5773502691896257),
                           volmdlr.Point3D(-2.074836558779022, 0.6096275102008377, 0.6096275102008377),
                           volmdlr.Point3D(-2.440348177278158, 2.138739714140839, 1.3275800167324001),
                           volmdlr.Point3D(-0.6809780318545704, 2.129619885853125, 1.76359462600285),
                           volmdlr.Point3D(-2.0495252060301175, 1.1397182920111604, 2.421430993147704),
                           volmdlr.Point3D(-2.302089202598333, 0.9429901935690035, 0.9429901935690039),
                           volmdlr.Point3D(-0.9882719284073431, -0.3481685602407296, -0.3481685602407292),
                           volmdlr.Point3D(0.10025582212028905, -0.5961084181997304, 0.7966200624403085),
                           volmdlr.Point3D(0.3711135994842809, 0.9127888199039815, -0.1705616209354197),
                           volmdlr.Point3D(-0.5773502691896258, -0.5773502691896258, -0.5773502691896258)]
        for point, expected_point in zip(rev_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_point_belongs(self):
        for point in self.list_points:
            self.assertTrue(circle3d.point_belongs(point))
        self.assertFalse(circle3d.point_belongs(volmdlr.Point3D(2, 2, 2)))

    def test_trim(self):
        trim = circle3d.trim(self.list_points[2], self.list_points[5])
        self.assertAlmostEqual(trim.length(), 2.3561944901923444)
        trim = circle3d.trim(self.list_points[5], self.list_points[2])
        self.assertAlmostEqual(trim.length(), 3.9269908169872423)
        trim = circle3d.trim(self.list_points[5], self.list_points[5])
        self.assertAlmostEqual(trim.length(), 6.283185307179586)

    def test_linsegment_intersections(self):
        circle = curves.Circle3D(volmdlr.OXYZ, 1)
        lineseg = edges.LineSegment3D(volmdlr.Point3D(0, 0, -1), volmdlr.Point3D(2, 0, 1))
        line_seg2 = edges.LineSegment3D(volmdlr.Point3D(1, 2, 0), volmdlr.Point3D(-1, -2, 0))
        circle_linseg_intersections = circle.linesegment_intersections(lineseg) + circle.linesegment_intersections(
            line_seg2)
        self.assertEqual(len(circle_linseg_intersections), 3)
        expected_intersections = [volmdlr.Point3D(1.0, 0.0, 0.0),
                                  volmdlr.Point3D(-0.447213595499958, -0.894427190999916, 0.0),
                                  volmdlr.Point3D(0.447213595499958, 0.8944271909999157, 0.0)]
        for expected_point, point in zip(expected_intersections, circle_linseg_intersections):
            self.assertTrue(expected_point.is_close(point))
        circle2 = curves.Circle3D(volmdlr.OYZX, 1)
        lineseg2 = edges.LineSegment3D(volmdlr.Point3D(-1, 0, -2), volmdlr.Point3D(2, 0, 1))
        circle_linseg_intersections = circle2.linesegment_intersections(lineseg2)
        self.assertEqual(len(circle_linseg_intersections), 1)
        self.assertTrue(circle_linseg_intersections[0].is_close(volmdlr.Point3D(0.0, 0.0, -1.0)))
        circle3 = curves.Circle3D(
            volmdlr.Frame3D(volmdlr.O3D, volmdlr.Vector3D(0.2672612419124244, 0.5345224838248488, 0.8017837257372732),
                            volmdlr.Vector3D(0.9636241116594316, -0.14824986333222026, -0.22237479499833038),
                            volmdlr.Vector3D(0, 0.8320502943378438, -0.5547001962252291)), 1)
        circle_linseg_intersections1 = circle3.linesegment_intersections(lineseg2)
        self.assertEqual(len(circle_linseg_intersections1), 1)
        self.assertTrue(circle_linseg_intersections1[0].is_close(volmdlr.Point3D(1, 0.0, 0.0)))
        circle, lineseg = DessiaObject.from_json(os.path.join(folder,
            'test_circle_linesegment_intersections221223.json')).primitives
        intersections = circle.linesegment_intersections(lineseg)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0].is_close(volmdlr.Point3D(-1.9999993561471823, -0.5135128860482583, 0.9978935668376178)))

    def test_circle_intersections(self):
        circle1 = curves.Circle3D.from_3_points(volmdlr.Point3D(-3, -3, 0),
                                                volmdlr.Point3D(6.324555320336761, -5.692099788303083,
                                                                -0.8973665961010275),
                                                volmdlr.Point3D(3, 3, 2))
        circle2 = curves.Circle3D.from_3_points(
            volmdlr.Point3D(1.2067665579541171, -1.246879774203074, -0.4359328108960321),
            volmdlr.Point3D(-1.2905737351068276, -5.961765089244547, -0.9872550297481824),
            volmdlr.Point3D(7.33260591629263, -4.272128323147327, -0.4240427743824422))
        intersections = circle1.circle_intersections(circle2)
        expected_results = [volmdlr.Point3D(-1.2905737351057338, -5.961765089245487, -0.9872550297484957),
                            volmdlr.Point3D(7.332605916292026, -4.272128323148586, -0.42404277438286175)]
        for intersection, expected_result in zip(intersections, expected_results):
            self.assertTrue(intersection.is_close(expected_result))

        circle1 = curves.Circle3D(volmdlr.Frame3D(origin=volmdlr.Point3D(1.0, 1.0, -0.9595959595959596),
                                                  u=volmdlr.Vector3D(1.0, 0.0, 0.0),
                                                  v=volmdlr.Vector3D(0.0, 1.0, 0.0),
                                                  w=volmdlr.Vector3D(0.0, 0.0, 1.0)), 1)
        circle2 = curves.Circle3D(volmdlr.Frame3D(
            origin=volmdlr.Point3D(-3.3306690738754696e-16, 0.9999999999986168, -0.959595959596),
            u=volmdlr.Vector3D(1.0, 0.0, 0.0),
            v=volmdlr.Vector3D(0.0, 1.0, 0.0),
            w=volmdlr.Vector3D(0.0, 0.0, 1.0)), 1.3977300414963283)
        circle_intersections = circle1.circle_intersections(circle2)
        self.assertTrue(circle_intersections[0].is_close(
            volmdlr.Point3D(0.976824634449, 1.999731415147, -0.959595959596)))
        self.assertTrue(circle_intersections[1].is_close(
            volmdlr.Point3D(0.976824634452, 0.000268584853, -0.959595959596)))

    def test_point_distance(self):
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)

        circle = curves.Circle3D(frame, 2)
        point = volmdlr.Point3D(3, -2.5, 2)

        point_distance = circle.point_distance(point)
        self.assertAlmostEqual(point_distance, 3.341699272287294)


if __name__ == '__main__':
    unittest.main()
