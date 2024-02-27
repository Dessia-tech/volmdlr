import math
import os
import unittest
from itertools import product

import volmdlr
from volmdlr import edges, curves


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'arcellipse_objects')


class TestArcEllipse2D(unittest.TestCase):
    u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
    v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
    ellipse2d = curves.Ellipse2D(2, 1, volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector))
    arc_ellipse2d = edges.ArcEllipse2D(ellipse2d, start=volmdlr.Point2D(0.5, 1.5), end=volmdlr.Point2D(1.5, 0.5))
    discretized_points = arc_ellipse2d.discretization_points(number_points=5)

    def test_init(self):
        list_points = self.ellipse2d.discretization_points(number_points=9)
        list_lengths = []
        for point1, point2 in product(list_points, repeat=2):
            if not point1.is_close(point2):
                arc_ellipse2d = edges.ArcEllipse2D(self.ellipse2d, start=point1, end=point2)
                list_lengths.append(round(arc_ellipse2d.length(), 7))
        list_expected_lengths = [0.9656637, 2.4221121, 3.8785604, 4.8442241, 5.8098879, 7.2663362, 8.7227845,
                                 8.7227845, 1.4564483, 2.9128966, 3.8785604, 4.8442241, 6.3006724, 7.7571207,
                                 8.7227845, 7.2663362, 8.2319999, 1.4564483, 2.4221121, 3.3877758, 4.8442241,
                                 6.3006724, 7.2663362, 5.8098879, 6.7755516, 8.2319999, 0.9656637, 1.9313275,
                                 3.3877758, 4.8442241, 5.8098879, 4.8442241, 5.8098879, 7.2663362, 8.7227845,
                                 0.9656637, 2.4221121, 3.8785604, 4.8442241, 3.8785604, 4.8442241, 6.3006724,
                                 7.7571207, 8.7227845, 1.4564483, 2.9128966, 3.8785604, 2.4221121, 3.3877758,
                                 4.8442241, 6.3006724, 7.2663362, 8.2319999, 1.4564483, 2.4221121, 0.9656637,
                                 1.9313275, 3.3877758, 4.8442241, 5.8098879, 6.7755516, 8.2319999, 0.9656637,
                                 0.9656637, 2.4221121, 3.8785604, 4.8442241, 5.8098879, 7.2663362, 8.7227845]
        for length, expected_length in zip(list_lengths, list_expected_lengths):
            self.assertAlmostEqual(length, expected_length)

    def test_length(self):
        self.assertAlmostEqual(self.arc_ellipse2d.length(), 7.757120732103266)

    def test_to_3d(self):
        expected_points = [volmdlr.Point3D(1.0, 1.5, 2.5),
                           volmdlr.Point3D(1.0, 0.0803671345838155, 1.477169381251236),
                           volmdlr.Point3D(1.0, -0.581093271630204, 0.060946250290002224),
                           volmdlr.Point3D(1.0, 0.06094625029000189, -0.5810932716302042),
                           volmdlr.Point3D(1.0, 1.4771693812512356, 0.08036713458381528),
                           volmdlr.Point3D(1.0, 2.5, 1.5)]
        arcellipse3d = self.arc_ellipse2d.to_3d(volmdlr.Point3D(1, 1, 1), volmdlr.Y3D, volmdlr.Z3D)
        points = arcellipse3d.discretization_points(number_points=6)
        for point, expected_point in zip(points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_discretization_points(self):
        expected_discretized_points = [volmdlr.Point2D(0.5000000000000001, 1.5),
                                       volmdlr.Point2D(-1.1944775825843852, 0.11208538229199128),
                                       volmdlr.Point2D(-1.414213562373095, -1.414213562373095),
                                       volmdlr.Point2D(0.11208538229199028, -1.1944775825843859),
                                       volmdlr.Point2D(1.5, 0.49999999999999967)]
        for expected_point, point in zip(expected_discretized_points, self.discretized_points):
            self.assertTrue(expected_point.is_close(point))

    def test_point_belongs(self):
        self.assertTrue(self.arc_ellipse2d.point_belongs(self.discretized_points[1]))
        self.assertFalse(self.arc_ellipse2d.point_belongs(volmdlr.Point2D(1.4142135623730954, 1.4142135623730947)))

    def test_abscissa(self):
        self.assertAlmostEqual(self.arc_ellipse2d.abscissa(self.discretized_points[2]), 3.878560366051633)

    def test_reverse(self):
        reversed_arcellipse2d = self.arc_ellipse2d.reverse()
        self.assertEqual(reversed_arcellipse2d.start, self.arc_ellipse2d.end)
        self.assertEqual(reversed_arcellipse2d.end, self.arc_ellipse2d.start)
        reversed_arcellipse2d_length = reversed_arcellipse2d.length()
        self.assertAlmostEqual(reversed_arcellipse2d_length, self.arc_ellipse2d.length())
        self.assertEqual(reversed_arcellipse2d_length - reversed_arcellipse2d.abscissa(self.discretized_points[1]),
                         self.arc_ellipse2d.abscissa(self.discretized_points[1]))

    def test_bounding_rectangle(self):
        expected_bounds = [-1.5745077023090706, 1.5, -1.5745077023090708, 1.5]
        for i, bound in enumerate(self.arc_ellipse2d.bounding_rectangle.bounds()):
            self.assertAlmostEqual(bound, expected_bounds[i])

    def test_straight_line_area(self):
        self.assertAlmostEqual(self.arc_ellipse2d.straight_line_area(), 5.71238898038469)

    def test_linesegment_intersections(self):
        lineseg = edges.LineSegment2D(volmdlr.Point2D(2, -1), volmdlr.Point2D(-2, 2))
        inters = self.arc_ellipse2d.linesegment_intersections(lineseg)
        self.assertEqual(len(inters), 2)
        self.assertTrue(inters[1].is_close(volmdlr.Point2D(-0.5154201866080642, 0.8865651399560482)))
        self.assertTrue(inters[0].is_close(volmdlr.Point2D(1.0636435368618713, -0.29773265264640336)))

    def test_translation(self):
        translated_ellipse = self.arc_ellipse2d.translation(volmdlr.Vector2D(1, 1))
        self.assertTrue(translated_ellipse.start.is_close(volmdlr.Point2D(1.5, 2.5)))
        self.assertTrue(translated_ellipse.end.is_close(volmdlr.Point2D(2.5, 1.5)))

    def test_point_distance(self):
        point = volmdlr.Point2D(0, 0)
        point_distance = self.arc_ellipse2d.point_distance(point)
        self.assertAlmostEqual(point_distance, 1, 3)

    def test_split(self):
        middle_point = self.arc_ellipse2d.point_at_abscissa(self.arc_ellipse2d.length()*.5)
        split1 = self.arc_ellipse2d.split(middle_point)
        split2 = self.arc_ellipse2d.split(self.arc_ellipse2d.start)
        split3 = self.arc_ellipse2d.split(self.arc_ellipse2d.end)
        self.assertTrue(split1[0].start.is_close(self.arc_ellipse2d.start))
        self.assertTrue(split1[0].end.is_close(middle_point))
        self.assertTrue(split1[1].start.is_close(middle_point))
        self.assertTrue(split1[1].end.is_close(self.arc_ellipse2d.end))
        self.assertIsNone(split2[0])
        self.assertEqual(split2[1], self.arc_ellipse2d)
        self.assertIsNone(split3[1])
        self.assertEqual(split3[0], self.arc_ellipse2d)

    def test_point_at_abscissa(self):
        list_abscissas = [0, self.arc_ellipse2d.length() / 4, self.arc_ellipse2d.length() / 3,
                          self.arc_ellipse2d.length() / 2, self.arc_ellipse2d.length() * 0.75,
                          self.arc_ellipse2d.length()]
        list_points = []
        for abscissa in list_abscissas:
            point_at_abscissa = self.arc_ellipse2d.point_at_abscissa(abscissa)
            list_points.append(point_at_abscissa)
        expected_points = [volmdlr.Point2D(0.5, 1.5),
                           volmdlr.Point2D(-1.0268603272197745, 0.34573464229361694),
                           volmdlr.Point2D(-1.373935276754997, -0.19836341255508938),
                           volmdlr.Point2D(-1.414213568389262, -1.4142135563569278),
                           volmdlr.Point2D(0.3457346310451254, -1.026860335985359),
                           volmdlr.Point2D(1.5, 0.5)]
        for point, expected_point in zip(list_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

        arcellipse = edges.ArcEllipse2D.from_json(os.path.join(folder, "ellipse2d_point_at_abscissa.json"))
        abscissa = 1.4594044224379008
        point = arcellipse.point_at_abscissa(abscissa)
        self.assertTrue(point.is_close(volmdlr.Point2D(0.2409700344646443, -0.10841585996396141)))

    def test_complementary(self):
        complementary = self.arc_ellipse2d.complementary()
        interior = complementary.point_at_abscissa(complementary.length()*0.5)
        self.assertTrue(interior.is_close(volmdlr.Point2D(1.4142135565836316, 1.4142135681625583)))
        self.assertAlmostEqual(complementary.length(), 1.9313274884444551)

    def test_normal_vector(self):
        direc_vector = self.arc_ellipse2d.normal_vector(self.arc_ellipse2d.length() * 0.4)
        self.assertTrue(direc_vector.is_close(volmdlr.Vector2D(-1.1622013637446964, 0.25201219862839874)))

    def test_direction_vector(self):
        direc_vector = self.arc_ellipse2d.direction_vector(self.arc_ellipse2d.length() * 0.4)
        self.assertTrue(direc_vector.is_close(volmdlr.Vector2D(0.25201219862839874, 1.1622013637446964)))

    def test_rotation(self):
        rotationed_ellipse2d = self.arc_ellipse2d.rotation(volmdlr.O2D, math.pi / 4)
        self.assertTrue(rotationed_ellipse2d.start.is_close(
            volmdlr.Point2D(-0.7071067811865475, 1.4142135623730951)))
        self.assertTrue(rotationed_ellipse2d.end.is_close(
            volmdlr.Point2D(0.7071067811865477, 1.414213562373095)))

    def test_frame_mapping(self):
        frame_mapped_arcelipsse1 = self.arc_ellipse2d.frame_mapping(self.ellipse2d.frame, 'new')
        self.assertTrue(frame_mapped_arcelipsse1.start.is_close(
            volmdlr.Point2D(1.4142135623730951, 0.7071067811865477)))
        self.assertTrue(frame_mapped_arcelipsse1.end.is_close(
            volmdlr.Point2D(1.4142135623730951, -0.7071067811865477)))

    def test_get_shared_section(self):
        #test1
        u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
        v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
        ellipse2d = curves.Ellipse2D(2, 1, volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector))

        arc_ellipse2d = edges.ArcEllipse2D(ellipse2d, start=ellipse2d.point_at_abscissa(0.25 * ellipse2d.length()),
                                           end=ellipse2d.point_at_abscissa(0.75 * ellipse2d.length()))

        arc_ellipse2d_2 = edges.ArcEllipse2D(ellipse2d, start=ellipse2d.point_at_abscissa(0.6 * ellipse2d.length()),
                                             end=ellipse2d.point_at_abscissa(0.9 * ellipse2d.length()))
        get_shared_section = arc_ellipse2d.get_shared_section(arc_ellipse2d_2)
        expected_ellipse = edges.ArcEllipse2D(ellipse2d, volmdlr.Point2D(-0.4969829723203407, -1.4989916288867593),
                                              volmdlr.Point2D(0.7071067817853691, -0.7071067805877258))
        self.assertTrue(get_shared_section[0].is_close(expected_ellipse))

        # test2
        arc_ellipse2d_2 = edges.ArcEllipse2D(ellipse2d, start=ellipse2d.point_at_abscissa(0.4 * ellipse2d.length()),
                                             end=ellipse2d.point_at_abscissa(0.6 * ellipse2d.length()))

        get_shared_section = arc_ellipse2d.get_shared_section(arc_ellipse2d_2)
        self.assertEqual(get_shared_section[0], arc_ellipse2d_2)

    def test_delete_shared_section(self):
        u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
        v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
        ellipse2d = curves.Ellipse2D(2, 1, volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector))

        arc_ellipse2d = edges.ArcEllipse2D(ellipse2d, start=ellipse2d.point_at_abscissa(0.25 * ellipse2d.length()),
                                           end=ellipse2d.point_at_abscissa(0.75 * ellipse2d.length()))
        arc_ellipse2d_2 = edges.ArcEllipse2D(ellipse2d, start=ellipse2d.point_at_abscissa(0.4 * ellipse2d.length()),
                                             end=ellipse2d.point_at_abscissa(0.6 * ellipse2d.length()))

        delete_shared_section = arc_ellipse2d.delete_shared_section(arc_ellipse2d_2)
        self.assertEqual(len(delete_shared_section), 2)
        self.assertEqual(delete_shared_section[0].start, arc_ellipse2d.start)
        self.assertEqual(delete_shared_section[0].end, arc_ellipse2d_2.start)
        self.assertEqual(delete_shared_section[1].end, arc_ellipse2d.end)
        self.assertEqual(delete_shared_section[1].start, arc_ellipse2d_2.end)

    def test_straight_line_center_of_mass(self):
        straight_line_center_of_mass = self.arc_ellipse2d.straight_line_center_of_mass()
        self.assertTrue(straight_line_center_of_mass, volmdlr.Point2D(-0.4100897136188068, -0.4100897136188069))

if __name__ == '__main__':
    unittest.main()
