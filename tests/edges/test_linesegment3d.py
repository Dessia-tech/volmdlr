import unittest

import volmdlr
from volmdlr import edges


class TestLineSegment3D(unittest.TestCase):
    linesegment1 = edges.LineSegment3D(volmdlr.Point3D(-0.2, -0.2, 0.2),
                                       volmdlr.Point3D(-0.2, -0.2, -0.2), name='linesegment1')
    linesegment2 = edges.LineSegment3D(volmdlr.Point3D(-0.225, -0.2, 0.125),
                                       volmdlr.Point3D(-0.225, -0.2, 0.275), name='linesegment2')
    linesegment3 = edges.LineSegment3D(volmdlr.Point3D(-0.225, -0.2, 0.275),
                                       volmdlr.Point3D(0.225, -0.2, 0.275), name='linesegment3')
    linesegment4 = edges.LineSegment3D(volmdlr.Point3D(0.225, -0.2, 0.275),
                                       volmdlr.Point3D(0.225, -0.2, 0.125), name='linesegment4')
    linesegment5 = edges.LineSegment3D(volmdlr.Point3D(0.225, -0.2, 0.125),
                                       volmdlr.Point3D(-0.225, -0.2, 0.125), name='linesegment5')

    def test_line_intersection(self):
        self.assertFalse(self.linesegment1.line_intersections(self.linesegment3.line))
        self.assertEqual(self.linesegment3.line_intersections(self.linesegment1.line),
                         [volmdlr.Point3D(-0.2, -0.2, 0.275)])

    def test_linesegment_intersection(self):
        for lineseg in [self.linesegment2, self.linesegment3, self.linesegment4]:
            self.assertFalse(self.linesegment1.linesegment_intersections(lineseg))
        self.assertTrue(self.linesegment1.linesegment_intersections(self.linesegment5)[0].is_close(
            volmdlr.Point3D(-0.2, -0.2, 0.125)))

    def test_matrix_distance(self):
        edge1 = volmdlr.edges.LineSegment3D(volmdlr.Point3D(1, 1, 1), volmdlr.Point3D(2, 4, 5))
        edge2 = volmdlr.edges.LineSegment3D(volmdlr.Point3D(1, -1, 1), volmdlr.Point3D(3, -2, 7))
        matrix_distance = edge1.matrix_distance(edge2)
        self.assertEqual(matrix_distance[0], volmdlr.Point3D(1, 1, 1))
        self.assertEqual(matrix_distance[1], volmdlr.Point3D(1, -1, 1))

    def test_minimum_distance(self):
        lineseg1 = edges.LineSegment3D(
            volmdlr.Point3D(-1.405755599081579, -5.5244351616677, -0.9525260021364759),
            volmdlr.Point3D(-1.2354684565606084, -6.038042700761348, -0.991655683450264))
        lineseg2 = edges.LineSegment3D(
            volmdlr.Point3D(-1.9316894719468425, -5.304938228054382, -0.7683127426847942),
            volmdlr.Point3D(0.15774316590678739, -6.857594985919942, -1.2858649953066474))
        dist, min_dist_point1, min_dist_point2 = lineseg1.minimum_distance(lineseg2, True)
        self.assertAlmostEqual(dist, 0.050974695310070525)
        self.assertTrue(min_dist_point1.is_close(
            volmdlr.Point3D(-1.3261157928570702, -5.764638880417798, -0.9708261520557401)))
        self.assertTrue(min_dist_point2.is_close(
            volmdlr.Point3D(-1.3135876697596944, -5.764249489316487, -0.9214164964388291)))


if __name__ == '__main__':
    unittest.main()
