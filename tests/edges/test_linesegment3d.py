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

    def test_matrix_distance(self):
        edge1 = volmdlr.edges.LineSegment3D(volmdlr.Point3D(1, 1, 1), volmdlr.Point3D(2, 4, 5))
        edge2 = volmdlr.edges.LineSegment3D(volmdlr.Point3D(1, -1, 1), volmdlr.Point3D(3, -2, 7))
        matrix_distance = edge1.matrix_distance(edge2)
        self.assertEqual(matrix_distance[0], volmdlr.Point3D(1, 1, 1))
        self.assertEqual(matrix_distance[1], volmdlr.Point3D(1, -1, 1))


if __name__ == '__main__':
    unittest.main()
