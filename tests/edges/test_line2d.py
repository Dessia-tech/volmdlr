import unittest

import volmdlr
from volmdlr import edges


class TestLine2D(unittest.TestCase):

    def test_point_distance(self):
        line_2d = edges.Line2D(volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1))
        self.assertEqual(line_2d.point_distance(volmdlr.Point2D(2, 2)), 0.0)
        self.assertEqual(line_2d.point_distance(volmdlr.Point2D(1, 2)), 0.7071067811865475)

    def test_point_belongs(self):
        line = edges.Line2D(volmdlr.O2D, volmdlr.Point2D(1, 1))
        point1 = volmdlr.Point2D(2, 2)
        point2 = volmdlr.Point2D(1, 2)
        self.assertTrue(line.point_belongs(point1))
        self.assertFalse(line.point_belongs(point2))


if __name__ == '__main__':
    unittest.main()