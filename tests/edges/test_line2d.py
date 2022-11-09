import unittest

import volmdlr
from volmdlr import edges


class TestLine2D(unittest.TestCase):

    def test_point_distance(self):
        line_2d = edges.Line2D(volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1))
        self.assertEqual(line_2d.point_distance(volmdlr.Point2D(2, 2)), 0.0)
        self.assertEqual(line_2d.point_distance(volmdlr.Point2D(1, 2)), 0.7071067811865475)


if __name__ == '__main__':
    unittest.main()