"""
Unittets for edges.Fullarc2D.

"""
import unittest
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import edges


class TestFullArc2D(unittest.TestCase):
    def test_split_between_two_points(self):
        fullarc, point1, point2 = DessiaObject.load_from_file(
            'edges/test_fullarc2d_split_between_two_points.json').primitives
        split_between_two_points = fullarc.split_between_two_points(point1, point2)
        self.assertEqual(split_between_two_points, edges.Arc2D(
            volmdlr.Point2D(0.5, 0.5499999999999999), volmdlr.Point2D(0.30000000000000004, 0.3499999999999999),
            volmdlr.Point2D(0.5, 0.1500000000000002), volmdlr.Point2D(0.49999999999999994, 0.3500000000000001)))


if __name__ == '__main__':
    unittest.main()
