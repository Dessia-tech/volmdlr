"""
Unittets for edges.Fullarc2D.

"""
import unittest
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import edges, curves


class TestFullArc2D(unittest.TestCase):
    circle2d = curves.Circle2D(volmdlr.OXY, 1)
    fullarc2d = edges.FullArc2D(circle2d, volmdlr.Point2D(-1, 0))

    def test_trim(self):
        split_point1 = volmdlr.Point2D(-0.7071067811865475, 0.7071067811865475)
        split_point2 = volmdlr.Point2D(-0.7071067811865475, -0.7071067811865475)
        split = self.fullarc2d.trim(split_point1, split_point2)
        self.assertEqual(split, edges.Arc2D(self.circle2d, volmdlr.Point2D(-0.7071067811865475, 0.7071067811865475),
                                            volmdlr.Point2D(-0.7071067811865475, -0.7071067811865475)))


if __name__ == '__main__':
    unittest.main()
