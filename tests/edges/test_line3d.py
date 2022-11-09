import unittest

import volmdlr
from volmdlr import edges


class TestLine3D(unittest.TestCase):
    def test_to_2d(self):
        line_3d = edges.Line3D(volmdlr.O3D, volmdlr.Point3D(1, 1, 1))
        line_2d = line_3d.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
        self.assertEqual(line_2d.point1, volmdlr.Point2D(0, 0))
        self.assertEqual(line_2d.point2, volmdlr.Point2D(1, 1))


if __name__ == '__main__':
    unittest.main()