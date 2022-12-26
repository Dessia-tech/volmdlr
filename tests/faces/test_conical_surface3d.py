import unittest
import math

import volmdlr
from volmdlr import faces, edges


class TestConicalSurface3D(unittest.TestCase):
    conical_surface = faces.ConicalSurface3D(volmdlr.OXYZ, math.pi/3)

    def test_arc3d_to_2d(self):
        arc1 = edges.Arc3D(volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 1 / math.sqrt(3)),
                           volmdlr.Point3D(-1, 0, 1 / math.sqrt(3)),
                           volmdlr.Point3D(-1 / math.sqrt(2), -1 / math.sqrt(2), 1 / math.sqrt(3)))
        arc2 = edges.Arc3D(volmdlr.Point3D(0, -1, 1 / math.sqrt(3)),
                           volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 1 / math.sqrt(3)),
                           volmdlr.Point3D(1, 0, 1 / math.sqrt(3)))
        test1 = self.conical_surface.arc3d_to_2d(arc3d=arc1)[0]
        test2 = self.conical_surface.arc3d_to_2d(arc3d=arc2)[0]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(test1, edges.LineSegment2D)
        self.assertIsInstance(test2, edges.LineSegment2D)

        # Assert that the returned object is right on the parametric domain (take into account periodicity)
        self.assertEqual(test1.start, volmdlr.Point2D(0.75 * math.pi, 0.5773502691896258))
        self.assertEqual(test1.end, volmdlr.Point2D(1.25 * math.pi, 0.5773502691896258))
        self.assertEqual(test2.start, volmdlr.Point2D(-0.5 * math.pi, 0.5773502691896258))
        self.assertEqual(test2.end, volmdlr.Point2D(-2 * math.pi, 0.5773502691896258))


if __name__ == '__main__':
    unittest.main()
