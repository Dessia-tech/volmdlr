"""
Unit tests for volmdlr.wires.Contour
"""
import unittest

import volmdlr.wires


class TestContour(unittest.TestCase):
    def test_is_overlapping(self):
        p1 = [
            volmdlr.Point2D(-0.2, -0.2),
            volmdlr.Point2D(0.2, -0.2),
            volmdlr.Point2D(0.2, 0.2),
            volmdlr.Point2D(-0.2, 0.2),
        ]

        contour1 = volmdlr.wires.ClosedPolygon2D(p1)
        contour2 = contour1.translation(volmdlr.Vector2D(0.5, 0))
        contour3 = contour1.translation(volmdlr.Vector2D(0.1, 0))

        self.assertFalse(contour1.is_overlapping(contour2))
        self.assertTrue(contour1.is_overlapping(contour3))


if __name__ == "__main__":
    unittest.main()
