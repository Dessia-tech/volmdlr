import unittest
from volmdlr.models.edges import bspline1, lineseg, arc, arc_ellipse2d


class TestEdge(unittest.TestCase):
    def test_direction_independent_is_close(self):
        self.assertFalse(bspline1.direction_independent_is_close(arc))
        self.assertTrue(bspline1.direction_independent_is_close(bspline1))
        self.assertTrue(bspline1.direction_independent_is_close(bspline1.reverse()))
        self.assertTrue(lineseg.direction_independent_is_close(lineseg))
        self.assertTrue(lineseg.direction_independent_is_close(lineseg.reverse()))
        self.assertTrue(arc.direction_independent_is_close(arc))
        self.assertTrue(arc.direction_independent_is_close(arc.reverse()))
        self.assertFalse(arc.direction_independent_is_close(arc.complementary()))
        self.assertTrue(arc_ellipse2d.direction_independent_is_close(arc_ellipse2d))
        self.assertTrue(arc_ellipse2d.direction_independent_is_close(arc_ellipse2d.reverse()))
        self.assertFalse(arc_ellipse2d.direction_independent_is_close(arc_ellipse2d.complementary()))


if __name__ == '__main__':
    unittest.main()
