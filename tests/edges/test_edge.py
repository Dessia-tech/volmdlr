import unittest

import volmdlr
from volmdlr import edges
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

    def test_split_between_two_points(self):
        point1 = volmdlr.Point3D(0.1744332430903422, 0.033444245563080795, 0.07798520478978595)
        point2 = volmdlr.Point3D(0.177922447446, 0.03351981629780013, 0.07827867754649165)
        arc3d = edges.Arc3D.load_from_file("edges/arc3d_split_between_two_points.json")
        new_arc3d = arc3d.split_between_two_points(point1, point2)
        self.assertTrue(new_arc3d)
        self.assertTrue(new_arc3d.start.is_close(point2))
        self.assertTrue(new_arc3d.end.is_close(point1))


if __name__ == '__main__':
    unittest.main()
