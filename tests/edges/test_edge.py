import unittest

import volmdlr
from volmdlr import edges
import volmdlr.step
from volmdlr.models.edges import bspline1, lineseg, arc, arc_ellipse2d, linesegment3d, arc3d, arc_ellipse3d


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

    def test_trim(self):
        point1 = volmdlr.Point3D(0.1744332430903422, 0.033444245563080795, 0.07798520478978595)
        point2 = volmdlr.Point3D(0.177922447446, 0.03351981629780013, 0.07827867754649165)
        arc3d = edges.Arc3D.from_json("edges/arc3d_split_between_two_points.json")
        new_arc3d = arc3d.trim(point1, point2)
        self.assertTrue(new_arc3d)
        self.assertTrue(new_arc3d.start.is_close(point2))
        self.assertTrue(new_arc3d.end.is_close(point1))

    def test_from_step(self):
        step = volmdlr.step.Step.from_file(filepath='edges/test_edge_from_step.stp')
        model = step.to_volume_model()
        self.assertTrue(model.primitives[0].faces[0].outer_contour3d.is_ordered())

    def test_to_step(self):
        current_id = 1
        content, current_id = linesegment3d().to_step(current_id)
        self.assertIn("LINE", content)
        self.assertIn("EDGE_CURVE", content)
        content, current_id = arc3d().to_step(current_id)
        self.assertIn("CIRCLE", content)
        self.assertIn("EDGE_CURVE", content)
        content, current_id = arc_ellipse3d().to_step(current_id)
        self.assertIn("ELLIPSE", content)
        self.assertIn("EDGE_CURVE", content)


if __name__ == '__main__':
    unittest.main()
