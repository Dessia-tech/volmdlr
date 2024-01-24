import unittest
import volmdlr
from volmdlr.core import BoundingRectangle


class TestBoundingRectangle(unittest.TestCase):
    def setUp(self):
        self.xmin = 0.0
        self.xmax = 5.0
        self.ymin = -2.0
        self.ymax = 8.0
        self.br = BoundingRectangle(self.xmin, self.xmax, self.ymin, self.ymax)

        self.b_rectangle1 = BoundingRectangle(-0.5, 1.0, -0.5, 1.0)
        self.b_rectangle2 = BoundingRectangle(-1.0, 1.0, -1.0, 1.0)
        self.b_rectangle3 = BoundingRectangle(4.0, 6.0, -1.0, 1.0)

    def test_getitem(self):
        self.assertEqual(self.br[0], self.xmin)
        self.assertEqual(self.br[1], self.xmax)
        self.assertEqual(self.br[2], self.ymin)
        self.assertEqual(self.br[3], self.ymax)

    def test_bounds(self):
        self.assertTupleEqual(self.br.bounds(), (self.xmin, self.xmax, self.ymin, self.ymax))

    def test_plot(self):
        ax = self.br.plot()
        x = [self.xmin, self.xmax, self.xmax, self.xmin, self.xmin]
        y = [self.ymin, self.ymin, self.ymax, self.ymax, self.ymin]

        self.assertListEqual(ax.lines[0].get_xydata().tolist(), [list(p) for p in zip(x, y)])

    def test_area(self):
        self.assertEqual(self.br.area(), 50.0)
        self.assertNotEqual(self.br.area(), 0.0)

        self.assertEqual(self.b_rectangle1.area(), 2.25)

    def test_center(self):
        self.assertEqual(self.br.center(), volmdlr.Point2D(2.5, 3.0))
        self.assertNotEqual(self.br.center(), volmdlr.O2D)

        self.assertEqual(self.b_rectangle1.center(), volmdlr.Point2D(0.25, 0.25))

    def test_is_intersecting(self):
        br2 = BoundingRectangle(-1.0, 2.0, -1.0, 2.0)
        self.assertTrue(self.br.is_intersecting(br2))
        br3 = BoundingRectangle(6.0, 7.0, 6.0, 7.0)
        self.assertFalse(self.br.is_intersecting(br3))

        self.assertTrue(self.b_rectangle1.is_intersecting(self.b_rectangle2))
        self.assertFalse(self.b_rectangle1.is_intersecting(self.b_rectangle3))

    def test_is_inside_b_rectangle(self):
        br2 = BoundingRectangle(-1.0, 6.0, -1.0, 6.0)
        self.assertFalse(self.br.is_inside_b_rectangle(br2))
        br3 = BoundingRectangle(1.0, 4.0, 1.0, 4.0)
        self.assertTrue(br3.is_inside_b_rectangle(self.br))
        br4 = BoundingRectangle(6.0, 7.0, 6.0, 7.0)
        self.assertFalse(self.br.is_inside_b_rectangle(br4))

        self.assertTrue(self.b_rectangle1.is_inside_b_rectangle(self.b_rectangle2))
        self.assertFalse(self.b_rectangle1.is_inside_b_rectangle(self.b_rectangle3))

    def test_point_belongs(self):
        self.assertTrue(self.br.point_inside(volmdlr.Point2D(3.0, 4.0)))
        self.assertFalse(self.br.point_inside(volmdlr.Point2D(6.0, 7.0)))

        self.assertTrue(self.b_rectangle1.point_inside(volmdlr.Point2D(0.25, 0.25)))
        self.assertFalse(self.b_rectangle1.point_inside(volmdlr.Point2D(1.0, 1.0)))

    def test_intersection_area(self):
        br2 = BoundingRectangle(1.0, 4.0, 1.0, 4.0)
        self.assertEqual(self.br.intersection_area(br2), 9.0)
        br3 = BoundingRectangle(6.0, 7.0, 6.0, 7.0)
        self.assertEqual(self.br.intersection_area(br3), 0.0)
        br4 = BoundingRectangle(-1.0, 6.0, -1.0, 6.0)
        self.assertEqual(self.br.intersection_area(br4), 35.0)

        self.assertEqual(self.b_rectangle1.intersection_area(self.b_rectangle2), 2.25)
        self.assertEqual(self.b_rectangle1.intersection_area(self.b_rectangle3), 0)

    def test_distance_to_b_rectangle(self):
        br2 = BoundingRectangle(1.0, 4.0, 1.0, 4.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br2), 0.0)
        br3 = BoundingRectangle(-4.0, -1.0, -5.0, -2.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br3), 1.0)
        br4 = BoundingRectangle(7.0, 10.0, 9.0, 10.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br4), 5**0.5)
        br5 = BoundingRectangle(-4.0, -1.0, 9.0, 10.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br5), 2**0.5)

        self.assertEqual(self.b_rectangle1.distance_to_b_rectangle(self.b_rectangle2), 0)
        self.assertEqual(self.b_rectangle2.distance_to_b_rectangle(self.b_rectangle3), 3)

    def test_distance_to_point(self):
        p0 = volmdlr.O2D
        self.assertEqual(self.br.distance_to_point(p0), 0.0)
        p1 = volmdlr.Point2D(1.0, 1.0)
        self.assertEqual(self.br.distance_to_point(p1), 1.0)
        p3 = volmdlr.Point2D(-3.0, 0.0)
        self.assertEqual(self.br.distance_to_point(p3), 3.0)
        p4 = volmdlr.Point2D(2.0, -5.0)
        self.assertEqual(self.br.distance_to_point(p4), 3.0)
        p5 = volmdlr.Point2D(6.0, 9.0)
        self.assertEqual(self.br.distance_to_point(p5), 2**0.5)

        self.assertEqual(self.b_rectangle1.distance_to_point(volmdlr.Point2D(0.25, 0.25)), 0.75)
        self.assertEqual(self.b_rectangle1.distance_to_point(volmdlr.Point2D(0.25, 2)), 1)


if __name__ == "__main__":
    unittest.main()
