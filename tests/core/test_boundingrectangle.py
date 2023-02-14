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

    def test_center(self):
        self.assertEqual(self.br.center(), volmdlr.Point2D(2.5, 3.0))
        self.assertNotEqual(self.br.center(), volmdlr.O2D)

    def test_b_rectangle_intersection(self):
        br2 = BoundingRectangle(-1.0, 2.0, -1.0, 2.0)
        self.assertTrue(self.br.b_rectangle_intersection(br2))
        br3 = BoundingRectangle(6.0, 7.0, 6.0, 7.0)
        self.assertFalse(self.br.b_rectangle_intersection(br3))

    def test_is_inside_b_rectangle(self):
        br2 = BoundingRectangle(-1.0, 6.0, -1.0, 6.0)
        self.assertFalse(self.br.is_inside_b_rectangle(br2))
        br3 = BoundingRectangle(1.0, 4.0, 1.0, 4.0)
        self.assertTrue(br3.is_inside_b_rectangle(self.br))
        br4 = BoundingRectangle(6.0, 7.0, 6.0, 7.0)
        self.assertFalse(self.br.is_inside_b_rectangle(br4))

    def test_point_belongs(self):
        self.assertTrue(self.br.point_belongs(volmdlr.Point2D(3.0, 4.0)))
        self.assertFalse(self.br.point_belongs(volmdlr.Point2D(6.0, 7.0)))

    def test_intersection_area(self):
        br2 = BoundingRectangle(1.0, 4.0, 1.0, 4.0)
        self.assertEqual(self.br.intersection_area(br2), 9.0)
        br3 = BoundingRectangle(6.0, 7.0, 6.0, 7.0)
        self.assertEqual(self.br.intersection_area(br3), 0.0)
        br4 = BoundingRectangle(-1.0, 6.0, -1.0, 6.0)
        self.assertEqual(self.br.intersection_area(br4), 35.0)

    def test_distance_to_b_rectangle(self):
        br2 = BoundingRectangle(1.0, 4.0, 1.0, 4.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br2), 0.0)
        br3 = BoundingRectangle(-4.0, -1.0, -5.0, -2.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br3), 1.0)
        br4 = BoundingRectangle(7.0, 10.0, 9.0, 10.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br4), 5**0.5)
        br5 = BoundingRectangle(-4.0, -1.0, 9.0, 10.0)
        self.assertEqual(self.br.distance_to_b_rectangle(br5), 2**0.5)

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


if __name__ == "__main__":
    unittest.main()
