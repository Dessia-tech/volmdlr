import unittest

import volmdlr

from volmdlr.core import BoundingRectangle


class TestBoundingRectangle(unittest.TestCase):
    b_rectangle1 = BoundingRectangle(-0.5, 1.0, -0.5, 1)
    b_rectangle2 = BoundingRectangle(-1.0, 1.0, -1.0, 1.0)
    b_rectangle3 = BoundingRectangle(4.0, 6.0, -1.0, 1.0)

    def test_area(self):
        self.assertEqual(self.b_rectangle1.area(), 2.25)  # add assertion here

    def test_center(self):
        self.assertEqual(self.b_rectangle1.center(), volmdlr.Point2D(0.25, 0.25))  # add assertion here

    def test_intersection(self):
        self.assertTrue(self.b_rectangle1.b_rectangle_intersection(self.b_rectangle2))
        self.assertFalse(self.b_rectangle1.b_rectangle_intersection(self.b_rectangle3))

    def test_is_inside(self):
        self.assertTrue(self.b_rectangle1.is_inside_b_rectangle(self.b_rectangle2))
        self.assertFalse(self.b_rectangle1.is_inside_b_rectangle(self.b_rectangle3))

    def test_point_belongs(self):
        self.assertTrue(self.b_rectangle1.point_belongs(volmdlr.Point2D(0.25, 0.25)))
        self.assertFalse(self.b_rectangle1.point_belongs(volmdlr.Point2D(1, 1)))

    def test_intersection_area(self):
        self.assertEqual(self.b_rectangle1.intersection_area(self.b_rectangle2), 2.25)
        self.assertEqual(self.b_rectangle1.intersection_area(self.b_rectangle3), 0)

    def test_distance_to_point(self):
        self.assertEqual(self.b_rectangle1.distance_to_point(volmdlr.Point2D(0.25, 0.25)), 0.75)
        self.assertEqual(self.b_rectangle1.distance_to_point(volmdlr.Point2D(0.25, 2)), 1)

    def test_distance_to_b_rectangle(self):
        self.assertEqual(self.b_rectangle1.distance_to_b_rectangle(self.b_rectangle2), 0)
        self.assertEqual(self.b_rectangle2.distance_to_b_rectangle(self.b_rectangle3), 3)


if __name__ == '__main__':
    unittest.main()
