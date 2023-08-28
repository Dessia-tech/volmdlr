"""
Unit testing for pixelization.
"""

import unittest

from volmdlr import Point2D, Vector2D
from volmdlr.discrete_representation import MatrixBasedPixelization, PointBasedPixelization
from volmdlr.edges import LineSegment2D
from volmdlr.wires import ClosedPolygon2D


class TestPixelization(unittest.TestCase):
    def setUp(self) -> None:
        points = [
            Point2D(5, 1),
            Point2D(5.25, 0.5),
            Point2D(6, 0.5),
            Point2D(5.45, 0),
            Point2D(6, -1),
            Point2D(5, -0.5),
            Point2D(4, -1),
            Point2D(4.55, 0),
            Point2D(4, 0.5),
            Point2D(4.75, 0.5),
            Point2D(5, 1),
        ]

        self.line_segment = LineSegment2D(points[0], points[1])
        self.closed_polygon_1 = ClosedPolygon2D(points)
        self.closed_polygon_2 = self.closed_polygon_1.translation(Vector2D(1.2, -0.5))

        self.pixel_size = 0.05

    def test_from_line_segment(self):
        pixelization_1 = PointBasedPixelization.from_line_segment(self.line_segment, self.pixel_size)
        pixelization_2 = MatrixBasedPixelization.from_line_segment(self.line_segment, self.pixel_size)

        self.assertEqual(24, len(pixelization_1))
        self.assertEqual(24, len(pixelization_2))

    def test_from_polygon(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)

        self.assertEqual(284, len(pixelization_1))
        self.assertEqual(284, len(pixelization_2))

    def test_to_point_based_pixelization(self):
        pixelization_1 = PointBasedPixelization.from_line_segment(self.line_segment, self.pixel_size)
        pixelization_2 = MatrixBasedPixelization.from_line_segment(
            self.line_segment, self.pixel_size
        ).to_point_based_pixelization()

        self.assertEqual(pixelization_1, pixelization_2)

    def test_to_point_matrix_based_pixelization(self):
        pixelization_1 = PointBasedPixelization.from_line_segment(
            self.line_segment, self.pixel_size
        ).to_matrix_based_pixelization()
        pixelization_2 = MatrixBasedPixelization.from_line_segment(self.line_segment, self.pixel_size)

        self.assertEqual(pixelization_1, pixelization_2)

    def test_union(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        pixelization_3 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_4 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        union_1 = pixelization_1.union(pixelization_2)
        union_2 = pixelization_3.union(pixelization_4)

        self.assertEqual(union_1, union_2.to_point_based_pixelization())
        self.assertEqual(union_2, union_1.to_matrix_based_pixelization())
        self.assertEqual(554, len(union_1))

    def test_difference(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        pixelization_3 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_4 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        difference_1 = pixelization_1.difference(pixelization_2)
        difference_2 = pixelization_3.difference(pixelization_4)

        self.assertEqual(difference_1, difference_2.to_point_based_pixelization())
        self.assertEqual(difference_2, difference_1.to_matrix_based_pixelization())
        self.assertEqual(270, len(difference_1))

    def test_intersection(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        pixelization_3 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_4 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        intersection_1 = pixelization_1.intersection(pixelization_2)
        intersection_2 = pixelization_3.intersection(pixelization_4)

        self.assertEqual(intersection_1, intersection_2.to_point_based_pixelization())
        self.assertEqual(intersection_2, intersection_1.to_matrix_based_pixelization())
        self.assertEqual(14, len(intersection_1))

    def test_symmetric_difference(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        pixelization_3 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_4 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_2, self.pixel_size)

        symmetric_difference_1 = pixelization_1.symmetric_difference(pixelization_2)
        symmetric_difference_2 = pixelization_3.symmetric_difference(pixelization_4)

        self.assertEqual(symmetric_difference_1, symmetric_difference_2.to_point_based_pixelization())
        self.assertEqual(symmetric_difference_2, symmetric_difference_1.to_matrix_based_pixelization())
        self.assertEqual(540, len(symmetric_difference_1))

    def test_inverse(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)

        inverse_1 = pixelization_1.inverse()
        inverse_2 = pixelization_2.inverse()

        self.assertEqual(inverse_1, inverse_2.to_point_based_pixelization())
        self.assertEqual(inverse_2, inverse_1.to_matrix_based_pixelization())
        self.assertEqual(1480, len(inverse_1))

    def test_serialization(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)

        self.assertEqual(pixelization_1, PointBasedPixelization.dict_to_object(pixelization_1.to_dict()))
        self.assertEqual(pixelization_2, MatrixBasedPixelization.dict_to_object(pixelization_2.to_dict()))

    def test_min_grid_center(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)

        self.assertEqual((3.975, -1.025), pixelization_1.min_grid_center)
        self.assertEqual((3.975, -1.025), pixelization_2.min_grid_center)

    def test_max_grid_center(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)
        pixelization_2 = MatrixBasedPixelization.from_closed_polygon(self.closed_polygon_1, self.pixel_size)

        self.assertEqual((6.025, 1.025), pixelization_1.max_grid_center)
        self.assertEqual((6.025, 1.025), pixelization_2.max_grid_center)

    def test_fill_outer_pixels(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(
            self.closed_polygon_1, self.pixel_size
        ).fill_outer_pixels()
        pixelization_2 = MatrixBasedPixelization.from_closed_polygon(
            self.closed_polygon_1, self.pixel_size
        ).fill_outer_pixels()

        self.assertEqual(pixelization_1, pixelization_2.to_point_based_pixelization())
        self.assertEqual(pixelization_2, pixelization_1.to_matrix_based_pixelization())
        self.assertEqual(1180, len(pixelization_1))

    def test_fill_enclosed_pixels(self):
        pixelization_1 = PointBasedPixelization.from_closed_polygon(
            self.closed_polygon_1, self.pixel_size
        ).fill_enclosed_pixels()
        pixelization_2 = MatrixBasedPixelization.from_closed_polygon(
            self.closed_polygon_1, self.pixel_size
        ).fill_enclosed_pixels()

        self.assertEqual(pixelization_1, pixelization_2.to_point_based_pixelization())
        self.assertEqual(pixelization_2, pixelization_1.to_matrix_based_pixelization())
        self.assertEqual(868, len(pixelization_1))
