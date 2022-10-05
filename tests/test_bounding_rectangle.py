import unittest
import volmdlr as vm
from volmdlr.core import BoundingRectangle
import volmdlr.wires as vmw
import volmdlr.edges as vme

line_seg1 = vme.LineSegment2D(vm.Point2D(-0.5, -0.2), vm.O2D)
line_seg2 = vme.LineSegment2D(vm.O2D, vm.Point2D(0.3, 1))
line_seg3 = vme.LineSegment2D(vm.Point2D(0.3, 1), vm.Point2D(1, 1))
line_seg4 = vme.LineSegment2D(vm.Point2D(1, 1), vm.Point2D(1, -0.5))
line_seg5 = vme.LineSegment2D(vm.Point2D(1, -0.5), vm.Point2D(-0.5, -0.2))

contour1 = vmw.Contour2D([line_seg1, line_seg2, line_seg3, line_seg4, line_seg5])
bd_rectangle = contour1.bounding_rectangle()
xmin = bd_rectangle[0]
xmax = bd_rectangle[1]
ymin = bd_rectangle[2]
ymax = bd_rectangle[3]

b_rectangle1 = BoundingRectangle(xmin, xmax, ymin, ymax)

line_fig2_seg1 = vme.LineSegment2D(vm.Point2D(0, 1), vm.Point2D(0.25, 0.5))
line_fig2_seg2 = vme.LineSegment2D(vm.Point2D(0.25, 0.5), vm.Point2D(1, 0.5))
line_fig2_seg3 = vme.LineSegment2D(vm.Point2D(1, 0.5), vm.Point2D(0.45, 0))
line_fig2_seg4 = vme.LineSegment2D(vm.Point2D(0.45, 0), vm.Point2D(1, -1))
line_fig2_seg5 = vme.LineSegment2D(vm.Point2D(1, -1), vm.Point2D(0, -0.5))
line_fig2_seg6 = vme.LineSegment2D(vm.Point2D(0, -0.5), vm.Point2D(-1, -1))
line_fig2_seg7 = vme.LineSegment2D(vm.Point2D(-1, -1), vm.Point2D(-0.45, 0))
line_fig2_seg8 = vme.LineSegment2D(vm.Point2D(-0.45, 0), vm.Point2D(-1, 0.5))
line_fig2_seg9 = vme.LineSegment2D(vm.Point2D(-1, 0.5), vm.Point2D(-0.25, 0.5))
line_fig2_seg10 = vme.LineSegment2D(vm.Point2D(-0.25, 0.5), vm.Point2D(0, 1))

contour2 = vmw.Contour2D([line_fig2_seg1, line_fig2_seg2, line_fig2_seg3, line_fig2_seg4,
                          line_fig2_seg5, line_fig2_seg6, line_fig2_seg7, line_fig2_seg8, line_fig2_seg9,
                          line_fig2_seg10])

bd_rectangle = contour2.bounding_rectangle()

xmin = bd_rectangle[0]
xmax = bd_rectangle[1]
ymin = bd_rectangle[2]
ymax = bd_rectangle[3]
b_rectangle2 = BoundingRectangle(xmin, xmax, ymin, ymax)

line_fig3_seg1 = vme.LineSegment2D(vm.Point2D(5, 1), vm.Point2D(5.25, 0.5))
line_fig3_seg2 = vme.LineSegment2D(vm.Point2D(5.25, 0.5), vm.Point2D(6, 0.5))
line_fig3_seg3 = vme.LineSegment2D(vm.Point2D(6, 0.5), vm.Point2D(5.45, 0))
line_fig3_seg4 = vme.LineSegment2D(vm.Point2D(5.45, 0), vm.Point2D(6, -1))
line_fig3_seg5 = vme.LineSegment2D(vm.Point2D(6, -1), vm.Point2D(5, -0.5))
line_fig3_seg6 = vme.LineSegment2D(vm.Point2D(5, -0.5), vm.Point2D(4, -1))
line_fig3_seg7 = vme.LineSegment2D(vm.Point2D(4, -1), vm.Point2D(3.55, 0))
line_fig3_seg8 = vme.LineSegment2D(vm.Point2D(3.55, 0), vm.Point2D(4, 0.5))
line_fig3_seg9 = vme.LineSegment2D(vm.Point2D(4, 0.5), vm.Point2D(3.75, 0.5))
line_fig3_seg10 = vme.LineSegment2D(vm.Point2D(3.75, 0.5), vm.Point2D(5, 1))

contour3 = vmw.Contour2D([line_fig3_seg1, line_fig3_seg2, line_fig3_seg3, line_fig3_seg4,
                          line_fig3_seg5, line_fig3_seg6, line_fig3_seg7, line_fig3_seg8, line_fig3_seg9,
                          line_fig3_seg10])
bd_rectangle = contour3.bounding_rectangle()

xmin = bd_rectangle[0]
xmax = bd_rectangle[1]
ymin = bd_rectangle[2]
ymax = bd_rectangle[3]
b_rectangle3 = BoundingRectangle(xmin, xmax, ymin, ymax)


class TestBoundingRectangle(unittest.TestCase):

    def test_area(self):
        self.assertEqual(b_rectangle1.area(), 2.25)  # add assertion here

    def test_center(self):
        self.assertEqual(b_rectangle1.center(), [0.25, 0.25])  # add assertion here

    def test_intersection_true(self):
        self.assertTrue(b_rectangle1.b_rectangle_intersection(b_rectangle2))

    def test_intersection_false(self):
        self.assertFalse(b_rectangle1.b_rectangle_intersection(b_rectangle3))

    def test_inside_true(self):
        self.assertTrue(b_rectangle1.is_inside_b_rectangle(b_rectangle2))

    def test_inside_false(self):
        self.assertFalse(b_rectangle1.is_inside_b_rectangle(b_rectangle3))

    def test_point_belongs_true(self):
        self.assertTrue(b_rectangle1.point_belongs([0.25, 0.25]))

    def test_point_belongs_false(self):
        self.assertFalse(b_rectangle1.point_belongs([1, 1]))

    def test_intersection_area(self):
        self.assertEqual(b_rectangle1.intersection_area(b_rectangle2), 2.25)

    def test_intersection_area_null(self):
        self.assertEqual(b_rectangle1.intersection_area(b_rectangle3), 0)

    def test_distance_to_point(self):
        self.assertEqual(b_rectangle1.distance_to_point([0.25, 0.25]), 0.75)

    def test_distance_to_point(self):
        self.assertEqual(b_rectangle1.distance_to_point([0.25, 2]), 1)

    def test_distance_to_b_rectangle_null(self):
        self.assertEqual(b_rectangle1.distance_to_b_rectangle(b_rectangle2), 0)

    def test_distance_to_b_rectangle(self):
        self.assertEqual(b_rectangle2.distance_to_b_rectangle(b_rectangle3), 3)


if __name__ == '__main__':
    unittest.main()
