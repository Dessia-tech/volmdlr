import unittest

from volmdlr.wires import Contour3D


class TestContour3D(unittest.TestCase):

    def test_order_contour(self):
        contour_to_order = Contour3D.load_from_file('wires/contour_order.json')
        self.assertFalse(contour_to_order.is_ordered())
        contour_to_order.order_contour()
        self.assertTrue(contour_to_order.is_ordered())

    def test_merge_with(self):
        contour1_to_merge = Contour3D.load_from_file('wires/contour3d_merge_with1.json')
        contour2_to_merge = Contour3D.load_from_file('wires/contour3d_merge_with2.json')
        expected_contour1 = Contour3D.load_from_file('wires/expected_contour_merge_with1.json')
        expected_contour2 = Contour3D.load_from_file('wires/expected_contour_merge_with2.json')
        merged_contours = contour1_to_merge.merge_with(contour2_to_merge)
        self.assertEqual(merged_contours[0], expected_contour1)
        self.assertEqual(merged_contours[1], expected_contour2)


if __name__ == '__main__':
    unittest.main()
