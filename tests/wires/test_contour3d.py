import unittest

from volmdlr.wires import Contour3D


class TestContour3D(unittest.TestCase):

    def test_order_contour(self):
        contour_to_order = Contour3D.load_from_file('wires/contour_order.json')
        self.assertFalse(contour_to_order.is_ordered())
        contour_to_order.order_contour()
        self.assertTrue(contour_to_order.is_ordered())


if __name__ == '__main__':
    unittest.main()
