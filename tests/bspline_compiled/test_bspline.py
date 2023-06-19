import unittest
import volmdlr.bspline_compiled as bc


class TestBspline(unittest.TestCase):
    def test_binomial_coefficient(self):
        self.assertEqual(bc.binomial_coefficient(4, 2), 6)


if __name__ == '__main__':
    unittest.main()
