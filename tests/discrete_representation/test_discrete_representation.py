"""
Unit testing of 'DiscreteRepresentation' class.
"""
import math
import unittest

from volmdlr.discrete_representation import DECIMALS, DiscreteRepresentation


class TestDiscreteRepresentation(unittest.TestCase):
    def test_check_center_is_in_implicit_grid(self):
        self.assertTrue(DiscreteRepresentation.check_center_is_in_implicit_grid((0.5, 0.5), 1))
        self.assertFalse(DiscreteRepresentation.check_center_is_in_implicit_grid((0.0, 0.0), 1))

    def test_check_element_size_number_of_decimals(self):
        self.assertWarns(
            Warning, DiscreteRepresentation._check_element_size_number_of_decimals, round(math.pi, DECIMALS)
        )
