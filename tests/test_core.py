#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing core module functions
"""

import math
import unittest

import numpy as npy

import volmdlr as vm
import volmdlr.wires


class TestContour(unittest.TestCase):

    def test_is_overlapping(self):

        v2D_1 = vm.Vector2D.random(-3, 3, -3, 3)

        # Testing if normalized vector has norm ==1 and is still colinear to original vector
        v2D_1_normalized = v2D_1.copy()
        v2D_1_normalized = v2D_1_normalized.unit_vector()
        self.assertAlmostEqual(v2D_1_normalized.norm(), 1)
        self.assertAlmostEqual(v2D_1_normalized.dot(v2D_1), v2D_1.norm())
        # Testing normal vector
        normal_v2D_1 = v2D_1.normal_vector()
        self.assertAlmostEqual(normal_v2D_1.dot(v2D_1), 0)
        normal_unit_v2D_1 = v2D_1.unit_normal_vector()
        self.assertAlmostEqual(normal_unit_v2D_1.norm(), 1.)


if __name__ == '__main__':
    unittest.main(verbosity=0)
