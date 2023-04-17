#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geometry unittest.
"""

import unittest

import numpy as npy

from volmdlr import Y3D, Point2D, geometry


class TestClosedShell3D(unittest.TestCase):
    def test_euler_angles(self):
        euler_angles = (0.1, 0.2, -0.15)
        M = geometry.euler_angles_to_transfer_matrix(*euler_angles)
        self.assertTrue(npy.isclose(geometry.transfer_matrix_to_euler_angles(M), euler_angles).all())

    def test_direction_to_euler(self):
        geometry.direction_to_euler_angles(Y3D)

    def test_images(self):
        self.assertAlmostEqual(geometry.cos_image(0, 1)[1], 1.)
        self.assertAlmostEqual(geometry.sin_image(-1, 0)[1], 0.)

    def test_huygens2d(self):
        true_res = (11.0, 492.0, -67.0)
        res = geometry.huygens2d(1, 2, 3, 10, Point2D(5, -2), Point2D(-2, -3))
        self.assertAlmostEqual(res[0], true_res[0])
        self.assertAlmostEqual(res[1], true_res[1])
        self.assertAlmostEqual(res[2], true_res[2])
        
        
if __name__ == '__main__':
    unittest.main()
