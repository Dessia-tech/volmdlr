#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

import numpy as npy
from volmdlr import geometry, Y3D

euler_angles = (0.1, 0.2, -0.15)
M = geometry.euler_angles_to_transfer_matrix(*euler_angles)

assert npy.isclose(geometry.transfer_matrix_to_euler_angles(M), euler_angles).all()

print(geometry.direction_to_euler_angles(Y3D))
    
print(geometry.cos_image(0, 1))

print(geometry.sin_image(-1, 0))
