#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 14:11:53 2021

@author: dasilva
"""

import volmdlr.stl as vmstl
import volmdlr.cloud
import volmdlr.core
# import volmdlr as vm
# import volmdlr.wires as vmw
# import volmdlr.faces as vmf
# import volmdlr.edges as vme
# import matplotlib.pyplot as plt

import os
import numpy as np
from scipy.spatial import ConvexHull
shells = []
path = os.getcwd()
stl_path = path
os.chdir(stl_path)
print('os.listdir(stl_path) :', os.listdir(stl_path))
# for stl_file in os.listdir(stl_path):
#     if len(stl_file.split('.')) > 1:
#         if stl_file.split('.')[1] == 'stl':
#             stl = volmdlr.stl.Stl.from_file(stl_file)
#             stl.name = stl_file
#             # list_points = stl.extract_points_BIS()
#             list_points = stl.extract_points()
#             pointcloud3d = volmdlr.cloud.PointCloud3D(list_points)
#             shell = pointcloud3d.alpha_shape(alpha = 0.2, number_point_samples= 500)
#             shell.color = (1, 0.1, 0.1)
#             shell.alpha = 0.6
#             volum = volmdlr.core.VolumeModel([shell])
#             volum.babylonjs()
#             shells.append(shell)

    
# volum = volmdlr.core.VolumeModel(shells)
# volum.babylonjs()
        
stl = volmdlr.stl.Stl.from_file('a320.stl')
stl.name = 'a320.stl'
list_points = stl.extract_points()
# list_points = stl.extract_points_BIS()
point_cloud = volmdlr.cloud.PointCloud3D(list_points)
shell = point_cloud.alpha_shape(alpha = 0.4, number_point_samples= 10000)
shell.color = (1, 0.1, 0.1)
shell.alpha = 0.6
volmdlr.core.VolumeModel([shell]).babylonjs()
