
import volmdlr.stl as vmstl
import volmdlr.cloud
# -*- coding: utf-8 -*-
"""
"""

# import volmdlr.stl as vmstl
# import volmdlr.cloud
import volmdlr.core
# import volmdlr as vm
# import volmdlr.wires as vmw
# import volmdlr.faces as vmf
# import volmdlr.edges as vme
# import matplotlib.pyplot as plt

import os
# import numpy as np
# from scipy.spatial import ConvexHull
shells = []
path = os.getcwd()


for stl_file in [
                'simple.stl',
                'cube_ascii.stl'
                  ]:
    stl_path = os.path.join('stl', stl_file)
    # print('start')
    # volum = volmdlr.core.VolumeModel(cloud_faces)
    # print('saving file' + stl_file)
    # volum.save_to_file(stl_file)
    # # print('len(cloud_faces)', len(cloud_faces))
    # faces.extend(cloud_faces)
    # print()

    stl = vmstl.Stl.from_file(stl_path)
    shell = stl.to_closed_shell()
    shell.alpha = 0.3
    assert len(shell.faces)
    # shell.babylonjs()
    shells.append(shell)
    # stl.extract_points()

    # cloud = volmdlr.cloud.PointCloud3D.from_stl(path + "/" + stl_file)
    # cloud_faces = cloud.subdescription_2d()
    # cloud_faces.babylonjs()
    # list_points = vmstl.Stl.from_file_points(stl_file)
    list_points = stl.extract_points_BIS()
    if len(list_points) > 1:
        pointcloud3d = volmdlr.cloud.PointCloud3D(list_points)
        # polygons2d = pointcloud3d.to_shell()
        # pointcloud3d.plot()
        shell2 = pointcloud3d.to_shell(resolution=15)
        shell2.alpha = 0.6
        shell2.color = (1, 0.1, 0.1)
        shells.append(shell2)


volume_model = volmdlr.core.VolumeModel(shells)
volume_model.babylonjs()



