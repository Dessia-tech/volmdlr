
import volmdlr.stl as vmstl
import volmdlr.cloud
# -*- coding: utf-8 -*-
"""
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
for stl_file in [
                # 'a320.stl',
                # 'a320_ENGINE_RIGHT.stl',
                # 'a320_FAN_RIGHT.stl',
                # 'a320_RIGHT_WING.stl',
                # 'a320_RUDDER.stl',
                'a320_STABILO_RIGHT.stl'
                  ]:
    # print('start')
    # volum = volmdlr.core.VolumeModel(cloud_faces)
    # print('saving file' + stl_file)
    # volum.save_to_file(stl_file)
    # # print('len(cloud_faces)', len(cloud_faces))
    # faces.extend(cloud_faces)
    # print()
    
    stl = vmstl.Stl.from_file(stl_file)
    # shell = stl.to_closed_shell()
    # shell.babylonjs()
    # shells.append(shell)
    # stl.extract_points()

    # cloud = volmdlr.cloud.PointCloud3D.from_stl(path + "/" + stl_file)
    # cloud_faces = cloud.subdescription_2d()
    # cloud_faces.babylonjs()
    # list_points = vmstl.Stl.from_file_points(stl_file)
    list_points = stl.extract_points_BIS()
    pointcloud3d = volmdlr.cloud.PointCloud3D(list_points)
    # polygons3d = pointcloud3d.to_shell()
    # pointcloud3d.plot()
    shells.append(pointcloud3d.to_shell())

    
volum = volmdlr.core.VolumeModel(shells)
volum.babylonjs()
        

# points = [(p.x, p.y) for p in subcloud2d[0].points]
# points = np.array(points)
# plt.plot(points[:,0], points[:,1], 'o')
# hull = ConvexHull(points)
# for simplex in hull.simplices:
#     plt.plot(points[simplex,0], points[simplex,1], 'k-')
# plt.show()

#