# -*- coding: utf-8 -*-
"""

"""

import volmdlr.cloud
import volmdlr.core
# import volmdlr as vm
# import volmdlr.wires as vmw
# import volmdlr.faces as vmf
# import volmdlr.edges as vme
# import matplotlib.pyplot as plt

import os

faces = []
path = os.getcwd()
for stl_file in [#'a320.stl',
                  # 'a320_ENGINE.stl',
                  # 'a320_FAN.stl',
                  # 'a320_LEFT_WING.stl',
                  # 'a320_RIGHT_WING.stl',
                  # 'a320_RUDDER.stl',
                   'a320_STABILO_LEFT.stl',
                  # 'a320_STABILO_RIGHT.stl'
                  ]:
    cloud = volmdlr.cloud.PointCloud3D.from_stl(path + "\\" + stl_file)
    cloud_faces = cloud.subdescription_2d()
    
    volum = volmdlr.core.VolumeModel(cloud_faces)
    # print('len(cloud_faces)', len(cloud_faces))
    faces.extend(cloud_faces)
    
    
volum = volmdlr.core.VolumeModel(faces)
volum.babylonjs()
        

