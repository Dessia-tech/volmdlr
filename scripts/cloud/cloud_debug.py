# -*- coding: utf-8 -*-
"""

"""

import volmdlr.cloud
import volmdlr.core
import volmdlr as vm
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.edges as vme
import matplotlib.pyplot as plt

# import os

# faces = []
# path = ''
# for file in os.listdir(path) :
    # cloud = volmdlr.cloud.PointCloud3D.from_stl(path + file)
#     cloud_faces = cloud.subdescription_2d()
    
#     volum = volmdlr.core.VolumeModel(cloud_faces)
#     faces.extend(cloud_faces)
#     volum.save_to_file(file)
    
# volum = volmdlr.core.VolumeModel(faces)
        
################### READ STL

# cloud3d = volmdlr.cloud.PointCloud3D.from_stl(path)
# faces = cloud3d.subdescription_2d(resolution=20)

# volum = volmdlr.core.VolumeModel(faces)

# volum.babylonjs()


# volumodel = volmdlr.core.VolumeModel.load_from_file('name.json')
