# -*- coding: utf-8 -*-
"""

"""

import volmdlr.cloud
import volmdlr.core
# import volmdlr as vm
# import volmdlr.wires as vmw
import volmdlr.faces as vmf
# import volmdlr.edges as vme
# import matplotlib.pyplot as plt

import os

volum = []
path = os.getcwd()
for stl_file in [#'a320.stl',
                  'a320_ENGINE_RIGHT.stl',
                  # 'a320_ENGINE_LEFT.stl',
                  'a320_FAN_RIGHT.stl',
                  # 'a320_FAN_LEFT.stl',
                  # 'a320_LEFT_WING.stl',
                  'a320_RIGHT_WING.stl',
                  # 'a320_RUDDER.stl',
                  # 'a320_STABILO_LEFT.stl',
                  # 'a320_STABILO_RIGHT.stl'
                  ]:
    print('start')
    cloud = volmdlr.cloud.PointCloud3D.from_stl(path + "/" + stl_file)
    shell = cloud.to_shell()
    
    print('saving file' + stl_file)
    shell.save_to_file(stl_file)
    # print('len(cloud_faces)', len(cloud_faces))
    volum.append(shell)
    print()
    
volumodel=volmdlr.core.VolumeModel(volum)    
    
    
# volum = volmdlr.core.VolumeModel(faces)
# volum.babylonjs()
        
# for json_file in ['a320_ENGINE_RIGHT.json',
                  # 'a320_FAN_RIGHT.json',
                   # 'a320_RIGHT_WING.json']:
    # vol = volmdlr.core.VolumeModel.load_from_file(json_file)
    # faces.extend(vol.primitives)
# 
# faces.pop(6230)
    # 
# volass = volmdlr.core.VolumeModel(faces)



# primitive_pb = volass.primitives[6230]