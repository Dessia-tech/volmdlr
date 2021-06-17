
import volmdlr.stl as vmstl
import volmdlr.cloud
import volmdlr.core
# import volmdlr as vm
# import volmdlr.wires as vmw
# import volmdlr.faces as vmf
# import volmdlr.edges as vme
# import matplotlib.pyplot as plt

# import os

# shells = []
# path = os.getcwd()
# # print(path)
# for stl_file in [
#         # 'a320.stl',
#                 # 'a320_ENGINE_RIGHT.stl',
#                 'a320_FAN_RIGHT.stl',
#         #         'a320_RIGHT_WING.stl',
#         #         'a320_RUDDER.stl',
#         #         'a320_STABILO_RIGHT.stl'
#                   ]:
#     # print('start')
#     # volum = volmdlr.core.VolumeModel(cloud_faces)
#     # print('saving file' + stl_file)
#     # volum.save_to_file(stl_file)
#     # print('len(cloud_faces)', len(cloud_faces))
#     # faces.extend(cloud_faces)
#     # print()
#     # points = volmdlr.stl.Stl.from_file_points(stl_file)
#     # print('size', len(points))
#     # pointcloud3d = volmdlr.cloud.PointCloud3D(points)
#     # shell3d = pointcloud3d.to_shell()
    
    
#     # stl = vmstl.Stl.from_file('/home/dasilva/volmdlr/scripts/stl/'+stl_file)
#     stl = vmstl.Stl.from_file(stl_file)
#     # print(stl)
#     shell = stl.to_closed_shell()
#     # shell.babylonjs()
#     shells.append(shell)
#     stl.extract_points()

#     # cloud = volmdlr.cloud.PointCloud3D.from_stl(path + "/" + stl_file)
#     # cloud_faces = cloud.subdescription_2d()
#     # cloud_faces.babylonjs()

    
# volum = volmdlr.core.VolumeModel(shells)
# volum.babylonjs()
        
# stl_file = vmstl.Stl.from_file("/home/dasilva/Downloads/sewingclean_304011549R--A.stl")
# stl = vmstl.Stl.from_file("/home/dasilva/Downloads/sewingclean_304011549R--A.stl")
# stl = vmstl.Stl.from_file("/home/dasilva/volmdlr/scripts/stl/a320_ENGINE_RIGHT.stl")
# shell = stl.to_closed_shell()
# volum = volmdlr.core.VolumeModel([shell])
# volum.babylonjs()

# stl_file = vmstl.Stl.from_file("/home/dasilva/Downloads/sewingclean_304011549R--A.stl")
# stl_file = vmstl.Stl.from_file("a320_ENGINE_RIGHT.stl")
points = volmdlr.stl.Stl.from_file_points('a320_ENGINE_RIGHT.stl')
print('size', len(points))
pointcloud3d = volmdlr.cloud.PointCloud3D(points)
shell3d = pointcloud3d.to_shell()
volum = volmdlr.core.VolumeModel([shell3d])
volum.babylonjs()