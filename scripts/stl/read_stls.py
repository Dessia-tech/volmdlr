import os

import volmdlr.cloud
import volmdlr.core
import volmdlr.stl as vmstl


path = os.getcwd()


for stl_file in [
                'simple.stl',
                'cube_ascii.stl',
                'double_space.stl'
                  ]:

    shells = []
    stl = vmstl.Stl.from_file(stl_file)
    shell = stl.to_closed_shell()
    shell.alpha = 0.3
    assert len(shell.faces)
    shells.append(shell)

    list_points = stl.extract_points_BIS()
    if len(list_points) > 50:
        pointcloud3d = volmdlr.cloud.PointCloud3D(list_points)
        shell2 = pointcloud3d.to_shell(resolution=15)
        shell2.alpha = 0.6
        shell2.color = (1, 0.1, 0.1)
        shells.append(shell2)

    volume_model = volmdlr.core.VolumeModel(shells)
    volume_model.babylonjs()
