"""
Shell intersection detection.
"""
import math

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.step import Step
from volmdlr.faces import Triangle3D
from volmdlr.shells import OpenTriangleShell3D
from volmdlr.discrete_representation import OctreeBasedVoxelization, PointBasedVoxelization

# Load
shell_1 = Step.from_file("../../step/block.step").to_volume_model().get_shells()[0]
shell_1.color = (0, 1, 0)
shell_2 = shell_1.translation(volmdlr.Vector3D(0.001, 0.001, 0.001))
shell_2 = shell_2.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
shell_2.color = (0, 0, 1)

volume_model = VolumeModel([shell_1, shell_2])
volume_model.babylonjs()

# face_combinations = shell1.intersecting_faces_combinations(shell2)

VOXEL_SIZE = 0.0001
voxelization_1 = OctreeBasedVoxelization.from_shell(shell_1, VOXEL_SIZE)
voxelization_2 = OctreeBasedVoxelization.from_shell(shell_2, VOXEL_SIZE)
intersection = voxelization_1.intersection(voxelization_2)

triangle_combinations = intersection._get_intersections_voxel_centers(len(voxelization_1._triangles))

for combination, voxel_centers in triangle_combinations.items():
    location = PointBasedVoxelization(voxel_centers, VOXEL_SIZE)
    location_shell = location.to_closed_triangle_shell()
    location_shell.color = (1, 0, 0)
    location_shell.alpha = 0.5

    triangle_face_1 = OpenTriangleShell3D([Triangle3D(
        volmdlr.Point3D(*intersection._triangles[combination[0]][0]),
        volmdlr.Point3D(*intersection._triangles[combination[0]][1]),
        volmdlr.Point3D(*intersection._triangles[combination[0]][2]),
    )])
    triangle_face_1.color = (0, 1, 0)
    triangle_face_2 = OpenTriangleShell3D([Triangle3D(
        volmdlr.Point3D(*intersection._triangles[combination[1]][0]),
        volmdlr.Point3D(*intersection._triangles[combination[1]][1]),
        volmdlr.Point3D(*intersection._triangles[combination[1]][2]),
    )])
    triangle_face_2.color = (0, 0, 1)

    intersection_volume_model = VolumeModel([location_shell, triangle_face_1, triangle_face_2])
    intersection_volume_model.babylonjs()
