"""
Shell intersection detection.
"""
import math
import time

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.step import Step
from volmdlr.discrete_representation import OctreeBasedVoxelization


# Load
shell_1 = Step.from_file("../../step/engine.step").to_volume_model().get_shells()[0]
shell_1.color = (0, 1, 0)
shell_2 = shell_1.translation(volmdlr.Vector3D(0.001, 0.001, 0.001))
shell_2 = shell_2.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
shell_2.color = (0, 0, 1)

volume_model = VolumeModel([shell_1, shell_2])
# volume_model.babylonjs()

# t0 = time.perf_counter()
# face_combinations_shell = shell_1.intersecting_faces_combinations(shell_2)
# t1 = time.perf_counter()
#
# print(f"Shell3D.intersecting_faces_combinations: {(t1 - t0) * 1000:.3f}ms")

VOXEL_SIZE = 0.00005

t0 = time.perf_counter()
face_combinations_octree = OctreeBasedVoxelization.intersecting_faces_combinations(shell_1, shell_2, VOXEL_SIZE)
t1 = time.perf_counter()

print(f"OctreeBasedVoxelization.intersecting_faces_combinations: {(t1 - t0) * 1000:.3f}ms")


def display_face_combination(face_combination):
    location = face_combination[1]
    location_shell = location.to_closed_triangle_shell()
    location_shell.color = (1, 0, 0)
    location_shell.alpha = 0.5

    face_1 = face_combination[0][0]
    face_1.color = (0, 1, 0)

    face_2 = face_combination[0][1]
    face_2.color = (0, 0, 1)

    intersection_volume_model = VolumeModel([location_shell, face_1, face_2])
    intersection_volume_model.babylonjs()
