"""
Shell intersection detection.
"""
import math
import time

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.step import Step
from volmdlr.faces import Triangle3D
from volmdlr.shells import OpenTriangleShell3D
from volmdlr.discrete_representation import OctreeBasedVoxelization, PointBasedVoxelization

def shell_to_triangles(shell):
    """
    Helper method to convert a Shell3D to a list of triangles.

    It uses the "triangulation" method to triangulate the Shell3D.

    :param shell: The Shell3D to convert to triangles.
    :type shell: Shell3D

    :return: The list of triangles extracted from the triangulated Shell3D.
    :rtype: List[Triangle]
    """

    face_id_by_triangle = {}

    for i, face in enumerate(shell.faces):
        triangulation = face.triangulation()

        triangles = [
            (
                (
                    float(triangulation.points[triangle[0]].x),
                    float(triangulation.points[triangle[0]].y),
                    float(triangulation.points[triangle[0]].z),
                ),
                (
                    float(triangulation.points[triangle[1]].x),
                    float(triangulation.points[triangle[1]].y),
                    float(triangulation.points[triangle[1]].z),
                ),
                (
                    float(triangulation.points[triangle[2]].x),
                    float(triangulation.points[triangle[2]].y),
                    float(triangulation.points[triangle[2]].z),
                ),
            )
            for triangle in triangulation.triangles
        ]

        for triangle in triangles:
            face_id_by_triangle[triangle] = face

    return face_id_by_triangle


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

VOXEL_SIZE = 0.01

t0 = time.perf_counter()
face_combinations_octree = OctreeBasedVoxelization.intersecting_faces_combinations(shell_1, shell_2, VOXEL_SIZE)
t1 = time.perf_counter()

print(f"OctreeBasedVoxelization.intersecting_faces_combinations: {(t1 - t0) * 1000:.3f}ms")

# voxelization_1 = OctreeBasedVoxelization.from_shell(shell_1, VOXEL_SIZE)
# voxelization_2 = OctreeBasedVoxelization.from_shell(shell_2, VOXEL_SIZE)
# intersection = voxelization_1.intersection(voxelization_2)
#
# triangle_combinations = intersection._get_intersections_voxel_centers(len(voxelization_1._triangles))
#
# face_id_by_triangle_1 = shell_to_triangles(shell_1)
# face_id_by_triangle_2 = shell_to_triangles(shell_2)
#
# face_combinations = {}

# for combination, voxel_centers in triangle_combinations.items():
#     face_1 = face_id_by_triangle_1[intersection._triangles[combination[0]]]
#     face_2 = face_id_by_triangle_2[intersection._triangles[combination[1]]]
#
#     if (face_1, face_2) not in face_combinations:
#         face_combinations[(face_1, face_2)] = set()
#     face_combinations[(face_1, face_2)].update(voxel_centers)
#
# for combination, voxel_centers in face_combinations.items():
#     location = PointBasedVoxelization(voxel_centers, VOXEL_SIZE)
#     location_shell = location.to_closed_triangle_shell()
#     location_shell.color = (1, 0, 0)
#     location_shell.alpha = 0.5
#
#     # face_1 = shell_1.faces[combination[0]]
#     face_1 = combination[0]
#     face_1.color = (0, 1, 0)
#     # face_2 = shell_2.faces[combination[1]]
#     face_2 = combination[1]
#     face_2.color = (0, 0, 1)
#
#     intersection_volume_model = VolumeModel([location_shell, face_1, face_2])
#     intersection_volume_model.babylonjs()

# for combination, voxel_centers in triangle_combinations.items():
#     location = PointBasedVoxelization(voxel_centers, VOXEL_SIZE)
#     location_shell = location.to_closed_triangle_shell()
#     location_shell.color = (1, 0, 0)
#     location_shell.alpha = 0.5
#
#     triangle_face_1 = OpenTriangleShell3D([Triangle3D(
#         volmdlr.Point3D(*intersection._triangles[combination[0]][0]),
#         volmdlr.Point3D(*intersection._triangles[combination[0]][1]),
#         volmdlr.Point3D(*intersection._triangles[combination[0]][2]),
#     )])
#     triangle_face_1.color = (0, 1, 0)
#     triangle_face_2 = OpenTriangleShell3D([Triangle3D(
#         volmdlr.Point3D(*intersection._triangles[combination[1]][0]),
#         volmdlr.Point3D(*intersection._triangles[combination[1]][1]),
#         volmdlr.Point3D(*intersection._triangles[combination[1]][2]),
#     )])
#     triangle_face_2.color = (0, 0, 1)
#
#     intersection_volume_model = VolumeModel([location_shell, triangle_face_1, triangle_face_2])
#     intersection_volume_model.babylonjs()
