import time

from volmdlr.primitives3d import Block
from volmdlr.discrete_representation import MatrixBasedVoxelization, EncodedOctreeBasedVoxelization
from volmdlr.core import VolumeModel
import volmdlr

block = Block(volmdlr.OXYZ)
mesh = block.triangulation()
# block.babylonjs()

VOXEL_SIZE = 0.01

t0 = time.perf_counter()
octree = EncodedOctreeBasedVoxelization.from_mesh_data(mesh.vertices, mesh.triangles, VOXEL_SIZE)
t1 = time.perf_counter()
matrix = MatrixBasedVoxelization.from_mesh_data(mesh.vertices, mesh.triangles, VOXEL_SIZE)
t2 = time.perf_counter()

print(f"Octree: {t1 - t0:.3f} s ; {octree._octree.nbytes / 1000} kB")
print(f"Matrix: {t2 - t1:.3f} s ; {matrix.matrix.nbytes / 1000} kB")

octree_old = matrix.to_octree_based_voxelization()

t0 = time.perf_counter()
a = octree.to_non_homogeneous_point_based_voxelizations()
t1 = time.perf_counter()
b = octree_old.to_non_homogeneous_point_based_voxelizations()
t2 = time.perf_counter()

print(f"New octree: {t1 - t0:.3f} s")
print(f"Old octree: {t2 - t1:.3f} s")


# volume_model = VolumeModel([_.to_mesh().split_shared_vertices() for _ in a])
# volume_model.babylonjs()
