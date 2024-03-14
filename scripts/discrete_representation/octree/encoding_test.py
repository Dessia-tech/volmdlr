import time

from volmdlr.primitives3d import Block
from volmdlr.discrete_representation_compiled import encode
from volmdlr.discrete_representation import MatrixBasedVoxelization
import volmdlr

block = Block(volmdlr.OXYZ)
mesh = block.triangulation()
# block.babylonjs()

VOXEL_SIZE = 0.001

t0 = time.perf_counter()
encoded = encode(mesh.vertices, mesh.triangles, VOXEL_SIZE)
t1 = time.perf_counter()
matrix = MatrixBasedVoxelization.from_mesh_data(mesh.vertices, mesh.triangles, VOXEL_SIZE)
t2 = time.perf_counter()

print(f"Octree: {t1 - t0:.3f} s ; {encoded.nbytes / 1000} kB")
print(f"Matrix: {t2 - t1:.3f} s ; {matrix.matrix.nbytes / 1000} kB")
