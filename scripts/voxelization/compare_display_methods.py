"""
Comparison of display methods
"""
import time

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.primitives3d import Cylinder, Sphere
from volmdlr.voxelization import Voxelization

VOXEL_SIZE = 0.02

# Create a volume model
sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
cylinder = Cylinder(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.Z3D, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([sphere, cylinder])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, name="Voxelization")

# Display the result -- fully triangulated
# With this method, each individual voxel is triangulated
start = time.perf_counter()
primitive_fully_triangulated = voxelization.to_closed_triangle_shell()
print(f"Fully triangulated shell created in {round((time.perf_counter() - start) * 1000)}ms")

start = time.perf_counter()
primitive_fully_triangulated.babylonjs()
print(f"Fully triangulated shell displayed in {round((time.perf_counter() - start) * 1000)}ms")

area = sum(primitive.area() for primitive in primitive_fully_triangulated.primitives)
print(f"Fully triangulated shell primitives area: {area}")

n_triangles = len(primitive_fully_triangulated.triangulation().triangles)
print(f"Fully triangulated shell number of triangles: {n_triangles}")

# Display the result -- simplified faces
# With this method, some PlaneFace3D are defined to represent the entire voxelization and avoid over triangulation
start = time.perf_counter()
primitive_simplified = voxelization.to_closed_shell()
print(f"\nSimplified shell created in {round((time.perf_counter() - start) * 1000)}ms")

start = time.perf_counter()
primitive_simplified.babylonjs()
print(f"Simplified shell displayed in {round((time.perf_counter() - start) * 1000)}ms")

area = sum(primitive.area() for primitive in primitive_simplified.primitives)
print(f"Simplified shell primitives area: {area}")

n_triangles = len(primitive_simplified.triangulation().triangles)
print(f"Simplified shell number of triangles: {n_triangles}")
