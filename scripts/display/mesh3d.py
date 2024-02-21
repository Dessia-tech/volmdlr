"""
Showcase of Mesh3D class.
"""
import numpy as np
from volmdlr.display import Mesh3D

positions3 = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 0.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 0.0],
        [1.0, 1.0, 1.0],
    ]
)
indices3 = np.array(
    [
        [2, 6, 7],
        [0, 4, 5],
        [1, 7, 5],
        [0, 2, 6],
        [4, 6, 7],
        [1, 3, 7],
        [0, 2, 3],
        [2, 7, 3],
        [0, 6, 4],
        [4, 7, 5],
        [0, 5, 1],
        [0, 3, 1],
    ]
)
mesh3 = Mesh3D(positions3, indices3, name="Mesh3")

positions4 = np.array(
    [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 2.0],
        [0.0, 1.0, 1.0],
        [0.0, 1.0, 2.0],
        [1.0, 0.0, 1.0],
        [1.0, 0.0, 2.0],
        [1.0, 1.0, 1.0],
        [1.0, 1.0, 2.0],
    ]
)
indices4 = np.array(
    [
        [2, 7, 3],
        [1, 7, 5],
        [0, 6, 4],
        [4, 7, 5],
        [0, 3, 1],
        [0, 2, 6],
        [4, 6, 7],
        [2, 6, 7],
        [0, 4, 5],
        [1, 3, 7],
        [0, 2, 3],
        [0, 5, 1],
    ]
)
mesh4 = Mesh3D(positions4, indices4, name="Mesh4")


# Plotting
mesh3.plot()
mesh4.plot()

# Distance
print(mesh3.minimum_distance(mesh4))

# 3D view
mesh3.babylonjs()
mesh4.babylonjs()

# Merge
mesh5 = mesh3 + mesh4  # without duplicated vertices / faces merging
# mesh5 = mesh3.merge(mesh4, merge_vertices=False, merge_triangles=False)  # equivalent code

mesh6 = mesh3 | mesh4  # with duplicated vertices / faces merging
# mesh6 = mesh3.merge(mesh4, merge_vertices=True, merge_triangles=True)  # equivalent code

print(mesh5.n_vertices, mesh5.n_triangles)  # Expected 16 / 24
print(mesh6.n_vertices, mesh6.n_triangles)  # Expected 12 / 22

mesh5.babylonjs()
mesh6.babylonjs()
