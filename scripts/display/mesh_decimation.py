"""
Showcase of decimation of a STL file.
"""
import os
import urllib.request

from volmdlr.display import Mesh3D

# Stanfrod Bunny model
MODEL_URL = "https://upload.wikimedia.org/wikipedia/commons/4/43/Stanford_Bunny.stl"
FILE_NAME = "Stanford_Bunny.stl"

# Check if the STL file already exists
if not os.path.exists(FILE_NAME):
    urllib.request.urlretrieve(MODEL_URL, FILE_NAME)
    print(f"File downloaded to {FILE_NAME}")
else:
    print(f"File already exists at {FILE_NAME}. Skipping download.")

# Load STL model using volmdlr
mesh = Mesh3D.from_stl_file(FILE_NAME)

# Decimate and show model
decimated_mesh = mesh.decimate(500, verbose=True)
decimated_mesh = decimated_mesh.split_shared_vertices()  # for display purpose
decimated_mesh.babylonjs()
