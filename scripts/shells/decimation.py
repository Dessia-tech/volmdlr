"""
Showcase of decimation of a STL file.
"""
import os
import urllib.request

from volmdlr.stl import Stl

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
stl_model = Stl.load_from_file(FILE_NAME)
closed_triangle_shell = stl_model.to_closed_shell()

# Decimate and show model
decimated_closed_triangle_shell = closed_triangle_shell.decimate(50, verbose=True)
decimated_closed_triangle_shell.babylonjs()
