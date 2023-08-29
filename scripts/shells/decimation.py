"""
Showcase of decimation of an STL file.
"""
from volmdlr.stl import Stl

STL_FILE = "../stl/Stanford_Bunny_sample.stl"
# STL_FILE = "../stl/simple.stl"

stl_model = Stl.load_from_file(STL_FILE)
closed_shell = stl_model.to_closed_shell()

decimated_shell = closed_shell.decimate(1000, verbose=True)

decimated_shell.babylonjs()
