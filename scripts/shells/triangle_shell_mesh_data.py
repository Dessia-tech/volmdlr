"""
Showcase of TriangleShell3D to_mesh_data and from_mesh_data methods.
"""

from volmdlr.stl import Stl
from volmdlr.shells import ClosedTriangleShell3D

triangle_shell = Stl.load_from_file("../stl/simple.stl").to_closed_shell()

mesh_data = triangle_shell.to_mesh_data()
triangle_shell_2 = ClosedTriangleShell3D.from_mesh_data(mesh_data[0], mesh_data[1])
