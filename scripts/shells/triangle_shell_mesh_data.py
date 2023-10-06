"""
Showcase of TriangleShell3D import / export.
"""

from volmdlr.stl import Stl
from volmdlr.shells import ClosedTriangleShell3D

triangle_shell_1 = Stl.load_from_file("../stl/simple.stl").to_closed_shell()

mesh_data = triangle_shell_1.to_mesh_data()
triangle_shell_2 = ClosedTriangleShell3D.from_mesh_data(mesh_data[0], mesh_data[1])

triangle_shell_1.save_to_file("triangle_shell.json")
triangle_shell_3 = ClosedTriangleShell3D.load_from_file("triangle_shell.json")

print(triangle_shell_3 == triangle_shell_1)
