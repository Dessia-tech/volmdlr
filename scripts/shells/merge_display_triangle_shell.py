"""
Testing merge of DisplayTriangleShell class
"""
from numpy import array
from volmdlr.shells import DisplayTriangleShell3D
from volmdlr.model import VolumeModel

indices_1 = array(
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

positions_1 = array(
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


indices_2 = array(
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

positions_2 = array(
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


display_triangle_shell_1 = DisplayTriangleShell3D(positions_1, indices_1, name="1")
display_triangle_shell_2 = DisplayTriangleShell3D(positions_2, indices_2, name="2")
display_triangle_shell_3 = display_triangle_shell_1 + display_triangle_shell_2

volume_model = VolumeModel([display_triangle_shell_1, display_triangle_shell_2, display_triangle_shell_3])
volume_model.babylonjs()
