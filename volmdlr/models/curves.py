"""
Some Volmdlr models for curves.

"""
import volmdlr
from volmdlr import curves

vector1 = volmdlr.Vector3D(1, 1, 1)
vector1 = vector1.unit_vector()
vector2 = vector1.deterministic_unit_normal_vector()
vector3 = vector1.cross(vector2)
frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)

circle3d = curves.Circle3D(frame, 1)
