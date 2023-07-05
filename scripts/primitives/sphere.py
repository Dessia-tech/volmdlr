"""
Demo script of creating and displaying a Sphere.
"""
import volmdlr
from volmdlr.primitives3d import Sphere

sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
sphere.babylonjs()
