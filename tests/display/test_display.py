"""
Unit testing of displaying various geomertry.
"""

import os
import unittest

import volmdlr.primitives3d as p3d
from volmdlr import O2D, O3D, OXYZ, X3D, Point2D, Point3D
from volmdlr.curves import Circle2D
from volmdlr.display import Mesh3D
from volmdlr.step import Step
from volmdlr.wires import Contour2D

FOLDER = os.path.dirname(os.path.realpath(__file__))


class TestDisplayPrimitives3D(unittest.TestCase):
    def setUp(self):
        self.contour_points = [Point2D(0.0, 0.1), Point2D(0.1, 2.1), Point2D(-1.4, 0.5)]
        self.contour2d = Contour2D.from_points(self.contour_points)

    def test_diplay_block(self):
        block = p3d.Block(frame=OXYZ)
        block.babylonjs()

    def test_display_extruded_profile(self):
        extruded_profile = p3d.ExtrudedProfile(
            frame=OXYZ, outer_contour2d=self.contour2d, inner_contours2d=[], extrusion_length=1.0
        )
        extruded_profile.babylonjs()

    def test_display_revolved_profile(self):
        revolved_profile = p3d.RevolvedProfile(frame=OXYZ, contour2d=self.contour2d, axis_point=O3D, axis=X3D)
        revolved_profile.babylonjs()

    def test_display_cylinder(self):
        cylinder = p3d.Cylinder(frame=OXYZ, radius=1.0, length=1.0)
        cylinder.babylonjs()

    def test_display_hollow_cylider(self):
        hollow_cylinder = p3d.HollowCylinder(frame=OXYZ, inner_radius=0.5, outer_radius=1.0, length=1.0)
        hollow_cylinder.babylonjs()

    def test_display_sphere(self):
        sphere = p3d.Sphere(O3D, 1.0)
        sphere.babylonjs()

    def test_display_sweep(self):
        sweep = p3d.Sweep(
            contour2d=Contour2D.from_circle(Circle2D.from_center_and_radius(O2D, 0.1)),
            wire3d=p3d.OpenRoundedLineSegments3D(
                points=[Point3D(0, 0, 0), Point3D(0, 0, 1), Point3D(0, 1, 1)], radius={"1": 0.2}
            ),
        )
        sweep.babylonjs()


class TestDisplayStep(unittest.TestCase):
    def setUp(self):
        self.block_path = os.path.join(FOLDER, "..", "..", "scripts", "step", "block.step")
        self.cheese_path = os.path.join(FOLDER, "..", "..", "scripts", "step", "cheese.step")
        self.bracket2_path = os.path.join(FOLDER, "..", "..", "scripts", "step", "bracket2.step")
        self.engine_path = os.path.join(FOLDER, "..", "..", "scripts", "step", "engine.step")

    def test_display_block(self):
        volume_model = Step.from_file(self.block_path).to_volume_model()
        volume_model.babylonjs()

    def test_display_cheese(self):
        volume_model = Step.from_file(self.cheese_path).to_volume_model()
        volume_model.babylonjs()

    def test_display_bracket2_path(self):
        volume_model = Step.from_file(self.bracket2_path).to_volume_model()
        volume_model.babylonjs()

    def test_display_engine_path(self):
        volume_model = Step.from_file(self.engine_path).to_volume_model()
        volume_model.babylonjs()


class TestDisplaySTL(unittest.TestCase):
    def setUp(self):
        self.cube_ascii_path = os.path.join(FOLDER, "..", "..", "scripts", "stl", "cube_ascii.stl")
        self.double_space_path = os.path.join(FOLDER, "..", "..", "scripts", "stl", "double_space.stl")
        self.simple_path = os.path.join(FOLDER, "..", "..", "scripts", "stl", "simple.stl")

    def test_display_cube_ascii(self):
        mesh = Mesh3D.from_stl_file(self.cube_ascii_path)
        mesh.babylonjs()

        # Split shared vertices for better shadow rendering
        mesh = mesh.split_shared_vertices()
        mesh.babylonjs()

    def test_display_double_space(self):
        mesh = Mesh3D.from_stl_file(self.double_space_path)
        mesh.babylonjs()

        # Split shared vertices for better shadow rendering
        mesh = mesh.split_shared_vertices()
        mesh.babylonjs()

    def test_display_simple(self):
        mesh = Mesh3D.from_stl_file(self.simple_path)
        mesh.babylonjs()

        # Split shared vertices for better shadow rendering
        mesh = mesh.split_shared_vertices()
        mesh.babylonjs()
