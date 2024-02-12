import math
import unittest
from copy import deepcopy
import volmdlr
from volmdlr.primitives3d import Block
from volmdlr.core import VolumeModel, BoundingBox


class TestVolumeModel(unittest.TestCase):
    def setUp(self):
        self.frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.O3D, volmdlr.Vector3D(1.0, 1.0, 1.0))
        self.block1 = Block(volmdlr.OXYZ)
        self.block2 = Block(self.frame)
        self.primitives = [self.block1, self.block2]
        self.volume_model = VolumeModel(deepcopy(self.primitives))

    def test_eq(self):
        different_volume_model = VolumeModel([deepcopy(self.block2)])
        self.assertNotEqual(self.volume_model, different_volume_model)
        self.assertEqual(self.volume_model, self.volume_model)

    def test_volume(self):
        self.assertEqual(self.volume_model.volume(), sum(p.volume() for p in self.primitives))

    def test_rotation(self):
        center = volmdlr.Point3D(0.0, 0.0, 0.0)
        axis = volmdlr.Vector3D(0.0, 0.0, 1.0)
        angle = math.pi / 4
        rotated_volume_model = self.volume_model.rotation(center, axis, angle)

        for p1, p2 in zip(rotated_volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.rotation(center, axis, angle)
            self.assertEqual(p1, p2)

    def test_translation(self):
        offset = volmdlr.Vector3D(1.0, 1.0, 1.0)
        translated_volume_model = self.volume_model.translation(offset)

        for p1, p2 in zip(translated_volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.translation(offset)
            self.assertEqual(p1, p2)

    def test_frame_mapping(self):
        frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.Point3D(1.0, 1.0, 1.0), volmdlr.Vector3D(1.0, 1.0, 1.0))
        side = "old"
        mapped_volume_model = self.volume_model.frame_mapping(frame, side)

        for p1, p2 in zip(mapped_volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.frame_mapping(frame, side)
            self.assertEqual(p1, p2)

    def test_bounding_box(self):
        self.assertEqual(
            self.volume_model.bounding_box, BoundingBox.from_bounding_boxes([p.bounding_box for p in self.primitives])
        )

    def test_plot(self):
        volume_model_plot_lines = self.volume_model.plot().lines
        primitives_plot_lines = [line for p in self.primitives for line in p.plot().lines]
        self.assertEqual(len(volume_model_plot_lines), len(primitives_plot_lines))

    def test_to_stl(self):
        self.volume_model.to_stl("volume_model")

    def test_to_mesh3d(self):
        mesh = self.volume_model.to_mesh()

        self.assertEqual(20, mesh.n_vertices)
        self.assertEqual(24, mesh.n_triangles)


if __name__ == "__main__":
    unittest.main()
