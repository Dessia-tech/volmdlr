import math
import unittest
import volmdlr
from volmdlr.primitives3d import Block
from volmdlr.core import VolumeModel, BoundingBox
from copy import deepcopy


class TestVolumeModel(unittest.TestCase):
    def setUp(self):
        self.primitives = [Block(volmdlr.OXYZ)]
        self.volume_model = VolumeModel(deepcopy(self.primitives))

    def test_volume(self):
        self.assertEqual(self.volume_model.volume(), sum([p.volume() for p in self.primitives]))

    def test_rotation(self):
        center = volmdlr.Point3D(0.0, 0.0, 0.0)
        axis = volmdlr.Vector3D(0.0, 0.0, 1.0)
        angle = math.pi / 4
        rotated_volume_model = self.volume_model.rotation(center, axis, angle)

        for p1, p2 in zip(rotated_volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.rotation(center, axis, angle)
            self.assertEqual(p1, p2)

    def test_rotation_inplace(self):
        center = volmdlr.Point3D(0.0, 0.0, 0.0)
        axis = volmdlr.Vector3D(0.0, 0.0, 1.0)
        angle = math.pi / 4
        self.volume_model.rotation_inplace(center, axis, angle)

        for p1, p2 in zip(self.volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2.rotation_inplace(center, axis, angle)
            self.assertEqual(p1, p2)

    def test_translation(self):
        offset = volmdlr.Vector3D(1.0, 1.0, 1.0)
        translated_volume_model = self.volume_model.translation(offset)

        for p1, p2 in zip(translated_volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.translation(offset)
            self.assertEqual(p1, p2)

    def test_translation_inplace(self):
        offset = volmdlr.Vector3D(1.0, 1.0, 1.0)
        self.volume_model.translation_inplace(offset)

        for p1, p2 in zip(self.volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2.translation_inplace(offset)
            self.assertEqual(p1, p2)

    def test_frame_mapping(self):
        frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.Point3D(1.0, 1.0, 1.0), volmdlr.Vector3D(1.0, 1.0, 1.0))
        side = "old"
        mapped_volume_model = self.volume_model.frame_mapping(frame, side)

        for p1, p2 in zip(mapped_volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.frame_mapping(frame, side)
            self.assertEqual(p1, p2)

    def test_frame_mapping_inplace(self):
        frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.Point3D(1.0, 1.0, 1.0), volmdlr.Vector3D(1.0, 1.0, 1.0))
        side = "old"
        self.volume_model.frame_mapping_inplace(frame, side)

        for p1, p2 in zip(self.volume_model.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2.frame_mapping_inplace(frame, side)
            self.assertEqual(p1, p2)


if __name__ == "__main__":
    unittest.main()
