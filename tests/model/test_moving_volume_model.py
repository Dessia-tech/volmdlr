import unittest
from copy import deepcopy

import volmdlr
from volmdlr.primitives3d import Block
from volmdlr.model import VolumeModel, MovingVolumeModel


class TestMovingVolumeModel(unittest.TestCase):
    def setUp(self):
        self.frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.O3D, volmdlr.Vector3D(1.0, 1.0, 1.0))
        self.block1 = Block(volmdlr.OXYZ)
        self.block2 = Block(self.frame)
        self.primitives = [self.block1, self.block2]
        self.volume_model = VolumeModel(deepcopy(self.primitives))
        self.step_frames = [[self.frame.translation(i * volmdlr.Vector3D(1.0, 1.0, 1.0)) for i in range(1, 3)]]
        self.moving_volume_model = MovingVolumeModel(self.primitives, self.step_frames)

    def test_is_consistent(self):
        self.assertTrue(self.moving_volume_model.is_consistent())

    def test_step_volume_model(self):
        step_model = self.moving_volume_model.step_volume_model(0)
        self.assertIsInstance(step_model, VolumeModel)

    def test_babylon_data(self):
        babylon_data = self.moving_volume_model.babylon_data()
        self.assertIsInstance(babylon_data, dict)
        self.assertIn('meshes', babylon_data)
        self.assertIn('max_length', babylon_data)
        self.assertIn('center', babylon_data)
        self.assertIn('steps', babylon_data)


if __name__ == "__main__":
    unittest.main()
