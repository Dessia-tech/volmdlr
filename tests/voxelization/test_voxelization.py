import unittest

import volmdlr
from volmdlr.primitives3d import Sphere, Block
from volmdlr.voxelization import Voxelization


class TestVoxelization(unittest.TestCase):
    def setUp(self):
        self.sphere = Sphere(center=volmdlr.O3D, radius=0.02)
        self.block = Block(volmdlr.OXYZ)

    def test_voxelize_sphere(self):
        voxelized_sphere = Voxelization.from_closed_shell(self.sphere, 0.005)
        self.assertEqual(len(voxelized_sphere), 284)

    def test_voxelize_block(self):
        voxelized_block = Voxelization.from_closed_shell(self.block, 0.2)
        self.assertEqual(len(voxelized_block), 152)
