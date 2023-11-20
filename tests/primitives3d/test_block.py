import unittest

import volmdlr
from volmdlr.primitives3d import Block
from volmdlr.core import BoundingBox


class TestBlock(unittest.TestCase):
    def setUp(self) -> None:
        self.block = Block(frame=volmdlr.OXYZ)

    def test_from_bounding_box(self):
        bounding_box = BoundingBox(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)

        block = Block.from_bounding_box(bounding_box)
        self.assertEqual(self.block.frame, block.frame)

    def test_get_bounding_box(self):
        bbox1 = BoundingBox(
            -0.0271021, 1.0271020999999998, -0.4024024, 0.8024024000000001, -0.0261011, 0.5261011000000001
        )
        block = Block.from_bounding_box(bbox1)
        bbox2 = block.bounding_box

        self.assertTrue(bbox1.is_close(bbox2))
