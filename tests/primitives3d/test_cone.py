import unittest

import volmdlr
from volmdlr.primitives3d import Cone


class TestCone(unittest.TestCase):
    def setUp(self) -> None:
        self.cone = Cone(frame=volmdlr.OXYZ, radius=0.1, length=0.1)
