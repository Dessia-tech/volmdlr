import unittest

import volmdlr


class TestFrame2D(unittest.TestCase):
    local_frame = volmdlr.Frame2D(volmdlr.Point2D(0.1, 0.3), volmdlr.Y2D, volmdlr.X2D)

    def test_global_to_local_coordinates(self):
        vector_global = volmdlr.Vector2D(3, 4)
        vector_local = self.local_frame.global_to_local_coordinates(vector_global)

        # Check that the converted vector has the expected coordinates
        self.assertEqual(vector_local.x, 3.7)
        self.assertEqual(vector_local.y, 2.9)

    def test_local_to_global_coordinates(self):
        vector_local = volmdlr.Vector2D(3.7, 2.9)
        vector_global = self.local_frame.local_to_global_coordinates(vector_local)

        # Check that the converted vector has the expected coordinates
        self.assertEqual(vector_global.x, 3)
        self.assertEqual(vector_global.y, 4)


if __name__ == "__main__":
    unittest.main()
