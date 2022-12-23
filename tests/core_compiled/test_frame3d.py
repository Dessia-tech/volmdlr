import unittest
import volmdlr


class TestFrame3D(unittest.TestCase):
    local_frame = volmdlr.Frame3D(
        volmdlr.Point3D(0.1, 0.3, 0.6),
        volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D,
    )

    def test_global_to_local_coordinates(self):
        vector_global = volmdlr.Vector3D(3, 4, 5)
        vector_local = self.local_frame.global_to_local_coordinates(vector_global)

        # Check that the converted vector has the expected coordinates
        self.assertEqual(vector_local.x, 4-0.3)
        self.assertEqual(vector_local.y, 5-0.6)
        self.assertEqual(vector_local.z, 3-0.1)

    def test_local_to_global_coordinates(self):
        vector_local = volmdlr.Vector3D(4-0.3, 5-0.6, 3-0.1)
        vector_global = self.local_frame.local_to_global_coordinates(vector_local)

        # Check that the converted vector has the expected coordinates
        self.assertEqual(vector_global.x, 3)
        self.assertEqual(vector_global.y, 4)
        self.assertEqual(vector_global.z, 5)


if __name__ == "__main__":
    unittest.main()
