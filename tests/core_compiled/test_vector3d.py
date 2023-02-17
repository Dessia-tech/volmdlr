"""
Unittest for volmdlr.Vector3D
"""
import unittest
import volmdlr


class TestVector3D(unittest.TestCase):
    def test_vector_projection(self):
        v1 = volmdlr.Z3D
        v2 = volmdlr.Vector3D(1, 1, 1)
        v3 = volmdlr.Vector3D(1, 1, -1)
        v4 = volmdlr.X3D
        p1 = v2.vector_projection(v1)
        p2 = v3.vector_projection(v1)
        p3 = v4.vector_projection(v1)
        self.assertEqual(p1, volmdlr.Vector3D(0, 0, 1))
        self.assertEqual(p2, volmdlr.Vector3D(0, 0, -1))
        self.assertEqual(p3, volmdlr.Vector3D(0, 0, 0))


if __name__ == '__main__':
    unittest.main()
