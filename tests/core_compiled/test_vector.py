import unittest
import volmdlr


class TestVector(unittest.TestCase):
    def test_is_colinear_to(self):
        vector_ref = volmdlr.Z3D
        vector_test = volmdlr.Vector3D(1e-3, 1e-3, 1)
        self.assertTrue(vector_ref.is_colinear_to(vector_test))
        self.assertFalse(vector_ref.is_colinear_to(vector_test, abs_tol=0))

    def test_is_perpendicular_to(self):
        vector_ref = volmdlr.Z3D
        vector_test = volmdlr.Vector3D(1-6, 1, 1e-6)
        self.assertTrue(vector_ref.is_perpendicular_to(vector_test))
        self.assertFalse(vector_ref.is_perpendicular_to(vector_test, abs_tol=0))


if __name__ == '__main__':
    unittest.main()
