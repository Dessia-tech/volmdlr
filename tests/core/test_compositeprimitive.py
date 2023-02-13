import unittest
import volmdlr
import volmdlr.edges
import volmdlr.core


class TestCompositePrimitive(unittest.TestCase):
    def setUp(self):
        self.primitive1 = volmdlr.edges.LineSegment2D(volmdlr.O2D, volmdlr.Point2D(volmdlr.TWO_PI, 0))
        self.primitive2 = volmdlr.edges.LineSegment2D(
            volmdlr.Point2D(volmdlr.TWO_PI, 0), volmdlr.Point2D(volmdlr.TWO_PI, 0.003)
        )
        self.primitives = [self.primitive1, self.primitive2]
        self.composite = volmdlr.core.CompositePrimitive(self.primitives, name="Composite")

    def test_primitive_to_index(self):
        self.assertEqual(self.composite.primitive_to_index(self.primitive1), 0)
        self.assertEqual(self.composite.primitive_to_index(self.primitive2), 1)

    def test_update_basis_primitives(self):
        self.composite.update_basis_primitives()
        self.assertListEqual(self.composite.basis_primitives, self.primitives)


if __name__ == "__main__":
    unittest.main()
