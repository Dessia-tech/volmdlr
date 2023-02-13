import unittest
import matplotlib.pyplot as plt
import volmdlr
import volmdlr.edges


class TestCompositePrimitive2D(unittest.TestCase):
    def setUp(self):
        self.primitives = [
            volmdlr.edges.LineSegment2D(volmdlr.O2D, volmdlr.Point2D(volmdlr.TWO_PI, 0)),
            volmdlr.edges.LineSegment2D(
                volmdlr.Point2D(volmdlr.TWO_PI, 0),
                volmdlr.Point2D(volmdlr.TWO_PI, 0.003),
            ),
        ]
        self.composite_2d = volmdlr.core.CompositePrimitive2D(self.primitives, name="test")

    def test_plot(self):
        ax = self.composite_2d.plot()
        self.assertIsNotNone(ax)

    def test_plot_equal_aspect(self):
        ax = self.composite_2d.plot(equal_aspect=True)
        self.assertEqual(ax.get_aspect(), 1.0)

    def test_init(self):
        self.assertEqual(self.composite_2d.primitives, self.primitives)
        self.assertEqual(self.composite_2d.name, "test")

    def test_rotation(self):
        center = volmdlr.Point2D(1.2, -3.4)
        angle = 2.56
        rotated_composite_2d = self.composite_2d.rotation(center, angle)

        for p1, p2 in zip(rotated_composite_2d.primitives, self.primitives):
            self.assertEqual(p1, p2.rotation(center, angle))


if __name__ == "__main__":
    unittest.main()
