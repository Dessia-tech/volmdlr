import unittest
import volmdlr
import volmdlr.core
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
            self.assertNotEqual(p1, p2)
            p2 = p2.rotation(center, angle)
            self.assertEqual(p1, p2)

    def test_rotation_inplace(self):
        center = volmdlr.Point2D(-1.7, -0.3)
        angle = -8
        self.composite_2d.rotation_inplace(center, angle)

        for p1, p2 in zip(self.composite_2d.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2.rotation_inplace(center, angle)
            self.assertEqual(p1, p2)

    def test_translation(self):
        offset = volmdlr.Vector2D(0.56, -3.4)
        rotated_composite_2d = self.composite_2d.translation(offset)

        for p1, p2 in zip(rotated_composite_2d.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.translation(offset)
            self.assertEqual(p1, p2)

    def test_translation_inplace(self):
        offset = volmdlr.Vector2D(11, -0.04)
        self.composite_2d.translation_inplace(offset)

        for p1, p2 in zip(self.composite_2d.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2.translation_inplace(offset)
            self.assertEqual(p1, p2)

    def test_frame_mapping(self):
        frame = volmdlr.Frame2D(volmdlr.O2D, volmdlr.Vector2D(1, 1), volmdlr.Vector2D(1, -1))
        side = 'new'
        mapped_composite_2d = self.composite_2d.frame_mapping(frame, side)

        for p1, p2 in zip(mapped_composite_2d.primitives, self.primitives):
            self.assertNotEqual(p1, p2)
            p2 = p2.frame_mapping(frame, side)
            self.assertEqual(p1, p2)


if __name__ == "__main__":
    unittest.main()
