import unittest
import volmdlr
import volmdlr.edges
from volmdlr.core import CompositePrimitive3D


class TestCompositePrimitive3D(unittest.TestCase):
    def setUp(self):
        self.primitives = [volmdlr.edges.LineSegment3D(volmdlr.O3D, volmdlr.Point3D(4, 2, 1))]
        self.composite_3d = CompositePrimitive3D(self.primitives, name="test")

    def test_plot(self):
        ax = self.composite_3d.plot()

        for ls, line in zip(self.composite_3d.primitives, ax.lines):
            data = line.get_data_3d()

            for i in range(3):
                self.assertListEqual(data[i].tolist(), [ls.start[i], ls.end[i]])
