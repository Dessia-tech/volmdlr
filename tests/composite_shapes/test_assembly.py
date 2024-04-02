import unittest
import volmdlr
import volmdlr.primitives3d as p3d
import volmdlr.composite_shapes


class TestAssembly(unittest.TestCase):
    def test_bounding_box(self):
        box1 = p3d.Block(volmdlr.OXYZ)
        components = [box1, box1]
        positions = [volmdlr.OXYZ, volmdlr.Frame3D(volmdlr.Point3D(0, 0, 1), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)]
        assembly = volmdlr.composite_shapes.Assembly(components, positions)
        self.assertAlmostEqual(assembly.bounding_box.volume(), 2.0)  # add assertion here


if __name__ == '__main__':
    unittest.main()
