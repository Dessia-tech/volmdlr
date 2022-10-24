import unittest
from volmdlr import primitives3d as p3d
import volmdlr as vm
import volmdlr.core
import volmdlr.wires as vmw
import volmdlr.edges as vme

circle = vmw.Circle2D(vm.O2D, 0.3)
frame1 = volmdlr.Frame3D(vm.O3D, vm.X3D, vm.Y3D, vm.Z3D)
frame2 = volmdlr.Frame3D(vm.Point3D(2, 2, 2), vm.X3D, vm.Y3D, vm.Z3D)
frame3 = volmdlr.Frame3D(vm.O3D, vm.X3D, vm.Z3D, vm.Y3D)

points = [vm.Point2D(0.5, 0.5), vm.Point2D(-0.5, 0.5), vm.Point2D(-0.5, -0.5), vm.Point2D(0.5, -0.5)]
polygon = vmw.ClosedPolygon2D(points, "square")

ext1 = p3d.ExtrudedProfile(frame1, circle, [], 2)
ext2 = p3d.ExtrudedProfile(frame1, polygon, [], 1)

ext3 = p3d.ExtrudedProfile(frame2, circle, [], 2)

ext4 = p3d.ExtrudedProfile(frame3, polygon, [circle], -1)


class TestExtrudedProfile(unittest.TestCase):
    def test_number_faces(self):
        self.assertEqual(len(ext1.faces), 3)
        self.assertEqual(len(ext2.faces), 6)

    def test_extrusion_vector(self):
        self.assertEqual(ext1.extrusion_vector.z, 2)
        self.assertEqual(ext2.extrusion_vector.z, 1)

    def test_origin(self):
        self.assertEqual(ext3.plane_origin.x, 2)
        self.assertEqual(ext3.plane_origin.y, 2)
        self.assertEqual(ext3.plane_origin.y, 2)

    def test_orientation(self):
        self.assertEqual(ext4.y.z, 1)
        self.assertEqual(ext4.extrusion_vector.y, -1)

if __name__ == '__main__':
    unittest.main()
