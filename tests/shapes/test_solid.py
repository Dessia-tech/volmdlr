import math
import unittest
import volmdlr
from volmdlr import faces, surfaces, shapes, wires, curves


class TestSolid(unittest.TestCase):
    faces_list = []
    ax = None
    for vector in [volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D]:
        for direction in [1, -1]:
            normal = direction * vector
            center = direction * vector.to_point()
            plane = surfaces.Plane3D.from_normal(center, normal)
            face = faces.PlaneFace3D.from_surface_rectangular_cut(plane, -1, 1, -1, 1)
            faces_list.append(face)

    solid1 = shapes.Solid.make_solid(shapes.Shell(faces=faces_list))
    faces2 = [f.translation(volmdlr.Vector3D(1, 1, 1)) for f in faces_list]
    solid2 = shapes.Solid.make_solid(shapes.Shell(faces=faces2))

    def test_union(self):
        union = self.solid1.union(self.solid2)
        self.assertAlmostEqual(union.volume(), 15)

    def test_subtraction(self):
        subtraction = self.solid1.subtraction(self.solid2)
        self.assertAlmostEqual(subtraction.volume(), 7.0)

    def test_intersection(self):
        intersection = self.solid1.intersection(self.solid2)
        self.assertAlmostEqual(intersection.volume(), 1.0)
    
    def test_make_extrusion(self):
        length, width, height, radius = 0.4, 0.3, 0.08, 0.1
        outer_contour2d = wires.Contour2D.rectangle_from_center_and_sides(volmdlr.O2D, x_length=length, y_length=width,
                                                                          is_trigo=True)
        inner_contours2d = [wires.Contour2D.from_circle(
            circle=curves.Circle2D.from_center_and_radius(volmdlr.O2D, radius, is_trigo=False))]
        solid = shapes.Solid.make_extrusion_from_frame_and_wires(volmdlr.OXYZ, outer_contour2d,
                                                                 inner_contours2d, height)
        self.assertAlmostEqual(solid.volume(), (length * width - math.pi * (radius ** 2)) * height)

    def test_make_wedge(self):
        dx, dy, dz = 1, 2, 1
        solid = shapes.Solid.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 2, xmax=dx / 2, zmin=dz / 2, zmax=dz / 2,
                                        local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
                                        local_frame_direction=-volmdlr.Y3D,
                                        local_frame_x_direction=volmdlr.X3D)
        
        self.assertAlmostEqual(solid.volume(), (1 / 3) * dy * (1 + 0.5**2 + 0.5))


if __name__ == '__main__':
    unittest.main()
