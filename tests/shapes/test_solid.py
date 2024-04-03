import os
import math
import unittest

from dessia_common.core import DessiaObject

from OCP.TopoDS import TopoDS_Solid

import volmdlr
from volmdlr.models.contours import rim_contour, w, wb, R, Rb, th
from volmdlr import curves, faces, shapes, surfaces, wires, primitives3d


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)))
objects_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "shapes_objects")


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

    solid1 = shapes.Solid.make_solid(shapes.Shell.from_faces(faces=faces_list))
    faces2 = [f.translation(volmdlr.Vector3D(1, 1, 1)) for f in faces_list]
    solid2 = shapes.Solid.make_solid(shapes.Shell.from_faces(faces=faces2))

    def test_check_platform(self):
        self.assertIsNone(self.solid1._check_platform())

    def test_to_dict_dict_to_object(self):
        to_dict = self.solid1.to_dict()
        dict_to_obejct = DessiaObject.dict_to_object(to_dict)
        self.assertEqual(dict_to_obejct, self.solid1)

    def test_to_brep_from_brep(self):
        self.solid1.to_brep(objects_folder + "/test_to_brep.brep")
        from_brep = shapes.Solid.from_brep(objects_folder + "/test_to_brep.brep")
        self.assertEqual(from_brep, self.solid1)

    def test_distance(self):
        faces2 = [f.translation(volmdlr.Vector3D(2, 3, 4)) for f in self.faces_list]
        shell = shapes.Shell.from_faces(faces=faces2)
        solid2 = shapes.Solid.make_solid(shell)
        distance = self.solid1.distance(solid2)
        self.assertAlmostEqual(distance, 2.23606797749979)
        point1, point2 = self.solid1.distance_points(solid2)
        self.assertEqual(point2, volmdlr.Point3D(1.0, 2.0, 3.0))
        self.assertEqual(point1, volmdlr.Point3D(1.0, 1.0, 1.0))

    def test_union(self):
        union = self.solid1.union(self.solid2)[0]
        self.assertAlmostEqual(union.volume(), 15)

    def test_subtraction(self):
        subtraction = self.solid1.subtraction(self.solid2)[0]
        self.assertAlmostEqual(subtraction.volume(), 7.0)

    def test_intersection(self):
        intersection = self.solid1.intersection(self.solid2)[0]
        self.assertAlmostEqual(intersection.volume(), 1.0)

    def test_box(self):
        box = shapes.Solid.make_box(length=2, width=3, height=5)
        self.assertEqual(box.volume(), 2 * 3 * 5)

    def test_cone(self):
        cone = shapes.Solid.make_cone(radius1=0, radius2=5, height=5, direction=volmdlr.X3D, angle_degrees=270)
        self.assertAlmostEqual(cone.volume(), 98.1747704332957)

    def test_cylinder(self):
        cylinder = shapes.Solid.make_cylinder(radius=5, height=5, direction=volmdlr.X3D, angle_degrees=270)
        self.assertAlmostEqual(cylinder.volume(), 294.5243112989299)

    def test_sphere(self):
        sphere1 = shapes.Solid.make_sphere(radius=5, direction=volmdlr.X3D, angle_degrees1=-90,
                                           angle_degrees2=60, angle_degrees3=270)
        self.assertAlmostEqual(sphere1.volume(), 387.6486925590814)

    def test_torus(self):
        torus1 = shapes.Solid.make_torus(radius1=2, radius2=.5, direction=volmdlr.X3D, angle_degrees1=0,
                                         angle_degrees2=360)
        self.assertAlmostEqual(torus1.volume(), 9.869604401089358)

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

        self.assertAlmostEqual(solid.volume(), (1 / 3) * dy)

        solid = shapes.Solid.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 4, xmax=3 * dx / 4,
                                        zmin=dz / 4, zmax=3 * dz / 4,
                                        local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
                                        local_frame_direction=-volmdlr.Y3D,
                                        local_frame_x_direction=volmdlr.X3D)

        self.assertAlmostEqual(solid.volume(), (1 / 3) * dy * (1 + 0.5 ** 2 + 0.5))

    def test_sweep(self):
        point1, point2 = volmdlr.Point3D(1.0, 1.0, 0.0), volmdlr.Point3D(1.0, 0.5, 0.0)
        point3 = volmdlr.Point3D(0.5, 0.5, 0.0)
        path = primitives3d.OpenRoundedLineSegments3D([point1, point2, point3], {"1": 0.2})
        outer_contour = wires.Contour2D.from_circle(curves.Circle2D(volmdlr.OXY, 0.05))
        inner_contours = [wires.Contour2D.from_circle(curves.Circle2D(volmdlr.OXY, 0.045))]
        direction = (point2 - point1).unit_vector()
        frame = volmdlr.Frame3D.from_point_and_vector(point=point1, vector=direction, main_axis=volmdlr.Z3D)
        section = faces.PlaneFace3D(surface3d=surfaces.Plane3D(frame=frame),
                                    surface2d=surfaces.Surface2D(outer_contour=outer_contour, inner_contours=inner_contours))
        sweep = shapes.Solid.make_sweep_from_contour(outer_contour, path, inner_contours)
        self.assertEqual(len(sweep.primitives[0].primitives), 8)

        path = wires.Wire3D.from_points([point1, point2, point3])
        sweep1 = shapes.Solid.make_sweep(face=section, path=path, transition_mode="right")
        self.assertEqual(len(sweep1.primitives[0].primitives), 6)
        sweep2 = shapes.Solid.make_sweep(face=section, path=path, transition_mode="round")
        self.assertEqual(len(sweep2.primitives[0].primitives), 10)

    def test_make_revolve(self):
        inner_contours = [wires.Contour2D.from_circle(
            curves.Circle2D.from_center_and_radius(volmdlr.Point2D(-0.5 * (w - wb), Rb - 0.15 * (Rb - (R - th))),
                                                   radius=0.5 * th))]
        y = volmdlr.X3D.random_unit_normal_vector()
        z = volmdlr.X3D.cross(y)
        axis_point = 0.5 * volmdlr.X3D.to_point()
        frame = volmdlr.Frame3D(axis_point, volmdlr.X3D, z, y)
        revolution_shape = shapes.Solid.make_revolve_from_contour(frame=frame, contour2d=rim_contour,
                                                                  axis_point=axis_point, axis=volmdlr.X3D,
                                                                  inner_contours=inner_contours,
                                                                  angle=3.1415, name="Conical rim")
        self.assertEqual(revolution_shape.name, "Conical rim")
        self.assertEqual(len(revolution_shape.primitives[0].primitives), 11)
        self.assertIsInstance(revolution_shape.wrapped, TopoDS_Solid)

    def test_loft(self):
        diameter = 0.3
        circle1 = curves.Circle3D(frame=volmdlr.OXYZ, radius=diameter / 2)
        circle2 = curves.Circle3D(
            frame=volmdlr.Frame3D(volmdlr.Point3D(0.3, 0.0, 0.5), volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D),
            radius=circle1.radius / 2)
        circle3 = curves.Circle3D(
            frame=volmdlr.Frame3D(volmdlr.Point3D(0.6, 0.0, 0.3), volmdlr.Y3D, volmdlr.X3D, -volmdlr.Z3D),
            radius=circle1.radius * 0.6)
        sections = [wires.Contour3D.from_circle(circle1), wires.Contour3D.from_circle(circle2),
                    wires.Contour3D.from_circle(circle3)]
        loft = shapes.Solid.make_loft(sections=sections, name="loft")
        self.assertEqual(loft.name, "loft")
        self.assertEqual(len(loft.primitives[0].primitives), 4)
        self.assertIsInstance(loft.wrapped, TopoDS_Solid)

        section2 = wires.Contour3D.from_points([volmdlr.Point3D(0.0, -0.1, 0.3), volmdlr.Point3D(0.05, -0.05, 0.3),
                                                volmdlr.Point3D(0.05, 0.05, 0.3), volmdlr.Point3D(0.0, 0.1, 0.3),
                                                volmdlr.Point3D(-0.05, 0.05, 0.3), volmdlr.Point3D(-0.05, -0.05, 0.3)])

        sections = [wires.Contour3D.from_circle(circle1), section2, volmdlr.Point3D(0.0, 0.0, 0.45)]
        loft = shapes.Solid.make_loft(sections=sections, ruled=True)
        self.assertEqual(len(loft.primitives[0].primitives), 15)
        self.assertIsInstance(loft.wrapped, TopoDS_Solid)


if __name__ == '__main__':
    unittest.main()
