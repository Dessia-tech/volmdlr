"""
Tests for places faces
"""
import math
import os
import unittest

import dessia_common.core
import dessia_common.core as dc

import volmdlr
from volmdlr import edges, wires, faces, surfaces, curves


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "objects_planeface_tests")


class TestPlaneFace3D(unittest.TestCase):
    face_with_3holes = dc.DessiaObject.from_json(os.path.join(folder, 'face_with_3holes.json'))
    face = dc.DessiaObject.from_json(os.path.join(folder,
                                                       'face_to_cut_the_one_with_3holes.json'))
    plane_face_cylindricalface_intersec = dc.DessiaObject.from_json(
        os.path.join(folder, 'plane_face_cylindrical_face_intersec.json'))

    def test_area(self):
        self.assertAlmostEqual(self.face_with_3holes.area(), 0.12160000)

    def test_point_belongs(self):
        volume_model = dc.DessiaObject.from_json(os.path.join(folder, 'test_planeface_point_belongs.json'))
        self.assertTrue(volume_model.primitives[0].point_belongs(volume_model.primitives[1]))

    def test_plane_face_intersection(self):
        face1, face2 = dc.DessiaObject.from_json(os.path.join(folder, 'test_planeface_intersections.json')).primitives
        face_intersections = face1.face_intersections(face2)

        self.assertEqual(1, len(face_intersections))
        self.assertAlmostEqual(0.003600000000881293, face_intersections[0].length())
        face1, face2 = faces.PlaneFace3D.from_json(
            os.path.join(folder, 'test_planef_inters291123.json')).primitives
        face_intersections = face1.face_intersections(face2)
        line_seg = edges.LineSegment3D(
            volmdlr.Point3D(0.034031786272172244, 0.019018077708099396, 0.07196336225667499),
            volmdlr.Point3D(0.03756212352555941, 0.022324611455208112, 0.07196336225667499)
        )
        self.assertTrue(face_intersections[0].primitives[0], line_seg)

    def test_face_inside(self):
        face2 = self.face.frame_mapping(volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.Vector3D(0.5, 0, 0),
                                        volmdlr.Vector3D(0, 0.5, 0), volmdlr.Vector3D(0, 0, 0.5)), 'old')
        self.assertEqual(self.face.face_inside(face2), True)
        self.assertEqual(face2.face_inside(self.face), False)
        face1, face2 = dc.DessiaObject.from_json(os.path.join(folder, 'test_face_inside.json')).primitives
        self.assertTrue(face1.face_inside(face2))
        face1, face2 = dc.DessiaObject.from_json(os.path.join(folder, 'test_face3_face_inside.json')).primitives
        self.assertFalse(face1.face_inside(face2))

    def test_face_intersections_with_holes(self):
        face_intersections = self.face.face_intersections(self.face_with_3holes)
        self.assertEqual(len(face_intersections), 4)

    def test_line_intersections(self):
        line_inside_hole = curves.Line3D(volmdlr.Point3D(0.1, 0.0, -.3), volmdlr.Point3D(-0.1, 0.0, 0.3))
        line_inside_face = curves.Line3D(volmdlr.Point3D(-0.05, 0.0, -.3), volmdlr.Point3D(-0.1, 0.0, 0.3))
        self.assertEqual([], self.face_with_3holes.line_intersections(line_inside_hole))
        self.assertEqual(1, len(self.face_with_3holes.line_intersections(line_inside_face)))

    def test_divide_face(self):
        face_intersections = self.face.face_intersections(self.face_with_3holes)
        cutting_contours = self.face_with_3holes.get_face_cutting_contours(
            {self.face_with_3holes: face_intersections})
        new_faces = self.face_with_3holes.divide_face(cutting_contours)
        self.assertEqual(len(new_faces), 2)
        cutting_contour = wires.Wire2D.from_points([
            volmdlr.Point2D(0.5, 1.), volmdlr.Point2D(0.5, 0.75),
            volmdlr.Point2D(1, 0.75), volmdlr.Point2D(1, 1), volmdlr.Point2D(1.25, 1),
            volmdlr.Point2D(1.25, 1.5), volmdlr.Point2D(1, 1.5)
        ])
        face_tobe_divided = dc.DessiaObject.from_json(os.path.join(folder, 'face_tobe_divided.json'))
        divided_faces = face_tobe_divided.divide_face([cutting_contour])
        self.assertEqual(len(divided_faces), 4)
        expected_areas = [1.4320458460875176, 0.05704584608751772, 0.125, 0.125]
        for i, face in enumerate(divided_faces):
            self.assertAlmostEqual(expected_areas[i], face.area())
        source_folder = os.path.join(folder, 'test_planeface_divide_face_json_files')
        expected_faces_areas = [[0.005578804359415823, 0.0948396741089057, 0.23430978309161565, 0.0055788043593624215],
                                [0.0855613934860544, 0.032085522557644186, 0.01069517418574345],
                                [0.002005345159845676, 0.002005345159820638, 0.0033422419331328506,
                                 0.0006684483866419249], [0.3403070659192998, 0.005578804359415823],
                                [0.010695174185850198, 0.08556139348605463, 0.0427806967433878],
                                [0.07754001284661505, 0.002673793546460905]]
        file_names = ['test_face_divide_face5.json', 'test_face_divide_face2.json',
                      'test_planeface3d_divide_face.json', 'test_face_divide_face3.json', 'test_face_divide_face.json',
                      'test_face_divide_face6.json']
        faces_areas = []
        for filename in file_names:
            file_path = os.path.join(source_folder, filename)
            obj = dc.DessiaObject.from_json(file_path)
            face = obj.primitives[0]
            list_cutting_contours = obj.primitives[1:]
            divide_faces = face.divide_face(list_cutting_contours)
            areas = [face.area() for face in divide_faces]
            faces_areas.append(areas)
        for solution, expected_solution in zip(faces_areas, expected_faces_areas):
            self.assertEqual(len(solution), len(expected_solution))
            for solution_area, expected_solution_area in zip(solution, expected_solution):
                self.assertAlmostEqual(solution_area, expected_solution_area)

    def test_set_operations_new_faces(self):
        volumemodel = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, 'test_set_operations_new_faces.json'))
        plane_face, cutting_contours3d = volumemodel.primitives[0], volumemodel.primitives[1:]
        divide_face = plane_face.set_operations_new_faces({plane_face: cutting_contours3d})
        divide_face = sorted(divide_face, key=lambda face: face.area())
        self.assertEqual(len(divide_face), 7)
        expected_areas = [0.05000000000000002, 0.25, 0.4375, 0.5, 1.0621681469282045, 1.29875, 2.4002627208849607]
        areas = [f.area() for f in divide_face]
        for area, expected_area in zip(areas, expected_areas):
            self.assertAlmostEqual(area, expected_area, 6)

    def test_cylindricalface_intersections(self):
        R = 0.15
        cylindricalsurface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, R)
        face = faces.CylindricalFace3D.from_surface_rectangular_cut(cylindricalsurface, 0, volmdlr.TWO_PI, -.25, .25)
        """ ========== CIRCLE3D ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 2)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(len(face_intersections), 1)
        self.assertIsInstance(face_intersections[0].primitives[0], edges.FullArc3D)
        self.assertEqual(face_intersections[0].primitives[0].circle.center, volmdlr.O3D)
        self.assertEqual(face_intersections[0].primitives[0].circle.radius, 0.15)
        """ ========== FULL ELLIPSE3D ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 4)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(len(face_intersections), 1)
        self.assertIsInstance(face_intersections[0].primitives[0], edges.FullArcEllipse3D)
        self.assertEqual(face_intersections[0].primitives[0].ellipse.center, volmdlr.O3D)
        self.assertAlmostEqual(face_intersections[0].primitives[0].ellipse.major_axis, 0.21213203435596426)
        self.assertTrue(face_intersections[0].primitives[0].ellipse.major_dir.is_close(
            volmdlr.Vector3D(0, -0.7071067811865475, 0.7071067811865475)))
        """ ========== THREE ARC ELLIPSES ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 7)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(len(face_intersections), 3)
        for inter in face_intersections:
            self.assertIsInstance(inter.primitives[0], edges.ArcEllipse3D)
        self.assertEqual(face_intersections[0].primitives[0].ellipse.center, volmdlr.O3D)
        self.assertTrue(face_intersections[0].primitives[0].ellipse.major_dir.is_close(volmdlr.Point3D(
            2.6567716615652136e-17, -0.4338837391180807, 0.9009688679021675)))
        self.assertAlmostEqual(face_intersections[0].primitives[0].ellipse.major_axis, 0.3457147306439571)
        list_expected_points = [[volmdlr.Point3D(0.08947272158306664, -0.12039365470206077, 0.25),
                                 volmdlr.Point3D(0.136637076048, -0.061889493851, 0.128514858204),
                                 volmdlr.Point3D(0.15, -0.0, 0.0)],
                                [volmdlr.Point3D(0.15, -0.0, 0.0),
                                 volmdlr.Point3D(0.136637075473, 0.06188949512, -0.128514860839),
                                 volmdlr.Point3D(0.08947272158306664, 0.12039365470206077, -0.25)],
                                [volmdlr.Point3D(-0.08947272158306664, 0.12039365470206077, -0.25),
                                 volmdlr.Point3D(-0.15, 0, 0),
                                 volmdlr.Point3D(-0.08947272158306664, -0.12039365470206077, 0.25)]]
        for expected_points, wire in zip(list_expected_points, face_intersections):
            arcellipse = wire.primitives[0]
            self.assertTrue(expected_points[0].is_close(arcellipse.start))
            self.assertTrue(expected_points[1].is_close(arcellipse.middle_point()))
            self.assertTrue(expected_points[2].is_close(arcellipse.end))
        """ ========== TWO PARALLEL LINES ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertTrue(face_intersections[0].primitives[0].is_close(
            edges.LineSegment3D(volmdlr.Point3D(0.15, 0.0, -0.25),
                                volmdlr.Point3D(0.15, 0.0, 0.25))))
        self.assertTrue(face_intersections[1].primitives[0].is_close(
            edges.LineSegment3D(volmdlr.Point3D(-0.15, 0.0, -0.25),
                                volmdlr.Point3D(-0.15, 0.0, 0.25))))
        """ ========== ONE LINE ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.translation(R * volmdlr.Y3D)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(face_intersections[0].primitives[0], edges.LineSegment3D(volmdlr.Point3D(0.0, 0.15, -0.25),
                                                                                  volmdlr.Point3D(0.0, 0.15, 0.25)))
        """======= [] ==========="""
        planeface, cylface = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, 'test_planeface_cylindricalface_intersections_none.json')).primitives
        self.assertFalse(planeface.face_intersections(cylface))

    def test_conical_face_intersections(self):
        def get_face(plane, x1=-1, x2=1, y1=-1, y2=1):
            return faces.PlaneFace3D.from_surface_rectangular_cut(plane, x1, x2, y1, y2)
        conical_surface = surfaces.ConicalSurface3D(volmdlr.OXYZ, math.pi / 6)
        conical_face = faces.ConicalFace3D.from_surface_rectangular_cut(conical_surface,
                                                                        0, 1.5 * math.pi, 0., 1)
        """======== Arc3D ========="""
        plane1 = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0.5),
                                                  volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D))
        face = get_face(plane1)
        intersections = face.face_intersections(conical_face)
        self.assertEqual(len(intersections), 1)
        self.assertIsInstance(intersections[0].primitives[0], edges.Arc3D)
        self.assertAlmostEqual(intersections[0].length(), 1.360349523175663)
        """========= Hyperbola (BSpline3D) ======="""
        plane2 = surfaces.Plane3D(
            volmdlr.Frame3D(volmdlr.Point3D(0, 0.25, 0.5), volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D))
        face = get_face(plane2)
        intersections = face.face_intersections(conical_face)
        self.assertEqual(len(intersections), 1)
        self.assertIsInstance(intersections[0].primitives[0], edges.BSplineCurve3D)
        self.assertAlmostEqual(intersections[0].length(), 1.5797706721593943)

        """============== Two LineSegments ==================="""
        plane3 = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(0, 0.0, 0.5),
                                                  volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D))
        face = get_face(plane3)
        intersections = face.face_intersections(conical_face)
        self.assertEqual(len(intersections), 2)
        self.assertIsInstance(intersections[0].primitives[0], edges.LineSegment3D)
        self.assertAlmostEqual(intersections[0].length(), 1.1547005383794386)
        self.assertAlmostEqual(intersections[1].length(), 1.1547005383794386)
        """===================== Ellipse 3D ============================="""
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0.5), vector1, vector2, vector3)
        plane4 = surfaces.Plane3D(frame)
        face = get_face(plane4)
        intersections = face.face_intersections(conical_face)
        self.assertEqual(len(intersections), 2)
        self.assertIsInstance(intersections[0].primitives[0], edges.ArcEllipse3D)
        self.assertIsInstance(intersections[1].primitives[0], edges.ArcEllipse3D)
        self.assertAlmostEqual(intersections[0].length(), 0.7233023692399578)
        self.assertAlmostEqual(intersections[1].length(), 1.1339664625179093, 6)
        """================== Parabola ===================="""
        point1 = conical_surface.frame.origin
        point2 = conical_surface.frame.local_to_global_coordinates(
            volmdlr.Point3D(10 * math.tan(conical_surface.semi_angle), 0, 10))
        generatrix = edges.LineSegment3D(point1, point2)
        normal = generatrix.unit_normal_vector()
        plane_origin = frame.origin - normal * .5
        plane5 = surfaces.Plane3D.from_normal(plane_origin, normal)
        face = get_face(plane5)
        intersections = face.face_intersections(conical_face)
        self.assertEqual(len(intersections), 1)
        self.assertIsInstance(intersections[0].primitives[0], edges.BSplineCurve3D)
        self.assertAlmostEqual(intersections[0].length(), 1.19097843887217)

    def test_linesegment_inside(self):
        lineseg = volmdlr.edges.LineSegment3D(volmdlr.Point3D(0.2, 0, -0.2), volmdlr.Point3D(0.1, 0.0, 0.2))
        self.assertTrue(self.plane_face_cylindricalface_intersec.linesegment_inside(lineseg))
        lineseg1 = volmdlr.edges.LineSegment3D(volmdlr.Point3D(0.2, 0, -0.2), volmdlr.Point3D(0.1, 0.1, 0.2))
        self.assertFalse(self.plane_face_cylindricalface_intersec.linesegment_inside(lineseg1))

    def test_circle_inside(self):
        circle = curves.Circle3D(volmdlr.OZXY, 0.1)
        self.assertTrue(self.plane_face_cylindricalface_intersec.circle_inside(circle))
        circle2 = curves.Circle3D(volmdlr.OYZX, 0.1)
        self.assertFalse(self.plane_face_cylindricalface_intersec.circle_inside(circle2))

    def test_merges_faces(self):
        source_folder = os.path.join(folder, 'test_planeface3d_merge_faces_json_files')
        faces_areas = []
        file_names = ['test_merge_faces4.json', 'test_merge_faces5.json', 'faces_merge_faces2.json',
                      'faces_merge_faces3.json', 'faces_merge_faces4.json']
        for filename in file_names:
            file_path = os.path.join(source_folder, filename)
            obj = dc.DessiaObject.from_json(file_path)
            faces_ = obj.primitives
            merged_faces = faces.PlaneFace3D.merge_faces(faces_)
            areas = []
            for face in merged_faces:
                areas.append(face.area())
            faces_areas.append(areas)
        expected_faces_areas = [[0.1621764423452034], [0.15508002569387766],
                                [0.005347587092921799, 0.032085522557310564, 0.18181796115851334],
                                [0.08021380639307588],
                                [0.05578804359331602, 0.005578804359362449, 0.011157608718620371,
                                 0.022315217437398727, 0.05020923923410867]]
        for solution, expected_solution in zip(faces_areas, expected_faces_areas):
            self.assertEqual(len(solution), len(expected_solution))
            for solution_area, expected_solution_area in zip(solution, expected_solution):
                self.assertAlmostEqual(solution_area, expected_solution_area)

    def test_grid_points(self):
        surface3d = surfaces.Plane3D(volmdlr.OXYZ)
        outer_contour2d = wires.Contour2D.from_circle(curves.Circle2D.from_center_and_radius(volmdlr.O2D, 1.0))
        inner_contour2d = wires.Contour2D.from_circle(curves.Circle2D.from_center_and_radius(volmdlr.O2D, 0.25,
                                                                                             is_trigo=False))
        surface2d = surfaces.Surface2D(outer_contour2d, [inner_contour2d])
        face = faces.PlaneFace3D(surface3d, surface2d)
        grid_points = face.grid_points([10, 10])
        self.assertEqual(len(grid_points), 56)


if __name__ == '__main__':
    unittest.main()
