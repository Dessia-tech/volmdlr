"""
Unittest for ClosedPolygon2D class
"""
import unittest

import volmdlr
import volmdlr.display as vmd
import volmdlr.wires as vmw


class TestClosedPolygon2D(unittest.TestCase):
    def test_triangulation(self):
        # Create a ClosedPolygon2D object with a list of points
        points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)]
        polygon = vmw.ClosedPolygon2D(points)

        # Call the triangulation method with different options
        mesh1 = polygon.triangulation()
        mesh2 = polygon.triangulation("p")
        mesh3 = polygon.triangulation("pa0.25")  # No triangles with area greter than 0.25

        # Assert that the returned object is a vmd.DisplayMesh2D
        self.assertIsInstance(mesh1, vmd.DisplayMesh2D)
        self.assertIsInstance(mesh2, vmd.DisplayMesh2D)
        self.assertIsInstance(mesh3, vmd.DisplayMesh2D)

        # Assert that the mesh has the correct number of points and triangles
        self.assertEqual(len(mesh1.points), 4)
        self.assertEqual(len(mesh1.triangles), 2)
        self.assertEqual(len(mesh2.points), 4)
        self.assertEqual(len(mesh2.triangles), 2)
        self.assertEqual(len(mesh3.points), 5)
        self.assertEqual(len(mesh3.triangles), 4)

    def test_point_belongs(self):
        # create a star-shaped polygon with 5 points
        polygon = vmw.ClosedPolygon2D(
            [
                volmdlr.Point2D(0, 0),
                volmdlr.Point2D(1, 2),
                volmdlr.Point2D(3, 0),
                volmdlr.Point2D(2, -2),
                volmdlr.Point2D(0, -1),
            ]
        )

        # test points inside the polygon
        points_inside = [volmdlr.Point2D(1, 1), volmdlr.Point2D(1.5, 0), volmdlr.Point2D(2, -1)]
        for point in points_inside:
            self.assertTrue(polygon.point_belongs(point))
            self.assertTrue(polygon.point_belongs(point, include_edge_points=True))

        # test points on the edge of the polygon
        edge_points = [
            volmdlr.Point2D(0, 0),
            volmdlr.Point2D(0.5, 1),
            volmdlr.Point2D(1, 2),
            volmdlr.Point2D(2, 1),
            volmdlr.Point2D(3, 0),
            volmdlr.Point2D(2.5, -1),
            volmdlr.Point2D(2, -2),
            volmdlr.Point2D(1, -1.5),
            volmdlr.Point2D(0, -1),
            volmdlr.Point2D(0, -0.5),
        ]
        for point in edge_points:
            self.assertFalse(polygon.point_belongs(point))
            self.assertTrue(polygon.point_belongs(point, include_edge_points=True))

        # test a point outside the polygon
        point = volmdlr.Point2D(4, 0)
        self.assertFalse(polygon.point_belongs(point))
        self.assertFalse(polygon.point_belongs(point, include_edge_points=True))


if __name__ == "__main__":
    unittest.main()
