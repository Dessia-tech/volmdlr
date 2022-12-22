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
        mesh2 = polygon.triangulation('p')
        mesh3 = polygon.triangulation('pa0.25')

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


if __name__ == '__main__':
    unittest.main()