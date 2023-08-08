import unittest
from itertools import product

import volmdlr
import volmdlr.edges as vme
import volmdlr.models.edges as edges_models


class TestEdgesDistances(unittest.TestCase):
    def test_edges_distances_2d(self):
        expected_distances = [3.0, 3.1719292821844314, 1.413480226519909, 2.9154759474226504, 3.144603017460846,
                              3.293419303698591, 1.4282856857085704, 2.8735269769962803, 3.1832796902406435,
                              3.3501380577938424, 1.4456832294800965, 2.857575183838226, 1.000566410613742,
                              1.1850960991514135, 0.0, 1.2598108740840506]

        vector1 = volmdlr.Vector2D(1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        distances = []
        for edge1, edge2 in product([edges_models.bspline1, edges_models.lineseg,
                                     edges_models.arc, edges_models.arc_ellipse2d], repeat=2):
            edge2 = edge2.frame_mapping(volmdlr.Frame2D(volmdlr.Point2D(3, 3), vector1, vector2), 'new')
            dist, min_dist_point_arc, min_dist_point_lineseg = edge1.minimum_distance(edge2, True)
            distances.append(dist)
        for distance, expected_distance in zip(distances, expected_distances):
            self.assertAlmostEqual(distance, expected_distance)

    def test_edges_distances_3d(self):
        expected_distances = [3.2969896808200936, 0.6599120175960898, 0.019789865983643993, 0.6385899371099841,
                              0.6599120175960898, 6.060992532937859, 2.667416425120059, 1.2248182129344367,
                              0.019789865983643993, 2.667416425120059, 5.08435180389005, 0.8731527317520039,
                              0.638417217163073, 1.2247754640880462, 0.8731527317520039, 4.289628496553127]
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)

        distances = []
        for edge1, edge2 in product([edges_models.bspline_curve3d(), edges_models.linesegment3d(),
                                     edges_models.arc3d(), edges_models.arc_ellipse3d()], repeat=2):
            if edge1 == edge2:
                edge2 = edge2.frame_mapping(volmdlr.Frame3D(volmdlr.Point3D(2, 5, -3), vector1, vector2, vector3),
                                            'new')
            dist, min_dist_point_arc, min_dist_point_lineseg = edge1.minimum_distance(edge2, True)
            distances.append(dist)
        for distance, expected_distance in zip(distances, expected_distances):
            self.assertAlmostEqual(distance, expected_distance)


if __name__ == '__main__':
    unittest.main()
