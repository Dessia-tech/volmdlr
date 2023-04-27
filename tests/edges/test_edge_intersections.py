import unittest
from itertools import product

import volmdlr
from volmdlr.models.edges import bspline1, lineseg, arc, arc_ellipse2d


class TestEdge2DIntersections(unittest.TestCase):

    def test_edge_intersections(self):
        expected_results = [[volmdlr.Point2D(0.21547678763159617, 0.17126976164912053),
                             volmdlr.Point2D(1.4999999999999996, 5.551115123125783e-17),
                             volmdlr.Point2D(2.7845232123683274, -0.17126976164911029)],
                            [volmdlr.Point2D(0.13933375795115266, 0.120520844200026),
                             volmdlr.Point2D(1.768612407725106, -0.1299992525199197)],
                            [volmdlr.Point2D(0.8296845207584518, 0.2682270371882852)],
                            [volmdlr.Point2D(0.21547678763159617, 0.17126976164912053),
                             volmdlr.Point2D(1.4999999999999996, 5.551115123125783e-17),
                             volmdlr.Point2D(2.7845232123683274, -0.17126976164911029)],
                            [volmdlr.Point2D(1.8893801268948387, -0.051917350252645156),
                             volmdlr.Point2D(0.07997689764965771, 0.18933641364671233)],
                            [volmdlr.Point2D(0.7598562777382369, 0.0986858296349018)],
                            [volmdlr.Point2D(0.13933375795115266, 0.120520844200026),
                             volmdlr.Point2D(1.768612407725106, -0.1299992525199197)],
                            [volmdlr.Point2D(1.8893801268948387, -0.051917350252645156),
                             volmdlr.Point2D(0.07997689764965771, 0.18933641364671233)],
                            [volmdlr.Point2D(0.5949972004520672, -0.19981456124516805)],
                            [volmdlr.Point2D(0.8296845207584518, 0.2682270371882852)],
                            [volmdlr.Point2D(0.7598562777382369, 0.0986858296349018)],
                            [volmdlr.Point2D(0.5949972004520672, -0.19981456124516805)]]

        intersection_results = []
        for edge1, edge2 in product([bspline1, lineseg, arc, arc_ellipse2d], repeat=2):
            if edge1 == edge2:
                continue
            intersections = edge1.intersections(edge2)
            intersection_results.append(intersections)
        for intersections, expected_intersections in zip(intersection_results, expected_results):
            for intersection, expected_intersection in zip(intersections, expected_intersections):
                self.assertTrue(intersection.is_close(expected_intersection))


if __name__ == '__main__':
    unittest.main()
