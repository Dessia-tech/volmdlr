#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 2022.

@author: s.bendjebla
"""

# %% Libraries

import volmdlr.edges
import volmdlr.wires

# %% Contour2d_1

primitives = [
    volmdlr.edges.LineSegment2D(volmdlr.Point2D(0.001, 0.014),
                                volmdlr.Point2D(0.001, 0.0125)),

    volmdlr.edges.Arc2D(volmdlr.Point2D(0.001, 0.0125), 
                        volmdlr.Point2D(0.009862829911410362, 0.007744326060968065),
                        volmdlr.Point2D(0.012539936203984454, 0.0)),

    volmdlr.edges.Arc2D(volmdlr.Point2D(0.012539936203984454, 0.0), 
                        volmdlr.Point2D(0.0, -0.012539936203984454),
                        volmdlr.Point2D(-0.012539936203984454, 0.0)),

    volmdlr.edges.Arc2D(volmdlr.Point2D(-0.012539936203984454, 0.0), 
                        volmdlr.Point2D(-0.00921384654213387, 0.008506176103162205),
                        volmdlr.Point2D(-0.001, 0.0125)),

    volmdlr.edges.LineSegment2D(volmdlr.Point2D(-0.001, 0.0125),
                                volmdlr.Point2D(-0.001, 0.014)),

    volmdlr.edges.LineSegment2D(volmdlr.Point2D(-0.001, 0.014),
                                volmdlr.Point2D(0.001, 0.014))
    ]

contour2d_1 = volmdlr.wires.Contour2D(primitives)

# %% Contour2d_2

points2d = [volmdlr.Point2D(0.7557261768236382, 0.047840610095316816),
            volmdlr.Point2D(0.7557153778090784, 0.053529840943971466),
            volmdlr.Point2D(0.7548352148318592, 0.059379242604894106),
            volmdlr.Point2D(0.7509972673886837, 0.06812701209830616),
            volmdlr.Point2D(0.7463815515547884, 0.06876367708780329),
            volmdlr.Point2D(0.7426873841819802, 0.06850598841747624),
            volmdlr.Point2D(0.7386859822441963, 0.0690723949917957),
            volmdlr.Point2D(0.7344087842895368, 0.0660771254673235),
            volmdlr.Point2D(0.7316087522378829, 0.0572133737658937),
            volmdlr.Point2D(0.7306755405178993, 0.04995246337945117),
            volmdlr.Point2D(0.7310451707836978, 0.044346281865748856)]

curve2d_1 = volmdlr.edges.BSplineCurve2D(degree = 3,
                                         control_points = points2d,
                                         knot_multiplicities = [4, 1, 1, 1, 1, 1, 1, 1, 4],
                                         knots = [0.0, 0.26206414743620277, 0.3531225286567413,
                                                  0.42241139891509305, 0.4822254408401885,
                                                  0.5440692281614715, 0.6202023477514225,
                                                  0.7222748100879735, 1.0])

points2d = [volmdlr.Point2D(0.7310451707836978, 0.044346281865748856),
            volmdlr.Point2D(0.7308827017192874, 0.03840965331596735),
            volmdlr.Point2D(0.7318250470334302, 0.03228668394110015),
            volmdlr.Point2D(0.7354633489443133, 0.023006441447214778),
            volmdlr.Point2D(0.7402551819180077, 0.022102466430262133),
            volmdlr.Point2D(0.7440795569187786, 0.02242177670352842),
            volmdlr.Point2D(0.7482500898246996, 0.02178574616551938),
            volmdlr.Point2D(0.7525849938211369, 0.025193737559804485),
            volmdlr.Point2D(0.7553036624157348, 0.03437896214237193),
            volmdlr.Point2D(0.7561315600530795, 0.042073178367651974),
            volmdlr.Point2D(0.7557261768236382, 0.047840610095316816)]

curve2d_2 = volmdlr.edges.BSplineCurve2D(degree = 3,
                                         control_points = points2d,
                                         knot_multiplicities = [4, 1, 1, 1, 1, 1, 1, 1, 4],
                                         knots = [0.0, 0.26231860269644847, 0.35446984209379995,
                                                  0.4237193471508679, 0.48275286541966095,
                                                  0.5439383602451768, 0.62045686469776,
                                                  0.7233434357617774, 1.0])

contour2d_2 = volmdlr.wires.Contour2D([curve2d_1, curve2d_2])
