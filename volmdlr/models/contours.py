#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 2022

@author: s.bendjebla
"""

# %% Librairies

import volmdlr.wires
import volmdlr.edges

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

Contour2d_1 = volmdlr.wires.Contour2D(primitives)
