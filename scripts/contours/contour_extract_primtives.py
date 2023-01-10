#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 2022

@author: s.bendjebla
"""

# %% Librairies

import volmdlr.wires
import volmdlr.edges
import matplotlib.pyplot as plt

# %% Data

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
contour2d = volmdlr.wires.Contour2D(primitives)

point1 = volmdlr.Point2D(-0.007116025403784438, 0.010325317547305484)
point2 = volmdlr.Point2D(-0.005383974596215561, 0.011325317547305485)

# %% Wire.extract_without_primitives

extracted_primitives_inside_true = contour2d.extract_primitives(point1, contour2d.primitives[3],
                                                                point2, contour2d.primitives[3],
                                                                inside=True)

extracted_primitives_inside_false = contour2d.extract_primitives(point1, contour2d.primitives[3],
                                                                 point2, contour2d.primitives[3],
                                                                 inside=False)

extracted_primitives = [extracted_primitives_inside_true, extracted_primitives_inside_false]

# %% Plots

fig, axs = plt.subplots(1, 2)

titles = ["inside=True", "inside=False"]
colors = ['g', 'r']
for i in range(len(axs)):
    contour2d.plot(ax=axs[i])
    point1.plot(ax=axs[i])
    point2.plot(ax=axs[i])
    for prim in extracted_primitives[i]:
        prim.plot(ax=axs[i], color=colors[i])
    axs[i].set_title(titles[i])
