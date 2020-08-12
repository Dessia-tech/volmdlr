#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 09:53:05 2020

@author: Pierrem
"""
from volmdlr import plot_data
import volmdlr as vm

hatching = plot_data.HatchingSet(1)
plot_data_state = plot_data.PlotDataState(name='name', hatching=hatching, stroke_width=1)

size = 1
pt1 = vm.Point2D((0, 0))
pt2 = vm.Point2D((0, size))
pt3 = vm.Point2D((size, size))
pt4 = vm.Point2D((size, 0))
c1 = vm.Contour2D([vm.LineSegment2D(pt1, pt2),
                   vm.LineSegment2D(pt2, pt3),
                   vm.LineSegment2D(pt3, pt4),
                   vm.LineSegment2D(pt4, pt1)])

d = c1.plot_data(plot_data_states=[plot_data_state])
plot_data.plot_d3([d])