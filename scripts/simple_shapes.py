#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:32:37 2018

@author: Steven Masfaraud masfaraud@dessia.tech
"""
import volmdlr as vm
#import volmdlr.primitives2D as primitives2D
from volmdlr import plot_data
import numpy as npy

#for i in range(20):
triangle_points=[vm.Point2D(npy.random.random(2)) for i in range(3)]
triangle=vm.Polygon2D(triangle_points)


cog_triangle=triangle.CenterOfMass()
c1 = vm.CompositePrimitive2D([triangle, cog_triangle])
#c1.MPLPlot()

print(triangle.Area())

p0=vm.Point2D((-1,0))
p1=vm.Point2D((-npy.cos(npy.pi/4),npy.sin(npy.pi/4)))
p2=vm.Point2D((0,1))

a = vm.Arc2D(p2,p1,p0)
l = vm.LineSegment2D(p2,a.center)
#list_node = a.Discretise()

c = vm.Contour2D([a, l])
#print(c.plot_data())
c2 = vm.CompositePrimitive2D([c])
#c2.MPLPlot()
#print(c.Area())

hatching = plot_data.HatchingSet(0.5, 3)
color_surface = plot_data.ColorSurfaceSet(color='white')
plot_data_state = plot_data.PlotDataState(name='be_sup', hatching=hatching, stroke_width=1)
plot_datas = [c.plot_data(plot_data_states=[plot_data_state])]
sol = [plt.to_dict() for plt in plot_datas]
plot_data.plot_d3(sol)
