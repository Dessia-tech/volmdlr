import random
import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
from volmdlr import plot_data

#Points = []

#for i in range(50):
#    Points.append(vm.Point2D.random(0,400,0,400))

#ax = plt.subplots()[1]

#for p in Points:
#    p.MPLPlot(ax=ax)

#### Point test ####
plot_datas = []
for i in range(50):
    point = vm.Point2D.random(0,0.1,0,0.1)
    plot_datas += [point.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'))])]
sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)