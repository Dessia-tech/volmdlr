import random
import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
from volmdlr import plot_data

Points = []

#for i in range(50):
#    Points.append(vm.Point2D.random(0,400,0,400))

#ax = plt.subplots()[1]

#for p in Points:
#    p.MPLPlot(ax=ax)

center = vm.Point2D((0,0))
radius = 0.5

plot_datas = []
for i in range(50):
    point = vm.Point2D.random(0,400,0,400)
    contour = vm.Contour2D([point])
    plot_datas += [contour.plot_data([plot_data.PlotDataState()])] #PlotDataState doit être dans une liste []
sol = [plots.to_dict() for plots in plot_datas]
plot_data.plot_d3(sol)