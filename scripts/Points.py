import random
import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
from volmdlr import plot_data

Points = []

for i in range(50):
    Points.append(vm.Point2D.random(0,400,0,400))

#ax = plt.subplots()[1]

#for p in Points:
#    p.MPLPlot(ax=ax)

radius = 0.5

for p in Points:

    plot_datas = [plot_data.PlotDataCircle2D(cx=random.randrange(0,400), cy=random.randrange(0,400), r=radius, plot_data_states=[plot_data.PlotDataState(color_surface='black')]) for k in range(50)]
    sol = [plt.to_dict() for plt in plot_datas]
plot_data.plot_d3(sol)