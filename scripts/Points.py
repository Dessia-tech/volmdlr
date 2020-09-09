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

#### Circle test ####
#center = vm.Point2D((0,0))
#radius = 0.5
#circle = vm.Circle2D(center,radius)
#contour = vm.Contour2D([circle])
#contour1 = [contour.plot_data([plot_data.PlotDataState()])]
#print(contour1)

#### Point test ####
point = vm.Point2D.random(0,400,0,400)

contour1 = [point.plot_data([plot_data.PlotDataState()])]
print(contour1)
sol = [c.to_dict() for c in contour1]
plot_data.plot_d3(sol)