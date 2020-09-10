import random
import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
from volmdlr import plot_data

#### Point test ####
plot_datas = []
window_size = plot_data.WindowSizeSet(height=0.3,width=0.2)
shape_set = plot_data.PointShapeSet(shape='circle')
for i in range(50):
    point = vm.Point2D.random(0,window_size.height,0,window_size.width)
    plot_datas += [point.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set)])]
sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)