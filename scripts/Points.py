import random
import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
from volmdlr import plot_data

#### Point test ####
plot_datas = []
window_size = plot_data.WindowSizeSet(width=2,height=1)
shape_set = plot_data.PointShapeSet(shape='crux')
point_size = plot_data.PointSizeSet(size=3)
point_color = plot_data.PointColorSet(color_fill='black', color_stroke='red')
for i in range(50):
    point = vm.Point2D.random(0,window_size.width,0,window_size.height)
    plot_datas += [point.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)