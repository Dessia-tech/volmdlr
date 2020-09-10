import random
import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
from volmdlr import plot_data

#### Point test ####
plot_datas = []
window_size = plot_data.WindowSizeSet(width=0.1,height=0.1)
shape_set = plot_data.PointShapeSet(shape='square')
for i in range(50):
    point = vm.Point2D.random(0,window_size.width,0,window_size.height)
    plot_datas += [point.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set)])]

### LineSegment2D test ###

#line = vm.LineSegment2D(vm.Point2D([0,0]), vm.Point2D([0.1,0.1]))
#contour = vm.Contour2D([line])
#plot_datas = [contour.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), color_line='red')])]
sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)