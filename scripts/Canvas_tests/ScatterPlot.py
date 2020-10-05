import volmdlr as vm
from volmdlr import plot_data
import numpy as np

# Point test ####

# PARAMETERS #
# Window size
width = 2
height = 1

# Shape set (circle, square, crux)
shape = 'circle'

# Point size (1 to 4)
size = 2

# Points' color
colorfill = 'violet'
colorstroke = 'grey'

# PlotDataState
surface_color = 'black'
strokewidth = 0.5  # Points' stroke width

# Scatter plot
nb_points_x = 10
nb_points_y = 10
font_size = 12
graduation_color = 'grey'
axis_color = 'grey'
axis_width = 0.5
arrow_on = False
grid_on = True

# Tooltip
colorfill = 'lightblue'
font = '12px sans-serif'  # Font family : Arial, Helvetica, serif, sans-serif, Verdana, Times New Roman, Courier New
tp_width = 90
tp_radius = 10
to_plot_list = ['cx', 'cy']

plot_datas = []
point_list = []
for i in range(50):
    point = vm.Point2D.random(0, 2, 0, 1)
    point_list += [point]

# point_list += [vm.Point2D([0,0])]
# point_list += [vm.Point2D([1,0])]
ScatterPlot = vm.ScatterPlot(point_list, colorfill=colorfill, colorstroke=colorstroke, size=size, shape=shape, strokewidth=strokewidth)
plot_datas += [ScatterPlot.plot_data([plot_data.PlotDataState()])]

axis = vm.Axis(nb_points_x=nb_points_x, nb_points_y=nb_points_y,
                              font_size=font_size,
                              graduation_color=graduation_color,
                              axis_color=axis_color, arrow_on=arrow_on,
                              axis_width=axis_width, grid_on=grid_on)

plot_datas += [axis.plot_data([plot_data.PlotDataState()])]

tooltip = vm.Tooltip(colorfill=colorfill, font=font, tp_width=tp_width,
                     tp_radius=tp_radius, to_plot_list=to_plot_list)
plot_datas += [tooltip.plot_data([plot_data.PlotDataState()])]

sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)