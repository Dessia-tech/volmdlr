import volmdlr as vm
from volmdlr import plot_data
import numpy as np

# PlotDataState
plot_datas = []
surface_color = 'black'
stroke_width = 0.5  # Points' stroke width

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

# Graph2D
point_list = []
k = 0
while k < 20 * np.pi:
    point = vm.Point2D([k, np.sin(k)])
    point_list.append(point)
    k = k + np.pi/20
dashline=[]
graph_colorstroke = 'black'
graph_linewidth = 0.5
point_colorfill = 'violet'
point_colorstroke = 'grey'
point_strokewidth = 0.5
graph_point_size = 2;
point_shape = 'circle'


graph = vm.Graph2D(point_list=point_list, dashline=dashline,
                   graph_colorstroke=graph_colorstroke, graph_linewidth=graph_linewidth,
                   point_colorfill=point_colorfill, point_colorstroke=point_colorstroke,
                   point_strokewidth=point_strokewidth, graph_point_size=graph_point_size,
                   point_shape=point_shape)
plot_datas += [graph.plot_data([plot_data.PlotDataState()])]

point_list1 = []
k = 0
while k < 20 * np.pi:
    point = vm.Point2D([k, np.sin(k + np.pi/3)])
    point_list1.append(point)
    k = k + np.pi/20
graph1 = vm.Graph2D(point_list=point_list1, dashline=[10,10], graph_colorstroke='red', graph_linewidth=0.5, point_colorfill='green', point_colorstroke='black', point_strokewidth=0.5, graph_point_size=2, point_shape='square')
plot_datas += [graph1.plot_data([plot_data.PlotDataState()])]


point_list2 = []
k = 0
while k < 20 * np.pi:
    point = vm.Point2D([k, np.sin(k - np.pi/3)])
    point_list2.append(point)
    k = k + np.pi/20
graph2 = vm.Graph2D(point_list=point_list2, dashline=[5,3,1,3], graph_colorstroke='blue', graph_linewidth=0.5, point_colorfill='yellow', point_colorstroke='black', point_strokewidth=0.5, graph_point_size=2, point_shape='crux')
plot_datas += [graph2.plot_data([plot_data.PlotDataState()])]

scatter_plot = vm.ScatterPlot(nb_points_x=nb_points_x, nb_points_y=nb_points_y,
                              font_size=font_size,
                              graduation_color=graduation_color,
                              axis_color=axis_color, arrow_on=arrow_on,
                              axis_width=axis_width, grid_on=grid_on)
plot_datas += [scatter_plot.plot_data([plot_data.PlotDataState()])]

tooltip = vm.Tooltip(colorfill=colorfill, font=font, tp_width=tp_width,
                     tp_radius=tp_radius, to_plot_list=to_plot_list)
plot_datas += [tooltip.plot_data([plot_data.PlotDataState()])]

sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)