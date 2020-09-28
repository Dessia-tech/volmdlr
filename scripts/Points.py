import volmdlr as vm
from volmdlr import plot_data

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
color_fill = 'violet'
color_stroke = 'grey'

# PlotDataState
surface_color = 'black'
stroke_width = 2  # Points' stroke width

# Scatter plot
nb_points_x = 10
nb_points_y = 10
font_size = 15
graduation_color = 'grey'
axis_color = 'grey'
arrow_on = False

# Tooltip
colorfill = 'lightblue'
font = '15px'  # Font family : Arial, Helvetica, serif, sans-serif, Verdana, Times New Roman, Courier New
tp_width = 90
tp_height = 35
tp_radius = 10
to_plot_list = ['cx', 'cy']

plot_datas = []
window_size = plot_data.WindowSizeSet(width=width, height=height)
shape_set = plot_data.PointShapeSet(shape=shape)
point_size = plot_data.PointSizeSet(size=size)
point_color = plot_data.PointColorSet(color_fill=color_fill,
                                      color_stroke=color_stroke)
for i in range(50):
    point = vm.Point2D.random(0, window_size.width, 0,
                              window_size.height).to_canvas_style()
    plot_datas += [point.plot_data([plot_data.PlotDataState(
        color_surface=plot_data.ColorSurfaceSet(color=surface_color),
        window_size=window_size, stroke_width=stroke_width,
        shape_set=shape_set, point_size=point_size, point_color=point_color)])]
# point0 = vm.Point2D([0,0])
# plot_datas += [point0.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
# point1 = vm.Point2D([1,1])
# plot_datas += [point1.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
# point2 = vm.Point2D([2,2])
# plot_datas += [point2.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
# point3 = vm.Point2D([-1,1])
# plot_datas += [point1.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]

scatter_plot = vm.ScatterPlot(nb_points_x=nb_points_x, nb_points_y=nb_points_y,
                              font_size=font_size,
                              graduation_color=graduation_color, axis_color=axis_color, arrow_on=arrow_on)
plot_datas += [scatter_plot.plot_data([plot_data.PlotDataState()])]

tooltip = vm.Tooltip(colorfill=colorfill, font=font, tp_width=tp_width, tp_radius=tp_radius, to_plot_list=to_plot_list)
plot_datas += [tooltip.plot_data([plot_data.PlotDataState()])]

sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)
