import volmdlr as vm
from volmdlr import plot_data

#### Point test ####
plot_datas = []
window_size = plot_data.WindowSizeSet(width=2,height=1)
shape_set = plot_data.PointShapeSet(shape='crux')
point_size = plot_data.PointSizeSet(size=3)
point_color = plot_data.PointColorSet(color_fill='violet', color_stroke='black')
for i in range(50):
    point = vm.Point2D.random(0,window_size.width,0,window_size.height).to_canvas_style()
    plot_datas += [point.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size,stroke_width=2, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
#point0 = vm.Point2D([0,0])
#plot_datas += [point0.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
#point1 = vm.Point2D([1,1])
#plot_datas += [point1.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
#point2 = vm.Point2D([2,2])
#plot_datas += [point2.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]
#point3 = vm.Point2D([-1,1])
#plot_datas += [point1.plot_data([plot_data.PlotDataState(color_surface=plot_data.ColorSurfaceSet(color='black'), window_size=window_size, shape_set=shape_set, point_size=point_size,point_color=point_color)])]

scatter_plot = vm.ScatterPlot(nb_points_x=10, nb_points_y=10, font_size=15, graduation_color='black')
plot_datas += [scatter_plot.plot_data([plot_data.PlotDataState()])]
sol = [c.to_dict() for c in plot_datas]
plot_data.plot_d3(sol)