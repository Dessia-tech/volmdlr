import volmdlr as vm
import volmdlr.edges as vme
import matplotlib.pyplot as plt


degree = 5
control_points = [vm.Point3D(0, 0, 0),
                  vm.Point3D(3, 2, 1),
                  vm.Point3D(5, -1, 4),
                  vm.Point3D(5, -4, 0),
                  vm.Point3D(-1, -2, -3),
                  vm.Point3D(-3, 4, 1)]
knots = [0.0, 1.0]
knot_multiplicities = [6, 6]
weights = None  # [1, 2, 1, 2, 1, 2]
bspline_curve3d = vme.BSplineCurve3D(degree=degree,
                                     control_points=control_points,
                                     knot_multiplicities=knot_multiplicities,
                                     knots=knots,
                                     weights=weights,
                                     periodic=False,
                                     name='B Spline Curve 3D 1')


ax = bspline_curve3d.plot()

l = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
for absc in l:
    pt = bspline_curve3d.point_at_abscissa(absc*bspline_curve3d.length())
    if absc != 0:
        pt.plot(ax=ax)
    else:
        pt.plot(ax=ax, color='g')
trimmed_bspline = bspline_curve3d.trim_between_evaluations(0.5, 0.9)
trimmed_bspline.plot(ax=ax, color='r')

cut_bspline = bspline_curve3d.cut_after(0.45)
cut_bspline.plot(ax=ax, color='b')

cut_bspline2 = bspline_curve3d.cut_before(0.95)
cut_bspline2.plot(ax=ax, color='g')
