import math
import volmdlr as vm
import volmdlr.edges as vme
import matplotlib.pyplot as plt

# degree = 3
# points = [vm.Point2D(0, 0),
#           vm.Point2D(1, 1),
#           vm.Point2D(2, 0),
#           vm.Point2D(3, 0)]
# bezier_curve2d = vme.BsplineCurve2D(degree=degree,
#                                    control_points=points,
#                                    name='bezier curve 1')
# _, ax = plt.subplots()
# bezier_curve2d.plot(ax=ax)
# [p.plot(ax=ax) for p in points]


# points = ((0, 0), (3, 4), (-1, 4), (-4, 0), (-4, -3))
degree = 3  # cubic curve

points = [vm.Point3D(0, 0, 0),
          vm.Point3D(0.2, 0.6, 0.3),
          vm.Point3D(1, 1, 2),
          vm.Point3D(2, 1, 1),
          vm.Point3D(3, 0, 4)]

# Do global curve interpolation
bezier_curve3d = vme.BSplineCurve3D.interpolate(points, degree)

point_split = bezier_curve3d.point_at_abscissa(0.2)
bezier_curve3d_s1, bezier_curve3d_s2 =  bezier_curve3d.split(point_split)
abscissa_point_split = bezier_curve3d.abscissa(point_split)

assert math.isclose(abscissa_point_split, 0.2, abs_tol=1e-6)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
ax = bezier_curve3d.plot()
[p.plot(ax=ax) for p in points]
point_split.plot(ax=ax, color='r')

ax2 = bezier_curve3d_s1.plot()
bezier_curve3d_s2.plot(ax=ax2, color='b')
point_split.plot(ax=ax2, color='r')
