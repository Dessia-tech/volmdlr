import volmdlr as vm
import volmdlr.edges as vme
import matplotlib.pyplot as plt

degree = 2
points = [vm.Point2D(0, 0),
          vm.Point2D(1, 1),
          vm.Point2D(2, 1),
          vm.Point2D(3, 0)]
bezier_curve2d = vme.BezierCurve2D(degree=degree,
                                   control_points=points,
                                   name='bezier curve 1')
_, ax = plt.subplots()
bezier_curve2d.plot(ax=ax)
[p.plot(ax=ax) for p in points]

degree = 3
points = [vm.Point3D(0, 0, 0),
          vm.Point3D(1, 1, 2),
          vm.Point3D(2, 1, 1),
          vm.Point3D(3, 0, 4)]
bezier_curve3d = vme.BezierCurve3D(degree=degree,
                                   control_points=points,
                                   name='bezier curve 1')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
bezier_curve3d.plot(ax=ax)
[p.plot(ax=ax) for p in points]