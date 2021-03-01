import volmdlr as vm
import volmdlr.faces as vmf
import matplotlib.pyplot as plt

degree_u = 3
degree_v = 2

control_points = [vm.Point3D(0, 0, 0), vm.Point3D(0, 4, 0), vm.Point3D(0, 8, -3),
                  vm.Point3D(2, 0, 6), vm.Point3D(2, 4, 0), vm.Point3D(2, 8, 0),
                  vm.Point3D(4, 0, 0), vm.Point3D(4, 4, 0), vm.Point3D(4, 8, 3),
                  vm.Point3D(6, 0, 0), vm.Point3D(6, 4, -3), vm.Point3D(6, 8, 0)]

bezier_surface = vmf.BezierSurface3D(degree_u=degree_u, degree_v=degree_v,
                                     control_points=control_points, nb_u=4,
                                     nb_v=3, name='bezier curve 1')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
bezier_surface.plot(ax=ax)
[p.plot(ax=ax) for p in control_points]