"""
Test code for intersections between cylinders
Generate random cylinders and create the casing for them
"""

from volmdlr.primitives3d import Cylinder
import volmdlr as vm
from scipy.optimize import minimize
from time import time


def min_distance_between_cylinders_axis(cyl0: Cylinder, cyl1: Cylinder):
    """
    :param cyl0: volmdlr Cylinder
    :param cyl1: volmdlr Cylinder
    :return: minimal distance between those 3D cylinders
    """
    def dist_cyl(x, axis_0: vm.Vector3D, axis_1: vm.Vector3D, center_0: vm.Point3D, center_1: vm.Point3D):
        axis_0.normalize()
        axis_1.normalize()

        p0 = center_0 + axis_0 * x[0]
        p1 = center_1 + axis_1 * x[1]

        return p0.point_distance(p1)

    args = (cyl0.axis, cyl1.axis, cyl0.position, cyl1.position)
    bounds = [(-cyl0.length/2, cyl0.length/2), (-cyl1.length/2, cyl1.length/2)]

    return minimize(fun=dist_cyl, x0=(0, 0), bounds=bounds, args=args).fun


def is_intersecting(cyl0: Cylinder, cyl1: Cylinder):
    """
    :param cyl0: volmdlr Cylinder
    :param cyl1: volmdlr Cylinder
    :return: boolean, True if cylinders are intersecting, False otherwise
    """
    dist = min_distance_between_cylinders_axis(cyl0, cyl1)

    return dist < cyl0.radius + cyl1.radius


cylinders = [Cylinder(position=vm.Point3D(0, 0.1, 0), axis=vm.Vector3D(1, 0, 0), radius=.01, length=.1),
             Cylinder(position=vm.Point3D(0, 0.05, 0), axis=vm.Vector3D(1, 1, 0), radius=.01, length=.1),
             Cylinder(position=vm.Point3D(0, 0.158, 0), axis=vm.Vector3D(1, 1, 0), radius=.01, length=.1)]

cylinders[0].color = (1, 0, 0)
cylinders[1].color = (0, 1, 0)
cylinders[2].color = (0, 0, 1)

# volume_model = vm.core.VolumeModel(cylinders)
# volume_model.babylonjs()

start = time()

print("Intersecting red & green:", is_intersecting(cylinders[0], cylinders[1]))
print("Intersecting red & blue:", is_intersecting(cylinders[0], cylinders[2]))
print("Intersecting green & blue:", is_intersecting(cylinders[1], cylinders[2]))

print(f"\nIntersection computing duration: {time()-start}s")
