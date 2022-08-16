"""
Test code for intersections between cylinders
Generate random cylinders and create the casing for them
"""

from time import time
from math import acos
from scipy.optimize import minimize, NonlinearConstraint
from volmdlr.primitives3d import Cylinder
import volmdlr as vm


def rotation_angle_z(vector):
    # Used to compute the angle between the global X-axis and a vector
    # This angle will be used for the rotation of the reference frame

    z = vm.Z3D

    dot = z.dot(vector)
    theta = acos(dot / (vector.norm() * z.norm()))

    return theta


def rotation_axis_z(vector):
    # Used to compute a normal vector to the plane containing global X axis and vector
    # This vector will be used as the axis of rotation of the reference frame

    z = vm.Z3D
    vector2 = vector - z

    normal = z.cross(vector2)

    normal.normalize()

    return normal


def new_frame(origin, axis):
    o = vm.O3D
    x = vm.X3D
    y = vm.Y3D
    z = vm.Z3D

    axis.normalize()

    if axis == z:
        # The local frame is oriented like the global frame
        frame = vm.Frame3D(origin, x, y, z)

    elif axis == -z:
        frame = vm.Frame3D(origin, -x, -y, -z)

    else:
        # The local frame is oriented differently from the global frame
        rot_axis = rotation_axis_z(axis)
        rot_angle = rotation_angle_z(axis)

        u = x.rotation(o, rot_axis, rot_angle)
        v = y.rotation(o, rot_axis, rot_angle)
        w = z.rotation(o, rot_axis, rot_angle)

        frame = vm.Frame3D(origin, u, v, w)

    return frame


def min_distance_between_cylinders_axis(cyl0: Cylinder, cyl1: Cylinder):
    """
    Compute the minimal distance between the axis of two volmdlr cylinders

    :param cyl0: volmdlr Cylinder
    :param cyl1: volmdlr Cylinder
    :return: minimal distance between those 3D cylinders axis
    """

    def dist_cyl_axis(
        x,
        axis_0: vm.Vector3D,
        axis_1: vm.Vector3D,
        center_0: vm.Point3D,
        center_1: vm.Point3D,
    ):
        """
        :param x: distance along cylinder 0 axis (from center point), distance along cylinder 0 axis (from center point
        :param axis_0: axis of cylinder 0
        :param axis_1: axis of cylinder 1
        :param center_0: center point of cylinder 0
        :param center_1: center point of cylinder 1
        :return: distance between the axis at given position along their axis
        """
        axis_0.normalize()
        axis_1.normalize()

        p0 = center_0 + axis_0 * x[0]
        p1 = center_1 + axis_1 * x[1]

        return p0.point_distance(p1)

    args = (cyl0.axis, cyl1.axis, cyl0.position, cyl1.position)
    bounds = [(-cyl0.length / 2, cyl0.length / 2), (-cyl1.length / 2, cyl1.length / 2)]

    return minimize(fun=dist_cyl_axis, x0=(0, 0), bounds=bounds, args=args).fun


def is_intersecting_1(cyl0: Cylinder, cyl1: Cylinder):
    """
    :param cyl0: volmdlr Cylinder
    :param cyl1: volmdlr Cylinder
    :return: boolean, True if cylinders are intersecting, False otherwise
    """
    dist = min_distance_between_cylinders_axis(cyl0, cyl1)

    return dist < cyl0.radius + cyl1.radius


def min_distance_between_cylinders(cyl0: Cylinder, cyl1: Cylinder):
    """
    Compute the minimal distance between two volmdlr cylinders

    :param cyl0: volmdlr Cylinder
    :param cyl1: volmdlr Cylinder
    :return: minimal distance between those 3D cylinders
    """

    # Objective function
    def dist_points(x, frame_0: vm.Frame3D, frame_1: vm.Frame3D):
        """
        :param x: coords of a point in cylinder 0 local frame, coords of a point in cylinder 1 local frame
        :param frame_0: local frame of cylinder 0
        :param frame_1: local frame of cylinder 1
        :return: distance between the two points
        """
        p0 = frame_0.old_coordinates(vm.Point3D(x[0], x[1], x[2]))
        p1 = frame_1.old_coordinates(vm.Point3D(x[3], x[4], x[5]))

        return p0.point_distance(p1)

    # Arguments: local frames of cylinders
    frame0 = new_frame(cyl0.position, cyl0.axis)
    frame1 = new_frame(cyl1.position, cyl1.axis)

    args = (frame0, frame1)

    # Initial vector
    p0 = frame0.old_coordinates(vm.O3D)
    p1 = frame1.old_coordinates(vm.O3D)
    x0 = (p0.x, p0.y, p0.z, p1.x, p1.y, p1.z)

    # Constraints
    def c0_0(x):
        # radius of cylinder 0
        return x[0] ** 2 + x[1] ** 2

    def c0_1(x):
        # height of cylinder 0
        return x[2]

    def c1_0(x):
        # radius of cylinder 1
        return x[3] ** 2 + x[4] ** 2

    def c1_1(x):
        # height of cylinder 1
        return x[5]

    constraints = [
        NonlinearConstraint(fun=c0_0, lb=0, ub=cyl0.radius ** 2),
        NonlinearConstraint(fun=c0_1, lb=-cyl0.length / 2, ub=cyl0.length / 2),
        NonlinearConstraint(fun=c1_0, lb=0, ub=cyl1.radius ** 2),
        NonlinearConstraint(fun=c1_1, lb=-cyl1.length / 2, ub=cyl1.length / 2),
    ]

    return minimize(fun=dist_points, x0=x0, constraints=constraints, args=args).fun


def is_intersecting_2(cyl0: Cylinder, cyl1: Cylinder):
    """
    :param cyl0: volmdlr Cylinder
    :param cyl1: volmdlr Cylinder
    :return: boolean, True if cylinders are intersecting, False otherwise
    """
    dist = min_distance_between_cylinders(cyl0, cyl1)

    return dist < 1e-4


cylinders = [
    Cylinder(
        position=vm.Point3D(0, 0.1, 0),
        axis=vm.Vector3D(1, 0, 0),
        radius=0.01,
        length=0.1,
    ),
    Cylinder(
        position=vm.Point3D(0, 0.05, 0),
        axis=vm.Vector3D(1, 1, 0),
        radius=0.005,
        length=0.1,
    ),
    Cylinder(
        position=vm.Point3D(0, 0.159, 0),
        axis=vm.Vector3D(1, 1, 0),
        radius=0.02,
        length=0.1,
    ),
]

cylinders[0].color = (1, 0, 0)
cylinders[1].color = (0, 1, 0)
cylinders[2].color = (0, 0, 1)

volume_model = vm.core.VolumeModel(cylinders)
volume_model.babylonjs()

print("Method 1: with axis and radius")
start = time()

print("Intersecting red & green:", is_intersecting_1(cylinders[0], cylinders[1]))
print("Intersecting red & blue:", is_intersecting_1(cylinders[0], cylinders[2]))
print("Intersecting green & blue:", is_intersecting_1(cylinders[1], cylinders[2]))

print(f"\nIntersection computing duration: {time() - start}s")

print("\n\nMethod 2: with constraints")
start = time()

print("Intersecting red & green:", is_intersecting_2(cylinders[0], cylinders[1]))
print("Intersecting red & blue:", is_intersecting_2(cylinders[0], cylinders[2]))
print("Intersecting green & blue:", is_intersecting_2(cylinders[1], cylinders[2]))

print(f"\nIntersection computing duration: {time() - start}s")
