#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import matplotlib.pyplot as plt
import dessia_common as dc
from geomdl import BSpline
import volmdlr
import volmdlr.faces3d



class Surface3D():
    """
    Abstract class
    """
    def face_from_contours3d(self, contours3d):
        contours2d = []
        max_area = 0.

        for ic, contour3d in contours3d:
            contour2d = self.contour3d_to_2d(contour3d)
            contour_area = contour2d.Area()
            if contour_area > max_area:
                max_area = contour_area
                outer_contour_index = ic
            contours2d.append(contour2d)

        outer_contour = contour2d[outer_contour_index]
        del contour2d[outer_contour_index]
        surface2d = (outer_contour, contours2d)

        self.face_class(self, surface2d)

class Plane3D(Surface3D):
    def __init__(self, frame: volmdlr.Frame3D, name: str = ''):
        """
        :param frame: u and v of frame describe the plane, w is the normal
        """
        self.frame = frame
        self.name = name

    def __hash__(self):
        return hash(self.frame)

    def __eq__(self, other_plane):
        return (self.frame.origin == other_plane.frame.origin and \
                self.frame.w.is_colinear_to(other_plane.frame.w))

    def to_dict(self):
        # improve the object structure ?
        dict_ = dc.DessiaObject.base_dict(self)
        dict_['frame'] = self.frame.to_dict()
        dict_['name'] = self.name
        dict_['object_class'] = 'volmdlr.core.Plane3D'
        return dict_

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]

        return cls(frame3d, arguments[0][1:-1])

    @classmethod
    def from_3_points(cls, point1, point2, point3):
        """
        Point 1 is used as origin of the plane
        """
        vector1 = point2 - point1
        vector2 = point3 - point1
        vector1.Normalize()
        vector2.Normalize()
        normal = vector1.Cross(vector2)
        normal.Normalize()
        vector = normal.Cross(vector1)
        frame = volmdlr.Frame3D(point1, vector1, normal.Cross(vector1), normal)
        return cls(frame)

    @classmethod
    def from_normal(cls, point, normal):
        v1 = normal.deterministic_unit_normal_vector()
        v2 = v1.Cross(normal)
        return cls(volmdlr.Frame3D(point, v1, v2, normal))

    @classmethod
    def from_points(cls, points):
        if len(points) < 3:
            raise ValueError
        elif len(points) == 3:
            return cls.from_3_points(volmdlr.Point3D(points[0].vector),
                                     volmdlr.Vector3D(points[1].vector),
                                     volmdlr.Vector3D(points[2].vector))
        else:
            points = [p.copy() for p in points]
            indexes_to_del = []
            for i, point in enumerate(points[1:]):
                if point == points[0]:
                    indexes_to_del.append(i)
            for index in indexes_to_del[::-1]:
                del points[index + 1]

            origin = volmdlr.Point3D(points[0].vector)
            vector1 = volmdlr.Vector3D(points[1] - origin)
            vector1.Normalize()
            vector2_min = volmdlr.Vector3D(points[2] - origin)
            vector2_min.Normalize()
            dot_min = abs(vector1.Dot(vector2_min))
            for point in points[3:]:
                vector2 = volmdlr.Vector3D(point - origin)
                vector2.Normalize()
                dot = abs(vector1.Dot(vector2))
                if dot < dot_min:
                    vector2_min = vector2
                    dot_min = dot
            return cls.from_3_points(origin, vector1 + origin,
                                     vector2_min + origin)

    def point_on_plane(self, point):
        if math.isclose(self.frame.w.Dot(point - self.frame.origin), 0,
                        abs_tol=1e-6):
            return True
        return False

    def line_intersections(self, line):
        u = line.points[1] - line.points[0]
        w = line.points[0] - self.frame.origin
        if math.isclose(self.frame.w.Dot(u), 0, abs_tol=1e-08):
            return []
        intersection_abscissea = - self.frame.w.Dot(w) / self.frame.w.Dot(u)
        return [line.points[0] + intersection_abscissea * u]

    def linesegment_intersection(self, linesegment, abscissea=False):
        u = linesegment.points[1] - linesegment.points[0]
        w = linesegment.points[0] - self.frame.origin
        normalDotu = self.frame.w.Dot(u)
        if math.isclose(normalDotu, 0, abs_tol=1e-08):
            if abscissea:
                return None, None
            return None
        intersection_abscissea = - self.frame.w.Dot(w) / normalDotu
        if intersection_abscissea < 0 or intersection_abscissea > 1:
            if abscissea:
                return None, None
            return None
        if abscissea:
            return linesegment.points[
                       0] + intersection_abscissea * u, intersection_abscissea
        return linesegment.points[0] + intersection_abscissea * u

    def equation_coefficients(self):
        """
        returns the a,b,c,d coefficient from equation ax+by+cz+d = 0
        """
        a, b, c = self.frame.w
        d = -self.frame.origin.Dot(self.frame.w)
        return (a, b, c, d)

    def plane_intersection(self, other_plane):
        line_direction = self.frame.w.Cross(other_plane.frame.w)

        if line_direction.Norm() < 1e-6:
            return None

        a1, b1, c1, d1 = self.equation_coefficients()
        a2, b2, c2, d2 = other_plane.equation_coefficients()

        if a1 * b2 - a2 * b1 != 0.:
            x0 = (b1 * d2 - b2 * d1) / (a1 * b2 - a2 * b1)
            y0 = (a2 * d1 - a1 * d2) / (a1 * b2 - a2 * b1)
            point1 = volmdlr.Point3D((x0, y0, 0))
        else:
            y0 = (b2 * d2 - c2 * d1) / (b1 * c2 - c1 * b2)
            z0 = (c1 * d1 - b1 * d2) / (b1 * c2 - c1 * b2)
            point1 = volmdlr.Point3D((0, y0, z0))

        point2 = point1 + line_direction
        return volmdlr.Line3D(point1, point2)

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_frame = self.frame.Rotation(center, axis, angle, copy=True)
            return Plane3D(new_frame)
        else:
            self.frame.Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_frame = self.frame.Translation(offset, True)
            return Plane3D(new_frame)
        else:
            self.origin.Translation(offset, False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            new_origin = frame.OldCoordinates(self.origin)
            new_vector1 = frame.Basis().OldCoordinates(self.vectors[0])
            new_vector2 = frame.Basis().OldCoordinates(self.vectors[1])
            if copy:
                return Plane3D(new_origin, new_vector1, new_vector2, self.name)
            else:
                self.origin = new_origin
                self.vectors = [new_vector1, new_vector2]
                self.normal = frame.Basis().OldCoordinates(self.normal)
                self.normal.Normalize()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.origin)
            new_vector1 = frame.Basis().NewCoordinates(self.vectors[0])
            new_vector2 = frame.Basis().NewCoordinates(self.vectors[1])
            if copy:
                return Plane3D(new_origin, new_vector1, new_vector2, self.name)
            else:
                self.origin = new_origin
                self.vectors = [new_vector1, new_vector2]
                self.normal = frame.Basis().NewCoordinates(self.normal)
                self.normal.Normalize()

    def copy(self):
        new_origin = self.origin.copy()
        new_vector1 = self.vectors[0].copy()
        new_vector2 = self.vectors[1].copy()
        return Plane3D(new_origin, new_vector1, new_vector2, self.name)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        self.origin.MPLPlot(ax)
        self.vectors[0].MPLPlot(ax, starting_point=self.origin, color='r')
        self.vectors[1].MPLPlot(ax, starting_point=self.origin, color='g')
        return ax

    def babylon_script(self):
        s = 'var myPlane = BABYLON.MeshBuilder.CreatePlane("myPlane", {width: 0.5, height: 0.5, sideOrientation: BABYLON.Mesh.DOUBLESIDE}, scene);\n'
        s += 'myPlane.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n'.format(
            self.origin[0], self.origin[1], self.origin[2])

        s += 'var axis1 = new BABYLON.Vector3({}, {}, {});\n'.format(
            self.vectors[0][0], self.vectors[0][1], self.vectors[0][2])
        s += 'var axis2 = new BABYLON.Vector3({}, {}, {});\n'.format(
            self.vectors[1][0], self.vectors[1][1], self.vectors[1][2])
        s += 'var axis3 = new BABYLON.Vector3({}, {}, {});\n'.format(
            self.normal[0], self.normal[1], self.normal[2])
        s += 'var orientation = BABYLON.Vector3.RotationFromAxis(axis1, axis2, axis3);\n'
        s += 'myPlane.rotation = orientation;\n'

        s += 'var planemat = new BABYLON.StandardMaterial("planemat", scene);\n'
        s += 'planemat.alpha = 0.4;\n'
        s += 'myPlane.material = planemat;\n'

        return s

    def point2d_to_3d(self, point2d):
        return point2d.To3D(self.frame.origin, self.frame.u, self.frame.v)

    def point3d_to_2d(self, point3d):
        return point3d.To2D(self.frame.origin, self.frame.u, self.frame.v)

    def rectangular_cut(self, x1:float, x2:float,
                        y1:float, y2:float, name:str=''):

        p1 = volmdlr.volmdlr.Point2D((x1, y1))
        p2 = volmdlr.volmdlr.Point2D((x2, y1))
        p3 = volmdlr.volmdlr.Point2D((x2, y2))
        p4 = volmdlr.volmdlr.Point2D((x1, y2))
        outer_contour = volmdlr.Polygon2D([p1, p2, p3, p4])
        surface = volmdlr.Surface2D(outer_contour, [])
        return volmdlr.faces3d.Planevolmdlr.faces3d.Face3D(self, surface, name)

XYZ = volmdlr.Basis3D(volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
YZX = volmdlr.Basis3D(volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D)
ZXY = volmdlr.Basis3D(volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D)

OXYZ = volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
OYZX = volmdlr.Frame3D(volmdlr.O3D, volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D)
OZXY = volmdlr.Frame3D(volmdlr.O3D, volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D)

PLANE3D_OXY = Plane3D(OXYZ)
PLANE3D_OYZ = Plane3D(OYZX)
PLANE3D_OZX = Plane3D(OZXY)

class CylindricalSurface3D(Surface3D):
    face_class = volmdlr.faces3d.Cylindricalvolmdlr.faces3d.Face3D
    """
    :param frame: frame.w is axis, frame.u is theta=0 frame.v theta=pi/2
    :type frame: volmdlr.Frame3D
    :param radius: Cylinder's radius
    :type radius: float
    """

    def __init__(self, frame, radius, name=''):
        self.frame = frame
        self.radius = radius
        self.name = name

    def point2d_to_3d(self, point2d):
        p = volmdlr.Point3D(volmdlr.Vector3D([self.radius * math.cos(point2d[0]),
                              self.radius * math.sin(point2d[0]),
                              point2d[1]]))
        return self.frame.OldCoordinates(p)

    def point3d_to_2d(self, point3d):
        x, y, z = point3d

        u1, u2 = x / self.radius, y / self.radius
        theta = volmdlr.sin_cos_angle(u1, u2)
        return volmdlr.Point2D([theta, z])

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, -frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = volmdlr.Frame3D(frame3d.origin, U, V, W)
        radius = float(arguments[2]) / 1000
        return cls(frame_direct, radius, arguments[0][1:-1])

    def frame_mapping(self, frame, side, copy=True):
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return CylindricalSurface3D(new_frame, self.radius,
                                            name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return CylindricalSurface3D(new_frame, self.radius,
                                            name=self.name)
            else:
                self.frame = new_frame

    def rectangular_cut(self, theta1:float, theta2:float,
                        z1:float, z2:float, name:str=''):
        print(theta1, theta2, z1, z2)
        # theta1 = angle_principal_measure(theta1)
        # theta2 = angle_principal_measure(theta2)

        if theta1 == theta2:
            theta2 += volmdlr.volmdlr.two_pi

        p1 = volmdlr.Point2D((theta1, z1))
        p2 = volmdlr.Point2D((theta2, z1))
        p3 = volmdlr.Point2D((theta2, z2))
        p4 = volmdlr.Point2D((theta1, z2))
        outer_contour = volmdlr.Polygon2D([p1, p2, p3, p4])
        surface2d = volmdlr.Surface2D(outer_contour, [])
        return volmdlr.faces3d.Cylindricalvolmdlr.faces3d.Face3D(self, surface2d, name)


class ToroidalSurface3D(Surface3D):
    face_class = volmdlr.faces3d.ToroidalFace3D
    """
    :param frame: Tore's frame: origin is thet center, u is pointing at
                    theta=0
    :type frame: volmdlr.Frame3D
    :param R: Tore's radius
    :type R: float
    :param r: Circle to revolute radius
    :type r: float
    Definitions of R and r according to https://en.wikipedia.org/wiki/Torus
    """

    def __init__(self, frame: volmdlr.Frame3D,
                 R: float, r: float, name: str = ''):
        self.frame = frame
        self.R = R
        self.r = r
        self.name = name
        self.frame = frame

    def _bounding_box(self):
        d = self.R + self.r
        p1 = self.frame.origin + self.frame.u * d + self.frame.v * d + self.frame.w * self.r
        p2 = self.frame.origin + self.frame.u * d + self.frame.v * d - self.frame.w * self.r
        p3 = self.frame.origin + self.frame.u * d - self.frame.v * d + self.frame.w * self.r
        p4 = self.frame.origin + self.frame.u * d - self.frame.v * d - self.frame.w * self.r
        p5 = self.frame.origin - self.frame.u * d + self.frame.v * d + self.frame.w * self.r
        p6 = self.frame.origin - self.frame.u * d + self.frame.v * d - self.frame.w * self.r
        p7 = self.frame.origin - self.frame.u * d - self.frame.v * d + self.frame.w * self.r
        p8 = self.frame.origin - self.frame.u * d - self.frame.v * d - self.frame.w * self.r

        return volmdlr.BoundingBox.from_points([p1, p2, p3, p4, p5, p6, p7, p8])

    def point2d_to_3d(self, point2d):
        theta, phi = point2d
        x = (self.R + self.r * math.cos(phi)) * math.cos(theta)
        y = (self.R + self.r * math.cos(phi)) * math.sin(theta)
        z = self.r * math.sin(phi)
        return self.frame.OldCoordinates(volmdlr.Point3D([x, y, z]))

    def point3d_to_2d(self, point3d):
        # points_2D = []
        x, y, z = point3d
        if z < -self.r:
            z = -self.r
        elif z > self.r:
            z = self.r

        zr = z / self.r
        phi = math.asin(zr)

        u = self.R + math.sqrt((self.r ** 2) - (z ** 2))
        u1, u2 = round(x / u, 5), round(y / u, 5)
        theta = volmdlr.sin_cos_angle(u1, u2)

        return volmdlr.volmdlr.Point2D([theta, phi])

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, -frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = volmdlr.Frame3D(frame3d.origin, U, V, W)
        rcenter = float(arguments[2]) / 1000
        rcircle = float(arguments[3]) / 1000
        return cls(frame_direct, rcenter, rcircle, arguments[0][1:-1])

    def frame_mapping(self, frame, side, copy=True):
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ToroidalSurface3D(new_frame,
                                         self.R, self.r,
                                         name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ToroidalSurface3D(new_frame,
                                         self.R, self.r,
                                         name=self.name)
            else:
                self.frame = new_frame

    def rectangular_cut(self, theta1, theta2, phi1, phi2, name=''):
        # theta1 = angle_principal_measure(theta1)
        # theta2 = angle_principal_measure(theta2)
        # phi1 = angle_principal_measure(phi1)
        # phi2 = angle_principal_measure(phi2)

        if phi1 == phi2:
            phi2 += volmdlr.two_pi
        elif phi2 < phi1:
            phi2 += volmdlr.two_pi
        if theta1 == theta2:
            theta2 += volmdlr.two_pi
        elif theta2 < theta1:
            theta2 += volmdlr.two_pi

        # print(theta1, theta2, phi1, phi2)
        p1 = volmdlr.volmdlr.Point2D((theta1, phi1))
        p2 = volmdlr.volmdlr.Point2D((theta1, phi2))
        p3 = volmdlr.volmdlr.Point2D((theta2, phi2))
        p4 = volmdlr.volmdlr.Point2D((theta2, phi1))
        outer_contour = volmdlr.Polygon2D([p1, p2, p3, p4])
        return self.face_class(self,
                               volmdlr.Surface2D(outer_contour, []),
                               name)

    def contour2d_to3d(self, contour2d):
        edges3d = []
        # contour2d.MPLPlot()
        for edge in contour2d.primitives:
            if isinstance(edge, volmdlr.volmdlr.LineSegment2D):
                if (edge.points[0][0] == edge.points[1][0]) \
                        or (edge.points[0][1] == edge.points[1][1]):
                    # Y progression: it's an arc
                    edges3d.append(volmdlr.Arc3D(self.point2d_to_3d(edge.points[0]),
                                         self.point2d_to_3d(
                                             0.5 * (edge.points[0] \
                                                    + edge.points[1])),
                                         self.point2d_to_3d(edge.points[1])))
                else:
                    edges3d.append(volmdlr.Arc3D(self.point2d_to_3d(edge.points[0]),
                                         self.point2d_to_3d(
                                             0.5 * (edge.points[0] \
                                                    + edge.points[1])),
                                         self.point2d_to_3d(edge.points[1])))
            else:
                raise NotImplementedError(
                    'The primitive {} is not supported in 2D->3D'.format(edge))

        return volmdlr.Contour3D(edges3d)

    def contour3d_to2d(self, contour3d, toroidalsurface3d):
        frame = toroidalsurface3d.frame
        n = frame.w
        # center = frame.origin
        rcenter, rcircle = toroidalsurface3d.rcenter, toroidalsurface3d.rcircle

        primitives, start_end, all_points = [], [], []
        for edge in contour3d.edges:
            new_points = [frame.NewCoordinates(pt) for pt in edge.points]
            if edge.__class__ is volmdlr.Arc3D:
                if edge.normal == n or edge.normal == -n:
                    start2d, end2d = self.points3d_to2d(new_points,
                                                        rcenter,
                                                        rcircle)

                    angle2d = abs(end2d[0] - start2d[0])
                    if math.isclose(edge.angle, volmdlr.two_pi, abs_tol=1e-6):
                        if start2d == end2d:
                            if math.isclose(start2d.vector[0], volmdlr.two_pi,
                                            abs_tol=1e-6):
                                end2d = end2d - volmdlr.volmdlr.Point2D((volmdlr.two_pi, 0))
                            else:
                                end2d = end2d + volmdlr.volmdlr.Point2D((volmdlr.two_pi, 0))
                    elif not (math.isclose(edge.angle, angle2d, abs_tol=1e-2)):
                        # if math.isclose(angle2d, volmdlr.two_pi, abs_tol=1e-2) :
                        if start2d[0] < end2d[0]:
                            end2d = start2d + volmdlr.volmdlr.Point2D((edge.angle, 0))
                        else:
                            end2d = start2d - volmdlr.volmdlr.Point2D((edge.angle, 0))
                    #####################

                    ls_toadd = volmdlr.volmdlr.LineSegment2D(start2d, end2d)
                    same = volmdlr.faces3d.volmdlr.faces3d.Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False:
                        primitives.append([ls_toadd])
                        all_points.extend(ls_toadd.points)
                        start_end.append(ls_toadd.points)

                else:
                    points2d = self.points3d_to2d(new_points,
                                                  rcenter, rcircle)
                    lines = []
                    for k in range(0, len(points2d) - 1):
                        lines.append(
                            volmdlr.LineSegment2D(points2d[k], points2d[k + 1]))
                    points, prim_list = [], []
                    for ls_toadd in lines:
                        same = volmdlr.faces3d.Face3D.LS2D_inprimitives(ls_toadd, primitives)
                        if same is False:
                            prim_list.append(ls_toadd)
                            points.extend(ls_toadd.points)
                    if len(points) > 0:
                        all_points.extend(points)
                        primitives.append(prim_list)
                        start_end.append([points[0], points[-1]])

            elif edge.__class__ is volmdlr.LineSegment3D:
                start2d, end2d = self.points3d_to2d(new_points,
                                                              rcenter, rcircle)
                ls_toadd = volmdlr.LineSegment2D(start2d, end2d)
                same = volmdlr.faces3d.Face3D.LS2D_inprimitives(ls_toadd, primitives)
                if same is False:
                    primitives.append([ls_toadd])
                    all_points.extend(ls_toadd.points)
                    start_end.append(ls_toadd.points)

            else:
                points2d = volmdlr.faces3d.volmdlr.points3d_to2d(new_points, rcenter,
                                                        rcircle)
                lines = []
                for k in range(0, len(points2d) - 1):
                    lines.append(volmdlr.LineSegment2D(points2d[k], points2d[k + 1]))
                points, prim_list = [], []
                for ls_toadd in lines:
                    same = volmdlr.faces3d.Face3D.LS2D_inprimitives(ls_toadd, primitives)
                    if same is False:
                        prim_list.append(ls_toadd)
                        points.extend(ls_toadd.points)
                if len(points) > 0:
                    all_points.extend(points)
                    primitives.append(prim_list)
                    start_end.append([points[0], points[-1]])

        points_se, primitives_se = [], []
        for double in start_end:
            primitives_se.append(volmdlr.LineSegment2D(double[0], double[1]))
            points_se.extend(double)
        poly_se = volmdlr.Polygon2D(points_se)

        xmax, xmin = max(pt[0] for pt in points_se), min(
            pt[0] for pt in points_se)
        ymax, ymin = max(pt[1] for pt in points_se), min(
            pt[1] for pt in points_se)
        pt1, pt2, pt3, pt4 = volmdlr.volmdlr.Point2D((xmin, ymin)), volmdlr.volmdlr.Point2D(
            (xmin, ymax)), volmdlr.volmdlr.Point2D((xmax, ymin)), volmdlr.volmdlr.Point2D((xmax, ymax))
        diag1, diag2 = volmdlr.LineSegment2D(pt1, pt4), volmdlr.LineSegment2D(pt2, pt3)
        diag1_cut, diag2_cut = [], []
        diag1_pointcut, diag2_pointcut = [], []
        for enum, l in enumerate(primitives_se):
            cut1 = diag1.line_intersections(l)
            cut2 = diag2.line_intersections(l)
            if cut1 is not None:
                diag1_cut.append(enum)
                diag1_pointcut.append(cut1)
            if cut2 is not None:
                diag2_cut.append(enum)
                diag2_pointcut.append(cut2)

        points_common = []
        for enum1, pos1 in enumerate(diag1_cut):
            for enum2, pos2 in enumerate(diag2_cut):
                if pos1 == pos2:
                    points_common.append(primitives_se[pos1].points)

        if len(points_common) >= 1:
            solve = False
            for couple in points_common:
                if solve:
                    break
                check1, check2 = poly_se.PointBelongs(
                    couple[0]), poly_se.PointBelongs(couple[1])
                start, end = couple[0].vector[0], couple[1].vector[0]
                if math.isclose(start, end, abs_tol=5e-2):
                    intersect = min(start, end)
                    if math.isclose(intersect, math.pi, abs_tol=5e-2):
                        # all_points = check_singularity(all_points)
                        # intersect = 0

                        ##################### NEW
                        points_sing = volmdlr.check_singularity(all_points)
                        pt0, pt2pi = 0, 0
                        for pt in points_sing:
                            if math.isclose(pt.vector[0], 0, abs_tol=1e-2):
                                pt0 += 1
                            elif math.isclose(pt.vector[0], volmdlr.two_pi,
                                              abs_tol=1e-2):
                                pt2pi += 1
                        points_sing.sort(key=lambda pt: pt[1])
                        points_sing.sort(key=lambda pt: pt[0])
                        if pt2pi != 0 and pt0 == 0:
                            points = [pt.copy() for pt in points_sing[::-1]]
                            points_sing = points
                        points_range = volmdlr.faces3d.Face3D.range_closest(points_sing)
                        all_points = volmdlr.delete_double_point(points_range)
                        break
                        #######################

                    # if math.isclose(intersect, 0, abs_tol = 1e-6) or math.isclose(intersect, volmdlr.two_pi, abs_tol = 1e-6) or (not check1 or not check2):
                    elif math.isclose(intersect, 0,
                                      abs_tol=1e-6) or math.isclose(intersect,
                                                                    volmdlr.two_pi,
                                                                    abs_tol=1e-6) or (
                            not check1 or not check2):
                        all_points = volmdlr.check_singularity(all_points)

                        points_cleaned = volmdlr.delete_double_point(all_points)
                        all_points = [pt.copy() for pt in points_cleaned]
                        all_points.sort(key=lambda pt: pt[0])
                        d1, d2 = (all_points[0] - all_points[-1]).Norm(), (
                                    all_points[0] - all_points[-2]).Norm()
                        if d2 < d1:
                            last = all_points[-1].copy()
                            all_points[-1] = all_points[-2].copy()
                            all_points[-2] = last
                        break
                    else:
                        points = []
                        for list_prim in primitives:
                            for k, prim in enumerate(list_prim):
                                new_list_points = []
                                change = 0
                                for pt in prim.points:
                                    if pt[0] < intersect:
                                        change += 1
                                        if math.isclose(pt[0], 0,
                                                        abs_tol=1e-1):
                                            new_list_points.append(volmdlr.volmdlr.Point2D((
                                                                           intersect + volmdlr.two_pi,
                                                                           pt[
                                                                               1])))
                                        else:
                                            new_list_points.append(volmdlr.volmdlr.Point2D(
                                                (volmdlr.two_pi + pt[0], pt[1])))
                                    elif math.isclose(pt[0], intersect,
                                                      abs_tol=1e-1):
                                        change += 1
                                        new_list_points.append(volmdlr.volmdlr.Point2D(
                                            (volmdlr.two_pi + pt[0], pt[1])))
                                    else:
                                        new_list_points.append(pt)
                                if change > 0:
                                    points.extend(new_list_points)
                                    # list_prim[k] = volmdlr.LineSegment2D(new_list_points[0], new_list_points[1])
                                else:
                                    points.extend(prim.points)
                                    continue
                        points_cleaned = volmdlr.delete_double_point(points)
                        all_points = volmdlr.faces3d.Face3D.range_trigo(points_cleaned)
                    solve = True
                else:
                    points_cleaned = volmdlr.delete_double_point(all_points)
                    all_points = [pt.copy() for pt in points_cleaned]
                    all_points.sort(key=lambda pt: pt[0])
                    d1, d2 = (all_points[0] - all_points[-1]).Norm(), (
                                all_points[0] - all_points[-2]).Norm()
                    if d2 < d1:
                        last = all_points[-1].copy()
                        all_points[-1] = all_points[-2].copy()
                        all_points[-2] = last
        else:
            points_cleaned = volmdlr.delete_double_point(all_points)
            all_points = [pt.copy() for pt in points_cleaned]
            all_points.sort(key=lambda pt: pt[0])
            d1, d2 = (all_points[0] - all_points[-1]).Norm(), (
                        all_points[0] - all_points[-2]).Norm()
            if d2 < d1:
                last = all_points[-1].copy()
                all_points[-1] = all_points[-2].copy()
                all_points[-2] = last

        primitives = volmdlr.faces3d.Face3D.create_primitives(all_points)

        # fig, ax = plt.subplots()
        # [pt.MPLPlot(ax=ax, color='g') for pt in all_points]
        # all_points[0].MPLPlot(ax=ax, color='m')
        # all_points[1].MPLPlot(ax=ax, color='r')
        # all_points[-1].MPLPlot(ax=ax)
        # all_points[-2].MPLPlot(ax=ax, color='b')
        # [p.MPLPlot(ax=ax) for p in primitives]

        l_vert = volmdlr.LineSegment2D((pt2 + pt4) / 2, (pt1 + pt3) / 2)
        solve = False
        for prim in primitives:
            if solve:
                break
            intersect = prim.line_intersections(l_vert)
            if intersect is not None:
                x_intersect = intersect.vector[0]
                y_intersect = intersect.vector[1]
                value1, value2 = ymax - 0.2 * (ymax - ymin), ymin + 0.2 * (
                            ymax - ymin)
                if y_intersect < max(value1, value2) and y_intersect > min(
                        value1, value2):
                    points = []
                    for k, prim in enumerate(primitives):
                        new_list_points, change = [], 0
                        for pt in prim.points:
                            if pt[0] < x_intersect:
                                change += 1
                                if math.isclose(pt[0], 0, abs_tol=1e-1):
                                    new_list_points.append(volmdlr.volmdlr.Point2D(
                                        (x_intersect + volmdlr.two_pi, pt[1])))
                                else:
                                    new_list_points.append(
                                        volmdlr.volmdlr.Point2D((volmdlr.two_pi + pt[0], pt[1])))
                            else:
                                new_list_points.append(pt)
                        if change > 0:
                            points.extend(new_list_points)
                            primitives[k] = volmdlr.LineSegment2D(new_list_points[0],
                                                          new_list_points[1])
                        else:
                            points.extend(prim.points)
                            continue
                    solve = True
                    points_cleaned = volmdlr.delete_double_point(points)
                    all_points = volmdlr.faces3d.Face3D.range_trigo(points_cleaned)
                    primitives = volmdlr.faces3d.Face3D.create_primitives(all_points)

        contour2d = [volmdlr.Contour2D(primitives)]
        return contour2d


class ConicalSurface3D(Surface3D):
    face_class = volmdlr.faces3d.ConicalFace3D
    """
    :param frame: Cone's frame to position it: frame.w is axis of cone
    :type frame: volmdlr.Frame3D
    :param r: Cone's bottom radius
    :type r: float
    :param semi_angle: Cone's semi-angle
    :type semi_angle: float
    """

    def __init__(self, frame, semi_angle, name=''):
        self.frame = frame
        self.semi_angle = semi_angle
        self.name = name

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = volmdlr.Frame3D(frame3d.origin, U, V, W)
        semi_angle = float(arguments[3])
        return cls(frame_direct, semi_angle, arguments[0][1:-1])

    def frame_mapping(self, frame, side, copy=True):
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ConicalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ConicalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

    def point2d_to_3d(self, point):
        theta, z = point
        new_point = volmdlr.Point3D((z * self.semi_angle * math.cos(theta),
                             z * self.semi_angle * math.sin(theta),
                             z,
                             ))
        return self.frame.OldCoordinates(new_point)

    def point3d_to_2d(self, point):
        z = self.frame.w.Dot(point)
        x, y = point.PlaneProjection(self.frame.origin, self.frame.u,
                                     self.frame.v)
        theta = math.atan2(y, x)
        return volmdlr.volmdlr.Point2D([theta, z])

    def rectangular_cut(self, theta1: float, theta2: float,
                        z1: float, z2: float, name: str=''):
        # theta1 = angle_principal_measure(theta1)
        # theta2 = angle_principal_measure(theta2)
        if theta1 == theta2:
            theta2 += volmdlr.two_pi

        p1 = volmdlr.volmdlr.Point2D((theta1, z1))
        p2 = volmdlr.volmdlr.Point2D((theta2, z1))
        p3 = volmdlr.volmdlr.Point2D((theta2, z2))
        p4 = volmdlr.volmdlr.Point2D((theta1, z2))
        outer_contour = volmdlr.Polygon2D([p1, p2, p3, p4])
        return volmdlr.faces3d.ConicalFace3D(self, outer_contour, [], name)


class SphericalSurface3D(Surface3D):
    face_class = volmdlr.faces3d.SphericalFace3D
    """
    :param frame: Sphere's frame to position it
    :type frame: volmdlr.Frame3D
    :param radius: Sphere's radius
    :type radius: float
    """

    def __init__(self, frame, radius, name=''):
        self.frame = frame
        self.radius = radius
        self.name = name
        V = frame.v
        V.Normalize()
        W = frame.w
        W.Normalize()
        self.plane = Plane3D(frame.origin, V, W)

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, frame3d.u
        U.Normalize()
        W.Normalize()
        V = W.Cross(U)
        frame_direct = volmdlr.Frame3D(frame3d.origin, U, V, W)
        radius = float(arguments[2]) / 1000
        return cls(frame_direct, radius, arguments[0][1:-1])


    def point2d_to3d(self, point2d):
        # source mathcurve.com/surfaces/sphere
        # -pi<theta<pi, -pi/2<phi<pi/2
        theta, phi = point2d
        x = self.radius * math.cos(phi) * math.cos(theta)
        y = self.radius * math.cos(phi) * math.sin(theta)
        z = self.radius * math.sin(phi)
        return self.frame3d.OldCoordinates(volmdlr.Point3D([x, y, z]))

    def point3d_to2d(self, point3d):
        x, y, z = point3d
        if z < -self.radius:
            z = -self.radius
        elif z > self.radius:
            z = self.radius

        zr = z / self.radius
        phi = math.asin(zr)

        u = math.sqrt((self.radius ** 2) - (z ** 2))
        if u == 0:
            u1, u2 = x, y
        else:
            u1, u2 = round(x / u, 5), round(y / u, 5)
        theta = volmdlr.sin_cos_angle(u1, u2)
        return volmdlr.volmdlr.Point2D((theta, phi))

class BSplineSurface3D(Surface3D):
    def __init__(self, degree_u, degree_v, control_points, nb_u, nb_v,
                 u_multiplicities, v_multiplicities, u_knots, v_knots,
                 weights=None, name=''):
        volmdlr.Primitive3D.__init__(self, basis_primitives=control_points, name=name)
        self.control_points = control_points
        self.degree_u = degree_u
        self.degree_v = degree_v
        self.nb_u = nb_u
        self.nb_v = nb_v
        u_knots = standardize_knot_vector(u_knots)
        v_knots = standardize_knot_vector(v_knots)
        self.u_knots = u_knots
        self.v_knots = v_knots
        self.u_multiplicities = u_multiplicities
        self.v_multiplicities = v_multiplicities
        self.weights = weights

        self.control_points_table = []
        points_row = []
        i = 1
        for pt in control_points:
            points_row.append(pt)
            if i == nb_v:
                self.control_points_table.append(points_row)
                points_row = []
                i = 1
            else:
                i += 1
        surface = BSpline.Surface()
        surface.degree_u = degree_u
        surface.degree_v = degree_v
        if weights is None:
            P = [(control_points[i][0], control_points[i][1],
                  control_points[i][2]) for i in range(len(control_points))]
            surface.set_ctrlpts(P, nb_u, nb_v)
        else:
            Pw = [(control_points[i][0] * weights[i],
                   control_points[i][1] * weights[i],
                   control_points[i][2] * weights[i],
                   weights[i]) for i in range(len(control_points))]
            surface.set_ctrlpts(Pw, nb_u, nb_v)
        knot_vector_u = []
        for i, u_knot in enumerate(u_knots):
            knot_vector_u.extend([u_knot] * u_multiplicities[i])
        knot_vector_v = []
        for i, v_knot in enumerate(v_knots):
            knot_vector_v.extend([v_knot] * v_multiplicities[i])
        surface.knotvector_u = knot_vector_u
        surface.knotvector_v = knot_vector_v
        surface.delta = 0.05
        surface_points = surface.evalpts

        self.surface = surface
        self.points = [volmdlr.Point3D((p[0], p[1], p[2])) for p in surface_points]

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        script = ""
        points = '['
        for i, pts_row in enumerate(self.control_points_table):
            pts = '['
            for j, pt in enumerate(pts_row):
                point = 'fc.Vector({},{},{}),'.format(pt[0], pt[1], pt[2])
                pts += point
            pts = pts[:-1] + '],'
            points += pts
        points = points[:-1] + ']'

        script += '{} = Part.BSplineSurface()\n'.format(name)
        if self.weights is None:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={})\n'.format(
                name, points, self.u_multiplicities, self.v_multiplicities,
                self.degree_u, self.degree_v, self.u_knots, self.v_knots)
        else:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={},weights={})\n'.format(
                name, points, self.u_multiplicities, self.v_multiplicities,
                self.degree_u, self.degree_v, self.u_knots, self.v_knots,
                self.weights)

        return script

    def Rotation(self, center, axis, angle, copy=True):
        new_control_points = [p.Rotation(center, axis, angle, True) for p in
                              self.control_points]
        new_BSplineSurface3D = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)
        if copy:
            return new_BSplineSurface3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineSurface3D.curve
            self.points = new_BSplineSurface3D.points

    def Translation(self, offset, copy=True):
        new_control_points = [p.Translation(offset, True) for p in
                              self.control_points]
        new_BSplineSurface3D = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)
        if copy:
            return new_BSplineSurface3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineSurface3D.curve
            self.points = new_BSplineSurface3D.points

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]

        degree_u = int(arguments[1])
        degree_v = int(arguments[2])
        points_sets = arguments[3][1:-1].split("),")
        points_sets = [elem + ")" for elem in points_sets[:-1]] + [
            points_sets[-1]]
        control_points = []
        for points_set in points_sets:
            points = [object_dict[int(i[1:])] for i in
                      points_set[1:-1].split(",")]
            nb_v = len(points)
            control_points.extend(points)
        nb_u = int(len(control_points) / nb_v)
        surface_form = arguments[4]
        if arguments[5] == '.F.':
            u_closed = False
        elif arguments[5] == '.T.':
            u_closed = True
        else:
            raise ValueError
        if arguments[6] == '.F.':
            v_closed = False
        elif arguments[6] == '.T.':
            v_closed = True
        else:
            raise ValueError
        self_intersect = arguments[7]
        u_multiplicities = [int(i) for i in arguments[8][1:-1].split(",")]
        v_multiplicities = [int(i) for i in arguments[9][1:-1].split(",")]
        u_knots = [float(i) for i in arguments[10][1:-1].split(",")]
        v_knots = [float(i) for i in arguments[11][1:-1].split(",")]
        knot_spec = arguments[12]

        if 13 in range(len(arguments)):
            weight_data = [float(i) for i in
                           arguments[13][1:-1].replace("(", "").replace(")",
                                                                        "").split(
                               ",")]
        else:
            weight_data = None

        return cls(degree_u, degree_v, control_points, nb_u, nb_v,
                   u_multiplicities, v_multiplicities, u_knots, v_knots,
                   weight_data, name)
