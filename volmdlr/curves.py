"""
Volmdlr curves.

"""
import math
import sys
from typing import List, Union

import matplotlib.pyplot as plt
import numpy as npy
import scipy.integrate as scipy_integrate
from matplotlib import __version__ as _mpl_version
from packaging import version

from dessia_common.core import DessiaObject

import plot_data.colors
import plot_data.core as plot_data
import volmdlr
from volmdlr import core, geometry
from volmdlr.nurbs.helpers import generate_knot_vector
import volmdlr.utils.common_operations as vm_common_operations
import volmdlr.utils.intersections as volmdlr_intersections
from volmdlr.core import EdgeStyle


def hyperbola_parabola_control_point_and_weight(start, start_tangent, end, end_tangent, point):
    """Gets control points and weights for hyperbola and parabola curves represented by bsplines."""
    start_tangent = start_tangent.unit_vector()
    end_tangent = end_tangent.unit_vector()
    line_class = globals()["Line"+start.__class__.__name__[-2:]]
    vector_02 = end - start
    line02 = line_class.from_point_and_vector(start, vector_02)

    line0 = line_class.from_point_and_vector(start, start_tangent)
    line2 = line_class.from_point_and_vector(end, end_tangent)
    point1 = line0.line_intersections(line2)[0]
    vector_p1 = point1 - point
    line1p = line_class.from_point_and_vector(point1, vector_p1)
    point_q = line02.line_intersections(line1p)[0]
    a = math.sqrt((start - point_q).norm()/(point_q - end).norm())
    u = a/(1.0 + a)
    num = ((1.0 - u)**2) * (point - start).dot(vector_p1) + u**2 * (point - end).dot(vector_p1)
    den = 2.0 * u * (1.0 - u) * vector_p1.dot(vector_p1)
    weight_1 = num/den
    return point1, weight_1


class Curve(DessiaObject):
    """Abstract class for a curve object."""

    def __init__(self, name: str = ''):
        self.periodic = False
        DessiaObject.__init__(self, name=name)

    def abscissa(self, point):
        """
        Calculate the abscissa of a point on the curve.
        """
        raise NotImplementedError(f'abscissa method not implemented by {self.__class__.__name__}')

    def sort_points_along_curve(self, points: List[Union[volmdlr.Point2D, volmdlr.Point3D]]):
        """
        Sort point along a curve.

        :param points: list of points to be sorted.
        :return: sorted points.
        """
        return sorted(points, key=self.abscissa)

    def line_intersections(self, line):
        """
        Calculate the line_intersections between line and curve.
        """
        raise NotImplementedError(f'line_intersections method not implemented by {self.__class__.__name__}')

    def linesegment_intersections(self, linesegment, abs_tol: float = 1e-6):
        """
        Calculates the intersections between a curve and a Line Segment.

        :param linesegment: the Line Segment.
        :param abs_tol: tolerance.
        :return:a list containing all intersections between the two objects, if any exists.
        """
        line_intersections = self.line_intersections(linesegment.line)
        intersections = []
        for intersection in line_intersections:
            if linesegment.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def curve_intersections(self, other_curve, abs_tol: float = 1e-6):
        """
        Gets the intersections between two curves.

        :param other_curve: other curve.
        :param abs_tol: tolerance.
        :return: list of intersection points.
        """
        method_name = f'{other_curve.__class__.__name__.lower()[:-2]}_intersections'
        if hasattr(self, method_name):
            intersections = getattr(self, method_name)(other_curve, abs_tol)
            return intersections
        method_name = f'{self.__class__.__name__.lower()[:-2]}_intersections'
        if hasattr(other_curve, method_name):
            intersections = getattr(other_curve, method_name)(self, abs_tol)
            return intersections
        raise NotImplementedError


class ClosedCurve(Curve):
    """Abstract class for defining closed curves (Circle, Ellipse) properties."""
    def __init__(self, name: str = ''):
        Curve.__init__(self, name=name)
        self.periodic = True

    def point_at_abscissa(self, abscissa):
        """
        Returns the point that corresponds to the given abscissa.

        :param abscissa: The abscissa
        :type abscissa: float
        :return: The point that corresponds to the given abscissa.
        :rtype: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        """
        raise NotImplementedError(f'point_at_abscissa method using abscissa'
                                  f'{abscissa} not implemented by {self.__class__.__name__}')

    def length(self):
        """
        Calculates the Closed Curve's length.
        """
        raise NotImplementedError(f'length method not implemented by {self.__class__.__name__}')

    def local_discretization(self, point1, point2, number_points: int = 10):
        """
        Gets n discretization points between two given points of the Curve.

        :param point1: point 1 on edge.
        :param point2: point 2 on edge.
        :param number_points: number of points to discretize locally.
        :return: list of locally discretized points.
        """
        abscissa1 = self.abscissa(point1)
        abscissa2 = self.abscissa(point2)
        if point1.is_close(point2) and point1.is_close(self.point_at_abscissa(0.0)):
            abscissa1 = 0.0
            abscissa2 = self.length()
            points = vm_common_operations.get_abscissa_discretization(self, abscissa1, abscissa2, number_points, False)
            return points + [points[0]]
        if abscissa1 > abscissa2 == 0.0:
            abscissa2 = self.length()
        return vm_common_operations.get_abscissa_discretization(self, abscissa1, abscissa2, number_points, False)


class Line(Curve):
    """
    Abstract class representing a line.

    :param point1: The first point defining the line
    :type point1: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    :param point2: The second point defining the line
    :type point2: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    :param name: Name of the line. Default value is an empty string
    :type name: str, optional
    """

    def __init__(self, point1, point2, name=''):
        self.point1 = point1
        self.point2 = point2
        self._direction_vector = None
        Curve.__init__(self, name=name)

    def __getitem__(self, key):
        """
        Get a point of the line by its index.
        """
        if key == 0:
            return self.point1
        if key == 1:
            return self.point2
        raise IndexError

    def unit_direction_vector(self, *args, **kwargs):
        """
        Get the unit direction vector of the line.

        :return: The unit direction vector of the line
        :rtype:  Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        vector = self.direction_vector()
        vector = vector.unit_vector()
        return vector

    def direction_vector(self, *args, **kwargs):
        """
        Get the direction vector of the line.

        :return: The direction vector of the line
        :rtype: Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        if not self._direction_vector:
            direction_vector = self.point2 - self.point1
            if isinstance(direction_vector, volmdlr.Point3D):
                direction_vector = direction_vector.to_vector()
            self._direction_vector = direction_vector
        return self._direction_vector

    def normal_vector(self, *args, **kwargs):
        """
        Get the normal vector of the line.

        :return: The normal vector of the line
        :rtype: Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        return self.direction_vector().normal_vector()

    def unit_normal_vector(self, *args, **kwargs):
        """
        Get the unit normal vector of the line.

        :return: The unit normal vector of the line
        :rtype: Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        return self.unit_direction_vector().normal_vector()

    def point_projection(self, point):
        """
        Calculate the projection of a point onto the line.

        :param point: The point to project
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :return: The projection of the point onto the line and the distance
            between the point and the projection
        :rtype: Tuple(Union[:class:`volmdlr.Point2D`,
            :class:`volmdlr.Point3D`], float)
        """
        vector = self.point2 - self.point1
        norm_u = vector.norm()
        projection_param_t = (point - self.point1).dot(vector) / norm_u ** 2
        projection = self.point1 + projection_param_t * vector
        projection = projection.to_point()
        return projection, projection_param_t * norm_u

    def abscissa(self, point):
        """
        Calculate the abscissa of a point on the line.

        :param point: The point for which to calculate the abscissa
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :return: The abscissa of the point
        :rtype: float
        """
        vector = self.point2 - self.point1
        norm_u = vector.norm()
        t_param = (point - self.point1).dot(vector) / norm_u
        return t_param

    def point_at_abscissa(self, abscissa):
        """
        Returns the point that corresponds to the given abscissa.

        :param abscissa: The abscissa
        :type abscissa: float
        :return: The point that corresponds to the given abscissa.
        :rtype: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        """
        return self.point1 + self.unit_direction_vector() * abscissa

    def split(self, split_point):
        """
        Split a line into two lines.

        :param split_point: The point where to split the line
        :type split_point: Union[:class:`volmdlr.Point2D`,
            :class:`volmdlr.Point3D`]
        :return: A list containing two lines
        """
        return [self.__class__(self.point1, split_point),
                self.__class__(split_point, self.point2)]

    def is_between_points(self, point1: Union[volmdlr.Point2D, volmdlr.Point3D],
                          point2: Union[volmdlr.Point2D, volmdlr.Point3D]):
        """
        Verifies if a line is between two points.

        :param point1: The first point
        :type point1: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param point2: The second point
        :type point2: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :return: True if the line is between the two points, False otherwise
        :rtype: bool
        """

        if point1.is_close(point2):
            return False

        line_segment = volmdlr.edges.LineSegment2D(point1, point2)
        if line_segment.line_intersections(self):
            return True
        return False

    def to_step(self, current_id, *args, **kwargs):
        """Exports to STEP format."""
        p1_content, p1_id = self.point1.to_step(current_id)
        # p2_content, p2_id = self.point2.to_step(current_id+1)
        current_id = p1_id + 1
        u_content, u_id = self.unit_direction_vector().to_step(current_id, vector=True)
        current_id = u_id + 1
        content = p1_content + u_content
        content += f"#{current_id} = LINE('{self.name}',#{p1_id},#{u_id});\n"
        return content, current_id

    def reverse(self):
        """Gets a line in the reverse direction."""
        return self.__class__(self.point2, self.point1, name=self.name + '_reverse')

    def closest_point_on_line(self, point):
        """
        Gets point on the line closest to given point.

        :param point: Other point.
        """
        segment_vector = self.direction_vector()
        p_vector = point - self.point1
        p_vector = p_vector.to_vector()
        t_param = p_vector.dot(segment_vector) / segment_vector.dot(segment_vector)
        point = self.point1 + t_param * segment_vector
        return point

    @classmethod
    def from_point_and_vector(cls, point: Union[volmdlr.Point2D, volmdlr.Point3D],
                              direction_vector: Union[volmdlr.Vector2D, volmdlr.Vector3D], name: str = ''):
        """
        Creates a Line object using only a point and a direction vector.

        :param point: line's origin point.
        :param direction_vector: line's direction vector.
        :param name: line's name.
        :return:
        """
        point2 = point + direction_vector
        return cls(point, point2, name=name)


class Line2D(Line):
    """
    Define an infinite line given by two points in 2D.

    """

    def __init__(self, point1: volmdlr.Point2D,
                 point2: volmdlr.Point2D, *, name=''):
        Line.__init__(self, point1, point2, name=name)

    def __hash__(self):
        return hash(('line2d', self.point1, self.point2))

    def to_3d(self, plane_origin, x1, x2):
        """
        Convert the line to a 3D line.

        :param plane_origin: Origin of the plane in which the line is.
        :type plane_origin: :class:`volmdlr.Point3D`
        :param x1: First direction of the plane in which the line is.
        :type x1: :class:`volmdlr.Vector3D`
        :param x2: Second direction of the plane in which the line is.
        :type x2: :class:`volmdlr.Vector3D`
        :return: The 3D line.
        :rtype: :class:`Line3D`
        """
        points_3d = [point.to_3d(plane_origin, x1, x2) for point in [self.point1, self.point2]]
        return Line3D(*points_3d, self.name)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Line2D rotation.

        :param center: rotation center.
        :param angle: angle rotation.
        :return: a new rotated Line2D.
        """
        return Line2D(*[point.rotation(center, angle)
                        for point in [self.point1, self.point2]])

    def translation(self, offset: volmdlr.Vector2D):
        """
        Line2D translation.

        :param offset: translation vector.
        :return: A new translated Line2D.
        """
        return Line2D(*[point.translation(offset) for point in [self.point1, self.point2]])

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Map the line to a new coordinate frame.

        :param frame: The new coordinate frame.
        :type frame: :class:`volmdlr.Frame2D`
        :param side: The side to which the mapping is made. 'old' for the
            original coordinate frame, 'new' for the new one.
        :type side: str
        :return: The mapped line.
        :rtype: :class:`Line2D`
        """
        return Line2D(*[point.frame_mapping(frame, side) for point in [self.point1, self.point2]])

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Plot the line.

        :param ax: Matplotlib axis on which to plot the line. If none,
            a new figure is created.
        :type ax: matplotlib.axes._subplots.AxesSubplot, optional
        :param edge_style: data class instance, containing all parameters needed to plot Line 2D.
        :return: The Matplotlib axis.
        :rtype: matplotlib.axes._subplots.AxesSubplot
        """
        if ax is None:
            _, ax = plt.subplots()

        if version.parse(_mpl_version) >= version.parse('3.3.2'):
            if edge_style.dashed:
                ax.axline((self.point1.x, self.point1.y),
                          (self.point2.x, self.point2.y),
                          dashes=[30, 5, 10, 5],
                          color=edge_style.color)
            else:
                ax.axline((self.point1.x, self.point1.y),
                          (self.point2.x, self.point2.y),
                          color=edge_style.color)
        else:
            direction_vector = self.direction_vector()
            point3 = self.point1 - 3 * direction_vector
            point4 = self.point2 + 4 * direction_vector
            if edge_style.dashed:
                ax.plot([point3[0], point4[0]], [point3[1], point4[1]], color=edge_style.color,
                        dashes=[30, 5, 10, 5])
            else:
                ax.plot([point3[0], point4[0]], [point3[1], point4[1]], color=edge_style.color)

        return ax

    def plot_data(self, edge_style=None):
        """
        Get plot data for the line.

        :param edge_style: Plotting style for the line.
        :type edge_style: :class:`plot_data.EdgeStyle`, optional
        :return: Plot data for the line.
        :rtype: :class:`plot_data.Line2D`
        """
        return plot_data.Line2D([self.point1.x, self.point1.y],
                                [self.point2.x, self.point2.y],
                                edge_style=edge_style)

    def line_intersections(self, line):
        """
        Calculate the intersection between the two lines.

        :param line: The line to calculate intersections with.
        :type line: :class:`volmdlr.Line2D`
        :return: A list of at most one intersection point between
            the two lines.
        :rtype: List[:class:`volmdlr.Point2D`]
        """

        point = volmdlr.Point2D.line_intersection(self, line)
        if point is not None:
            point_projection1, _ = self.point_projection(point)
            if point_projection1 is None:
                return []

            if line.__class__.__name__ == 'Line2D':
                point_projection2, _ = line.point_projection(point)
                if point_projection2 is None:
                    return []

            return [point_projection1]
        return []

    def linesegment_intersections(self, linesegment):
        """
        Calculate the intersection between a line and a line segment.

        :param linesegment: The line segment to calculate intersections with.
        :type linesegment: :class:`volmdlr.edges.LineSegment2D`
        :return: A list of at most one intersection point between the two lines.
        :rtype: List[:class:`volmdlr.Point2D`]
        """
        return linesegment.line_intersections(self)

    @staticmethod
    def _compute_data_create_tangent_circle(line, point, other_line):
        """
        Static helper method to compute some data used in create_tangent_circle method.
        """

        def vectors_from_line_and_point(line1, line2, point_):
            vector_i = volmdlr.Vector2D(point_.x, point_.y)
            vector_a = volmdlr.Vector2D(line1.point1.x, line1.point1.y)
            vector_b = volmdlr.Vector2D(line1.point2.x, line1.point2.y)
            vector_c = volmdlr.Vector2D(line2.point1.x, line2.point1.y)
            vector_d = volmdlr.Vector2D(line2.point2.x, line2.point2.y)
            return vector_i, vector_a, vector_b, vector_c, vector_d

        if math.isclose(line.point_distance(point), 0, abs_tol=1e-10):
            vectors = vectors_from_line_and_point(line, other_line, point)
        elif math.isclose(other_line.point_distance(point), 0, abs_tol=1e-10):
            vectors = vectors_from_line_and_point(other_line, line, point)
        else:
            raise AttributeError("The point isn't on any of the two lines")
        return vectors

    @staticmethod
    def _change_reference_frame(vector_i, vector_a, vector_b, vector_c, vector_d):
        new_u = volmdlr.Vector2D((vector_b - vector_a))
        new_u = new_u.unit_vector()
        new_v = new_u.unit_normal_vector()
        new_basis = volmdlr.Frame2D(vector_i, new_u, new_v)

        new_a = new_basis.global_to_local_coordinates(vector_a)
        new_b = new_basis.global_to_local_coordinates(vector_b)
        new_c = new_basis.global_to_local_coordinates(vector_c)
        new_d = new_basis.global_to_local_coordinates(vector_d)

        return new_basis, new_a, new_b, new_c, new_d

    @staticmethod
    def compute_tangent_circle_for_parallel_segments(new_basis, new_a, new_c):
        """
        Compute tangent circle between parallel segments.

        """
        segments_distance = abs(new_c[1] - new_a[1])
        radius = segments_distance / 2
        new_circle_center = volmdlr.Point2D((0, npy.sign(new_c[1] - new_a[1]) * radius))
        circle_center = new_basis.local_to_global_coordinates(new_circle_center)
        circle = Circle2D(circle_center, radius)
        return circle, None

    @staticmethod
    def compute_tangent_circles_for_perpendicular_segments(new_basis, new_a, new_b, new_c, new_d):
        """
        Computes tangent circle between perpendicular segments.

        """
        line_ab = Line2D(volmdlr.Point2D(new_a), volmdlr.Point2D(new_b))
        line_cd = Line2D(volmdlr.Point2D(new_c), volmdlr.Point2D(new_d))
        new_pt_k = volmdlr.Point2D.line_intersection(line_ab, line_cd)

        radius = abs(new_pt_k[0])
        new_circle_center1 = volmdlr.Point2D((0, radius))
        new_circle_center2 = volmdlr.Point2D((0, -radius))
        circle_center1 = new_basis.local_to_global_coordinates(new_circle_center1)
        circle_center2 = new_basis.local_to_global_coordinates(new_circle_center2)
        circle1 = Circle2D(circle_center1, radius)
        circle2 = Circle2D(circle_center2, radius)

        return circle1, circle2

    @staticmethod
    def get_concurrent_segments_tangent_circles(vector_i, vector_c, vector_d, new_point_k, new_basis):
        """Creates circles tangents to concurrent segments."""
        point_k = volmdlr.Point2D(new_basis.local_to_global_coordinates(new_point_k))

        if point_k.is_close(vector_i):
            return None, None

        # CHANGEMENT DE REPERE:
        new_u2 = volmdlr.Vector2D(point_k - vector_i).unit_vector()
        new_v2 = new_u2.unit_vector()
        new_basis2 = volmdlr.Frame2D(vector_i, new_u2, new_v2)
        new_vector_c = new_basis2.global_to_local_coordinates(vector_c)
        new_vector_d = new_basis2.global_to_local_coordinates(vector_d)
        new_point_k = new_basis2.global_to_local_coordinates(point_k)
        teta1 = math.atan2(new_vector_c[1], new_vector_c[0] - new_point_k[0])
        teta2 = math.atan2(new_vector_d[1], new_vector_d[0] - new_point_k[0])

        if teta1 < 0:
            teta1 += math.pi
        if teta2 < 0:
            teta2 += math.pi
        teta = teta1
        if not math.isclose(teta1, teta2, abs_tol=1e-08):
            if math.isclose(teta1, math.pi, abs_tol=1e-08) or math.isclose(
                    teta1, 0., abs_tol=1e-08):
                teta = teta2
            elif math.isclose(teta2, math.pi,
                              abs_tol=1e-08) or math.isclose(teta2, 0.,
                                                             abs_tol=1e-08):
                teta = teta1
        radius1 = new_point_k[0] * math.sin(teta) / (1 + math.cos(teta))
        radius2 = new_point_k[0] * math.sin(teta) / (1 - math.cos(teta))
        circle_center1 = new_basis2.local_to_global_coordinates(volmdlr.Point2D(0, -radius1))
        circle_center2 = new_basis2.local_to_global_coordinates(volmdlr.Point2D(0, radius2))

        if new_basis.global_to_local_coordinates(circle_center1)[1] > 0:
            return Circle2D(circle_center1, radius1), Circle2D(circle_center2, radius2)
        return Circle2D(circle_center2, radius2), Circle2D(circle_center1, radius1)

    def create_tangent_circle(self, point, other_line):
        """
        Computes the two circles that are tangent to 2 lines and intersect a point located on one of the two lines.
        """
        vector_i, vector_a, vector_b, vector_c, vector_d = self._compute_data_create_tangent_circle(
            self, point, other_line)
        # Basis change
        new_basis, new_a, new_b, new_c, new_d = self._change_reference_frame(vector_i, vector_a, vector_b,
                                                                             vector_c, vector_d)

        if new_c[1] == 0 and new_d[1] == 0:
            # Segments are on the same line: no solution
            return None, None

        if math.isclose(self.unit_direction_vector().dot(
                other_line.unit_normal_vector()), 0, abs_tol=1e-06):
            # Parallel segments: one solution
            return self.compute_tangent_circle_for_parallel_segments(new_basis, new_a, new_c)

        if math.isclose(self.unit_direction_vector().dot(
                other_line.unit_direction_vector()), 0, abs_tol=1e-06):
            # Perpendicular segments: 2 solution
            return self.compute_tangent_circles_for_perpendicular_segments(new_basis, new_a, new_b, new_c, new_d)

        # =============================================================================
        # LES SEGMENTS SONT QUELCONQUES
        #   => 2 SOLUTIONS
        # =============================================================================

        line_ab = Line2D(volmdlr.Point2D(new_a), volmdlr.Point2D(new_b))
        line_cd = Line2D(volmdlr.Point2D(new_c), volmdlr.Point2D(new_d))
        return self.get_concurrent_segments_tangent_circles(
            vector_i, vector_c, vector_d, volmdlr.Point2D.line_intersection(line_ab, line_cd), new_basis)

    def cut_between_two_points(self, point1: volmdlr.Point2D,
                               point2: volmdlr.Point2D):
        """
        Cut the line between two points to create a linesegment.

        :param point1: The first point defining the linesegment
        :type point1: :class:`volmdlr.Point2D`
        :param point2: The second point defining the linesegment
        :type point2: :class:`volmdlr.Point2D`
        :return: The created linesegment
        :rtype: :class:`volmdlr.edges.LineSegment2D`
        """
        return volmdlr.edges.LineSegment2D(point1, point2)

    def point_belongs(self, point2d, abs_tol: float = 1e-6):
        """
        Verifies if the point 2D belongs to the line.

        :param point2d: point to be verified.
        :param abs_tol: absolute tolerance to consider in calculus.
        :return: True if point belongs to line, False otherwise.
        """
        return math.isclose(self.point_distance(point2d), 0, abs_tol=abs_tol)

    def point_distance(self, point2d):
        """
        Calculate the shortest distance between a line and a point.

        :param point2d: Point to calculate distance.
        :type point2d: :class:`volmdlr.Point2D`.
        :return: Distance to point.
        :rtype: float.
        """
        vector_r = self.point1 - point2d
        vector_v = self.normal_vector()
        return abs(vector_v.dot(vector_r)) / vector_v.norm()

    def get_slope(self):
        """
        Gets the line's slope.
        """
        if self.point1.x == self.point2.x:
            return math.inf
        return (self.point2.y - self.point1.y) / (self.point2.x - self.point1.x)

    def get_y_intersection(self):
        """
        Gets the intersection of the 2D line with the Y axis.

        :return: y-intersection value.
        """
        slope = self.get_slope()
        if slope == math.inf:
            return None
        return self.point1.y - slope * self.point1.x


class Line3D(Line):
    """
    Define an infinite line passing through the 2 points.

    """
    _non_data_eq_attributes = ['name', 'basis_primitives', 'bounding_box']

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 name: str = ''):
        Line.__init__(self, point1, point2, name=name)
        self._bbox = None

    @property
    def bounding_box(self):
        """Bounding Box getter."""
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        """Bounding Box setter."""
        self._bbox = new_bounding_box

    def _bounding_box(self):
        """Calculates the Bounding box."""
        xmin = min([self.point1[0], self.point2[0]])
        xmax = max([self.point1[0], self.point2[0]])
        ymin = min([self.point1[1], self.point2[1]])
        ymax = max([self.point1[1], self.point2[1]])
        zmin = min([self.point1[2], self.point2[2]])
        zmax = max([self.point1[2], self.point2[2]])

        return core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def point_belongs(self, point3d, tol: float = 1e-6):
        """
        Verifies if a point belongs to the Line 3D.

        :param point3d: point to be verified.
        :param tol: tolerance.
        :return: returns True if point belongs to the line, and False otherwise.
        """
        if point3d.is_close(self.point1):
            return True
        return self.point_distance(point3d) < tol

    def point_distance(self, point):
        """Returns the minimal distance to a point."""
        vector1 = point - self.point1
        vector1.to_vector()
        vector2 = self.point2 - self.point1
        vector2.to_vector()
        return vector1.cross(vector2).norm() / vector2.norm()

    def line_distance(self, line2):
        """
        Calculates the distance between two Line3D.

        :param line2: other Line3D.
        :return: The distance between the two lines.
        """
        direction_vector1 = self.direction_vector()
        direction_vector2 = line2.direction_vector()
        if direction_vector1.is_colinear_to(direction_vector2):
            return direction_vector1.cross(line2.point1 - self.point1).norm() / direction_vector1.norm()
        vector = line2.point1 - self.point1
        line_distance = abs(vector.dot(direction_vector1.cross(direction_vector2))) / direction_vector1.cross(
            direction_vector2).norm()
        return line_distance

    def skew_to(self, line):
        """
        Verifies if two Line3D are skew to each other, that is, they are not parallel and never intersect.

        :param line: other line.
        :return: True if they are skew, False otherwise.
        """
        if self.direction_vector().is_colinear_to(line.direction_vector()):
            return False
        if math.isclose(self.line_distance(line), 0, abs_tol=1e-6):
            return False
        return True

    def intersection(self, line2, tol: float = 1e-6):
        """
        Calculates the intersection between two Line3D, if there is an intersection.

        :param line2: other Line3D
        :param tol: maximum tolerance.
        :return: None if there is no intersection between Lines.
        A volmdlr.Point3D if there exists an intersection.
        """
        direction_vector1 = self.direction_vector()
        direction_vector2 = line2.direction_vector()
        distance_to_line = self.line_distance(line2)
        if direction_vector1.is_colinear_to(direction_vector2) or \
                not math.isclose(distance_to_line, 0, abs_tol=tol):
            return None
        if math.isclose(distance_to_line, 0, abs_tol=tol) and \
                math.isclose(direction_vector1.dot(direction_vector2), 0, abs_tol=tol):
            projected_point, _ = self.point_projection(line2.point1)
            return projected_point
        vector = self.point1 - line2.point1
        t_coefficient = (vector.dot(direction_vector2) * direction_vector2.dot(direction_vector1) -
                         vector.dot(direction_vector1) * direction_vector2.dot(direction_vector2)) / (
                                direction_vector1.dot(direction_vector1) * direction_vector2.dot(direction_vector2) -
                                direction_vector1.dot(direction_vector2) * direction_vector2.dot(direction_vector1))
        # u_coefficient = (vector.dot(direction_vector2) + t_coefficient * direction_vector1.dot(
        # direction_vector2)) / direction_vector2.dot(direction_vector2)
        intersection = self.point1 + t_coefficient * direction_vector1
        return intersection

    def line_intersections(self, line):
        """
        Gets the intersection between two Line3D, if there is an intersection.

        :param line: other Line3D
        :return: None if there is no intersection between Lines.
        A volmdlr.Point3D if there exists an intersection.
        """
        return [self.intersection(line)]

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plot method for Line 3D using Matplotlib."""
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        # Line segment
        ax.plot([self.point1.x, self.point2.x], [self.point1.y, self.point2.y],
                [self.point1.z, self.point2.z], color=edge_style.color, alpha=edge_style.alpha)

        # Drawing 3 times length of segment on each side
        u = self.point2 - self.point1
        v1 = self.point1 - u * 3
        x1, y1, z1 = v1.x, v1.y, v1.z
        v2 = self.point2 - u * 3
        x2, y2, z2 = v2.x, v2.y, v2.z
        if edge_style.dashed:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=edge_style.color,
                    dashes=[30, 5, 10, 5])
        else:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=edge_style.color)
        return ax

    def plane_projection2d(self, center, x, y):
        """
        Project the 3D line onto a 2D plane defined by the center point and two orthogonal vectors, x and y.

        :param center: The center point of the plane.
        :param x: A tuple representing the first orthogonal vector (x-component, y-component, z-component).
        :param y: A tuple representing the second orthogonal vector (x-component, y-component, z-component).

        :return: A new 2D line resulting from the projection of the current 3D line onto the specified plane.
        """
        return Line2D(self.point1.plane_projection2d(center, x, y),
                      self.point2.plane_projection2d(center, x, y))

    def minimum_distance_points(self, other_line):
        """
        Returns the points on this line and the other line that are the closest of lines.
        """
        if self.point_belongs(other_line.point1):
            return other_line.point1, other_line.point1
        if self.point_belongs(other_line.point2):
            return other_line.point2, other_line.point2
        u = self.point2 - self.point1
        v = other_line.point2 - other_line.point1
        w = self.point1 - other_line.point1
        u_dot_u = u.dot(u)
        u_dot_v = u.dot(v)
        v_dot_v = v.dot(v)
        u_dot_w = u.dot(w)
        v_dot_w = v.dot(w)

        s_param = (u_dot_v * v_dot_w - v_dot_v * u_dot_w) / (u_dot_u * v_dot_v - u_dot_v ** 2)
        t_param = (u_dot_u * v_dot_w - u_dot_v * u_dot_w) / (u_dot_u * v_dot_v - u_dot_v ** 2)
        point1 = self.point1 + s_param * u
        point2 = other_line.point1 + t_param * v
        return point1, point2

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Line3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Line3D
        """

        return Line3D(*[point.rotation(center, axis, angle) for point in
                        [self.point1, self.point2]])

    def translation(self, offset: volmdlr.Vector3D):
        """
        Line3D translation.

        :param offset: translation vector
        :return: A new translated Line3D
        """
        return Line3D(*[point.translation(offset) for point in
                        [self.point1, self.point2]])

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and return a new Line3D.

        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.local_to_global_coordinates(self.point1)
            new_end = frame.local_to_global_coordinates(self.point2)
        elif side == 'new':
            new_start = frame.global_to_local_coordinates(self.point1)
            new_end = frame.global_to_local_coordinates(self.point2)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        return Line3D(new_start, new_end)

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D, **kwargs):
        """
        Trims a line creating a line segment.

        :param point1: line segment start.
        :param point2: line segment end.
        :return: line segment.
        """
        if not self.point_belongs(point1) or not self.point_belongs(point2):
            print('Point not on curve')
        #     raise ValueError('Point not on curve')

        return volmdlr.edges.LineSegment3D(point1, point2)

    def copy(self, *args, **kwargs):
        """Creates a Copy of Line3D and returns it."""
        return Line3D(*[point.copy() for point in [self.point1, self.point2]])

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to an Line3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding Line3D object
        :rtype: :class:`Line3D`
        """
        point1 = object_dict[arguments[1]]
        direction = object_dict[arguments[2]]
        point2 = point1 + direction
        return cls(point1, point2, arguments[0][1:-1])

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Line3D into an Line2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Line2D.
        """
        p2d = [point.to_2d(plane_origin, x, y) for point in (self.point1, self.point2)]
        if p2d[0] == p2d[1]:
            return None
        return Line2D(*p2d, name=self.name)


class CircleMixin:
    """Circle abstract class."""

    def split_at_abscissa(self, abscissa):
        """
        Splits a Circle into two at a given fraction of its length (abscissa parameter).

        :param abscissa: The fraction of the circle length at which to perform the split.
                Value should be between 0.0 and circle.length(), where 0.0 represents the start of the circle and
                circle.length() represents the end of the arc.
        :type abscissa: float.

        :return: A list containing the two split Arc objects.
        :rtype: List[Arc].
        :raises: ValueError - If the abscissa value is outside the valid range [0.0, circle length].

        """
        if abscissa == 0.0:
            fullarc_class_ = getattr(volmdlr.edges, "FullArc" + self.__class__.__name__[-2:])
            return [fullarc_class_.from_curve(self)]
        start = self.point_at_abscissa(0.0)
        point_at_absccissa = self.point_at_abscissa(abscissa)
        return self.split(start, point_at_absccissa)

    def trim(self, point1: Union[volmdlr.Point2D, volmdlr.Point3D], point2: Union[volmdlr.Point2D, volmdlr.Point3D],
             same_sense: bool = True):
        """
        Trims a circle between two points.

        :param point1: point 1 used to trim circle.
        :param point2: point2 used to trim circle.
        :param same_sense: Used for periodical curves only. Indicates whether the curve direction agrees with (True)
            or is in the opposite direction (False) to the edge direction. By default, it's assumed True
        :return: arc between these two points.
        """
        fullar_arc_class_ = getattr(volmdlr.edges, 'FullArc' + self.__class__.__name__[-2:])
        arc_class_ = getattr(volmdlr.edges, 'Arc' + self.__class__.__name__[-2:])
        circle = self
        if not same_sense:
            circle = self.reverse()
        if not self.point_belongs(point1, 1e-5):
            angle = circle.get_arc_point_angle(point1)
            point1 = circle.point_at_abscissa(angle * self.radius)
        if not self.point_belongs(point2, 1e-5):
            angle = circle.get_arc_point_angle(point2)
            point2 = circle.point_at_abscissa(angle * self.radius)
        if point1.is_close(point2):
            return fullar_arc_class_(circle, point1)
        return arc_class_(circle, point1, point2)


class Circle2D(CircleMixin, ClosedCurve):
    """
    A class representing a 2D circle.

    This class inherits from `CircleMixin` and `Curve` classes,
    and provides methods to work with 2D circles.

    :param center: The center point of the circle.
    :type center: volmdlr.Point2D
    :param radius: The radius of the circle.
    :type radius: float.
    :param name: The name of the circle. Defaults to ''.
    :type name: str, optional
    """

    def __init__(self, center: volmdlr.Point2D, radius: float, name: str = ''):
        self.center = center
        self.radius = radius
        self._bounding_rectangle = None
        self.frame = volmdlr.Frame2D(center, volmdlr.X2D, volmdlr.Y2D)
        ClosedCurve.__init__(self, name=name)

    def __hash__(self):
        return int(round(1e6 * (self.center.x + self.center.y + self.radius)))

    def __eq__(self, other_circle):
        if self.__class__.__name__ != other_circle.__class__.__name__:
            return False

        return math.isclose(self.center.x,
                            other_circle.center.x, abs_tol=1e-06) \
            and math.isclose(self.center.y,
                             other_circle.center.y, abs_tol=1e-06) \
            and math.isclose(self.radius, other_circle.radius,
                             abs_tol=1e-06)

    def __getitem__(self, key):
        if key == 0:
            return self.center
        if key == 1:
            return self.radius
        raise IndexError

    @classmethod
    def from_3_points(cls, point1, point2, point3, name: str = ''):
        """
        Creates a circle 2d from 3 points.

        :return: circle 2d.
        """
        x_interior, y_interior = point2.x, point2.y
        x_end, y_end = point3.x, point3.y
        x_start, y_start = point1.x, point1.y
        matrix1 = [[2 * (x_start - x_interior), 2 * (y_start - y_interior)],
                   [2 * (x_start - x_end), 2 * (y_start - y_end)]]
        b_vector_components = [x_interior ** 2 + y_interior ** 2 - x_start ** 2 - y_start ** 2,
                               x_end ** 2 + y_end ** 2 - x_start ** 2 - y_start ** 2]
        try:
            matrix_a = volmdlr.Matrix22(*matrix1[0], *matrix1[1])
            b_vector = - volmdlr.Vector2D(*b_vector_components)
            inv_matrix_a = matrix_a.inverse()
            center = volmdlr.Point2D(*inv_matrix_a.vector_multiplication(b_vector))
        except ValueError:
            matrix_a = npy.array(matrix1)
            b_vector = - npy.array(b_vector_components)
            center = volmdlr.Point2D(*npy.linalg.solve(matrix_a, b_vector))
        circle = cls(center, point1.point_distance(center), name=name)
        return circle

    def area(self):
        """
        Calculates the area for a circle 2d.

        :return: circle area.
        """
        return math.pi * self.radius ** 2

    def second_moment_area(self, point):
        """Second moment area of part of disk."""
        sma = math.pi * self.radius ** 4 / 4
        return geometry.huygens2d(sma, sma, 0, self.area(), self.center, point)

    def center_of_mass(self):
        """Gets the circle's center of mass."""
        return self.center

    def length(self):
        """
        Calculates the length of the Circle 2D.

        :return: the circle's length.
        """

        return volmdlr.TWO_PI * self.radius

    def point_symmetric(self, point):
        """
        Creates a circle symmetrically from a point.

        :param point: symmetry point.
        :return: Circle 2D symmetric to point.
        """
        center = 2 * point - self.center
        return Circle2D(center, self.radius)

    def axial_symmetry(self, line):
        """
        Finds out the symmetric circle 2d according to a line.
        """
        return self.__class__(center=self.center.axial_symmetry(line),
                              radius=self.radius)

    def copy(self, *args, **kwargs):
        """
        Create a copy of the arc 2d.

        :return: copied circle 2d.
        """
        return Circle2D(self.center.copy(), self.radius)

    def point_at_abscissa(self, curvilinear_abscissa):
        """
        Gets the point at a given abscissa.

        :param curvilinear_abscissa: a portion of the circle's length - (0, length).
        :return: Point found at given abscissa.
        """
        start = self.center + self.radius * volmdlr.X3D
        return start.rotation(self.center, curvilinear_abscissa / self.radius)

    def abscissa(self, point: volmdlr.Point2D, tol=1e-6):
        """
        Returns the abscissa of a given point 2d.

        """
        if not math.isclose(point.point_distance(self.center), self.radius, abs_tol=tol):
            raise ValueError('Point not in arc')
        u1, u2 = point.x / self.radius, point.y / self.radius
        point_angle = geometry.sin_cos_angle(u1, u2)
        return self.radius * point_angle

    def point_belongs(self, point, include_edge_points: bool = True, tol: float = 1e-6):
        """
        Verifies if a point is inside the Circle 2D.

        :param point: A 2D point to check if it is inside the Circle 2D.
        :type point: `volmdlr.Point2D`
        :param include_edge_points: A Boolean indicating whether points on the edge of the Circle 2D
            should be considered inside the circle.
        :type include_edge_points: bool
        :param tol: tolerance.
        :return: True if point inside the circle or false otherwise.
        :rtype: bool
        """

        if include_edge_points:
            return point.point_distance(self.center) <= self.radius + tol
        return point.point_distance(self.center) < self.radius

    def point_distance(self, point):
        """
        Calculates the distance of given point to the circle.

        :param point: point to calculate distance.
        :return: the distance from the point to the circle 2D.
        """
        return abs(point.point_distance(self.center) - self.radius)

    @property
    def bounding_rectangle(self):
        """
        Gets the bounding rectangle for the circle.

        :return: bounding rectangle.
        """
        if not self._bounding_rectangle:
            self._bounding_rectangle = self.get_bounding_rectangle()
        return self._bounding_rectangle

    def get_bounding_rectangle(self):
        """
        Calculates the circle's bounding rectangle.

        :return: returns a Bounding Rectangle object.
        """
        x_min = self.center.x - self.radius
        x_max = self.center.x + self.radius
        y_min = self.center.y - self.radius
        y_max = self.center.y + self.radius
        return core.BoundingRectangle(x_min, x_max, y_min, y_max)

    def cut_by_line(self, line: Line2D):
        """
        Cuts a circle by a line and returns the resulting contours.

        :param line: The line used to cut the circle.
        :type line: (Line2D)
        :return: A list containing the resulting contours after the cut.
        :rtype: List[Union[self, Contour2D]]
        :raises: NotImplementedError - If there is only one intersection point, the method is not implemented.
                 ValueError: If there are more than two intersection points, the input is invalid.
        """
        intersection_points = self.line_intersections(line)
        if not intersection_points:
            return [self]
        if len(intersection_points) == 1:
            raise NotImplementedError
        if len(intersection_points) == 2:
            linesegment = volmdlr.edges.LineSegment2D(intersection_points[0],
                                                      intersection_points[1])
            arc1, arc2 = self.split(intersection_points[0],
                                    intersection_points[1])
            # from volmdlr import wires
            contour1 = volmdlr.wires.Contour2D([arc1, linesegment.copy()])
            contour2 = volmdlr.wires.Contour2D([arc2, linesegment.copy()])
            return [contour1, contour2]
        raise ValueError

    def line_intersections(self, line2d: Line2D, tol=1e-9):
        """
        Calculates the intersections between a circle 2D and Line 2D.

        :param line2d: line to calculate intersections
        :param tol: tolerance to consider in calculations.
        :return: circle and line intersections.
        """
        if line2d.point_distance(self.center) > self.radius + tol:
            return []
        if line2d.point_belongs(self.center):
            direction_vector = line2d.unit_direction_vector()
            return [self.center + self.radius * direction_vector, self.center - self.radius * direction_vector]
        if not self.center.is_close(volmdlr.O2D):
            local_line = line2d.frame_mapping(self.frame, 'new')
            local_circle = self.frame_mapping(self.frame, 'new')
            local_line_intersections = local_circle.line_intersections(local_line)
            return [self.frame.local_to_global_coordinates(point) for point in local_line_intersections]
        m = line2d.get_slope()
        c = line2d.get_y_intersection()
        if m == math.inf and c is None:
            x_line = line2d.point1.x
            y1 = - math.sqrt(self.radius**2 - x_line**2)
            y2 = math.sqrt(self.radius**2 - x_line**2)
            return [volmdlr.Point2D(x_line, y1), volmdlr.Point2D(x_line, y2)]
        quad_eq_a = 1 + m ** 2
        quad_eq_b = 2 * m * c
        quad_eq_c = c ** 2 - self.radius ** 2
        delta = quad_eq_b ** 2 - 4 * quad_eq_a * quad_eq_c
        if delta < 0.0 or quad_eq_a == 0.0:
            return []
        if math.isclose(delta, 0, abs_tol=1e-6):
            x1 = - quad_eq_b / 2 * quad_eq_a
            y1 = m * x1 + c
            return [volmdlr.Point2D(x1, y1)]
        x1 = (-quad_eq_b + math.sqrt(delta)) / (2 * quad_eq_a)
        x2 = (-quad_eq_b - math.sqrt(delta)) / (2 * quad_eq_a)
        y1 = m * x1 + c
        y2 = m * x2 + c
        return [volmdlr.Point2D(x1, y1), volmdlr.Point2D(x2, y2)]

    def linesegment_intersections(self, linesegment: 'volmdlr.edges.LineSegment2D', tol=1e-9):
        """
        Calculates the intersections between a circle 2D and line segment 2D.

        :param linesegment: line segment to calculate intersections
        :param tol: tolerance to consider in calculations.
        :return: circle and line segment intersections.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(linesegment.bounding_rectangle) > tol:
            return []
        line_intersections = self.line_intersections(linesegment.line, tol)
        linesegment_intersections = []
        for intersection in line_intersections:
            if linesegment.point_belongs(intersection):
                linesegment_intersections.append(intersection)
        return linesegment_intersections

    def circle_intersections(self, circle: 'Circle2D'):
        """
        Finds the intersection points between this circle and another circle.

        :param circle: The other circle to find intersections with.
        :type circle: (Circle2D).
        :return: A list of intersection points between the two circles.
        :rtype: List[Point2D].
        """
        return volmdlr_intersections.get_circle_intersections(self, circle)

    def arc_intersections(self, arc2d: 'volmdlr.edges.Arc2D', abs_tol: float = 1e-6):
        """
        Finds the intersection points between this circle and an arc 2d.

        :param arc2d: The arc 2d to find intersections with.
        :type arc2d: (edges.Arc2D).
        :param abs_tol: tolerance to be considered while validating an intersection.
        :return: A list of intersection points between the circle and the arc.
        :rtype: List[Point2D].
        """
        circle_intesections = self.circle_intersections(arc2d.circle)
        intersections = []
        for inter in circle_intesections:
            if arc2d.point_belongs(inter, abs_tol):
                intersections.append(inter)
        return intersections

    def ellipse_intersections(self, ellipse2d, abs_tol: float = 1e-7):
        """
        Finds the intersection points between this circle and an arc 2d.

        :param ellipse2d: The Ellipse 2d to find intersections with.
        :type ellipse2d: (Ellipse2D).
        :param abs_tol: Tolerance.
        :return: A list of intersection points between the circle and the arc.
        :rtype: List[Point2D].
        """
        if self.bounding_rectangle.distance_to_b_rectangle(ellipse2d.bounding_rectangle) > abs_tol:
            return []
        intersections = volmdlr_intersections.get_bsplinecurve_intersections(ellipse2d, self, abs_tol)
        return intersections

    def bsplinecurve_intersections(self, bsplinecurve: 'volmdlr.edges.BSplineCurve2D', abs_tol: float = 1e-6):
        """
        Calculates the intersections between a circle 2d and a BSpline Curve 2D.

        :param bsplinecurve: bsplinecurve to search for intersections.
        :param abs_tol: tolerance to be considered while validating an intersection.
        :return: a list with all intersections between circle and bsplinecurve.
        """
        return volmdlr_intersections.get_bsplinecurve_intersections(self, bsplinecurve, abs_tol)

    def hyperbola_intersections(self, hyperbola2d, abs_tol: float = 1e-6):
        """
        Calculates the intersections between a circle 2d and a Hyperbola 2D.

        :param hyperbola2d: hyperbola to search for intersections with.
        :param abs_tol: tolerance to be considered while validating an intersection.
        :return: a list with all intersections between circle and hyperbola.
        """
        b_rectangle = self.bounding_rectangle
        hyperbola_point1 = volmdlr.Point2D(hyperbola2d.get_x(b_rectangle.ymin), b_rectangle.ymin)
        hyperbola_point2 = volmdlr.Point2D(hyperbola2d.get_x(b_rectangle.ymax), b_rectangle.ymax)
        hyperbola_bspline = hyperbola2d.trim(hyperbola_point1, hyperbola_point2)
        return self.bsplinecurve_intersections(hyperbola_bspline, abs_tol)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plots the circle using Matplotlib."""
        return vm_common_operations.plot_circle(self, ax, edge_style)

    def plot_data(self, edge_style: plot_data.EdgeStyle = None, surface_style: plot_data.SurfaceStyle = None):
        """
        Get plot data for the circle 2d.

        :param edge_style: Plotting style for the line.
        :type edge_style: :class:`plot_data.EdgeStyle`, optional
        :return: Plot data for the line.
        :rtype: :class:`plot_data.Circle2D`
        """
        return plot_data.Circle2D(cx=self.center.x, cy=self.center.y,
                                  r=self.radius,
                                  edge_style=edge_style,
                                  surface_style=surface_style)

    def to_3d(self, plane_origin, x, y):
        """
        Transforms a Circle2D into an Circle3D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Circle3D.
        """
        normal = x.cross(y)
        center3d = self.center.to_3d(plane_origin, x, y)
        return Circle3D(volmdlr.Frame3D(center3d, x, y, normal), self.radius, self.name)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Circle2D rotation.

        :param center: rotation center.
        :param angle: angle rotation.
        :return: a new rotated Circle2D.
        """
        return Circle2D(self.center.rotation(center, angle), self.radius)

    def translation(self, offset: volmdlr.Vector2D):
        """
        Circle2D translation.

        :param offset: translation vector
        :return: A new translated Circle2D
        """
        return Circle2D(self.center.translation(offset), self.radius)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Circle2D.

        side = 'old' or 'new'
        """
        if side == 'old':
            return Circle2D(frame.local_to_global_coordinates(self.center),
                            self.radius)
        if side == 'new':
            return Circle2D(frame.global_to_local_coordinates(self.center),
                            self.radius)
        raise ValueError('Side should be \'new\' \'old\'')

    def split_by_line(self, line: Line2D):
        """
        Split the Circle with a line into two Arc2D.
        """
        split_points = self.line_intersections(line)
        return self.split(split_points[0], split_points[1])

    def split(self, split_start, split_end):
        """
        Splits the current object into two Arc2D edges.

        This method creates two Arc2D edges by splitting the current object between the specified start and end points.
        The new Arc2D edges will connect the split_start and split_end points, and split_end and split_start points
        respectively.

        :param (Point2D) split_start: The starting point of the split.
        :param (Point2D) split_end: The ending point of the split.

        :return: A list containing the two newly created Arc2D edges resulting from the split.
        """
        return [volmdlr.edges.Arc2D(self, split_start, split_end),
                volmdlr.edges.Arc2D(self, split_end, split_start)]

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 40):
        """
        Discretize a Contour to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc
        :return: a list of sampled points
        """
        if not number_points and angle_resolution:
            number_points = math.ceil(math.pi * angle_resolution) + 2
        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step) for i in range(number_points)]

    def get_geo_points(self):
        """
        Represents the circle in 3D space.
        """
        return [volmdlr.Point3D(self.radius, self.center.y, 0),
                volmdlr.Point3D(self.center.x, self.center.y, 0),
                volmdlr.Point3D(-self.radius, self.center.y, 0)]


class Circle3D(CircleMixin, ClosedCurve):
    """
    Defines a Circle in three dimensions, with a center and a radius.

    frame.u, frame.v define the plane, frame.w the normal
    """
    _non_serializable_attributes = ['point', 'edges', 'point_inside_contour']
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, frame: volmdlr.Frame3D, radius: float,
                 name: str = ''):
        self.radius = radius
        self.frame = frame
        self._bbox = None
        self.angle = 2 * math.pi
        ClosedCurve.__init__(self, name=name)

    @property
    def center(self):
        """Gets the center point of the circle 3d."""
        return self.frame.origin

    @property
    def normal(self):
        """
        Gets the circle's normal.
        """
        return self.frame.w

    def __hash__(self):
        return hash(self.frame.origin)

    def __eq__(self, other_circle):
        return self.frame.origin.is_close(other_circle.frame.origin) \
            and self.frame.w.is_colinear_to(other_circle.frame.w) \
            and math.isclose(self.radius,
                             other_circle.radius, abs_tol=1e-06)

    def __getitem__(self, key):
        if key == 0:
            return self.frame
        if key == 1:
            return self.radius
        raise IndexError

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Discretize a Circle to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc
        :return: a list of sampled points
        """
        if number_points:
            angle_resolution = number_points
        discretization_points_3d = [self.center + self.radius * math.cos(theta) * self.frame.u +
                                    self.radius * math.sin(theta) * self.frame.v for theta in
                                    npy.linspace(0, volmdlr.TWO_PI, angle_resolution)]
        return discretization_points_3d

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :param tol: tolerance.
        :return: abscissa
        """
        if not math.isclose(self.center.point_distance(point), self.radius, abs_tol=tol):
            raise ValueError('Point is not on circle')
        x, y, _ = self.frame.global_to_local_coordinates(point)
        u1 = x / self.radius
        u2 = y / self.radius
        theta = geometry.sin_cos_angle(u1, u2)

        return self.radius * abs(theta)

    def length(self):
        """Calculates the arc length of the circle."""
        return volmdlr.TWO_PI * self.radius

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Circle3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Circle3D
        """
        return Circle3D(self.frame.rotation(center, axis, angle),
                        self.radius, self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Circle3D translation.

        :param offset: translation vector
        :return: A new translated Circle3D
        """
        return Circle3D(self.frame.translation(offset), self.radius, self.name)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Circle3D.

        side = 'old' or 'new'.
        """
        return Circle3D(self.frame.frame_mapping(frame, side), self.radius)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plot method for Circle3D."""
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        return vm_common_operations.plot_from_discretization_points(ax, edge_style, self, close_plot=True)

    def point_at_abscissa(self, curvilinear_abscissa):
        """ Start point is at intersection of frame.u axis. """
        start = self.frame.origin + self.radius * self.frame.u
        return start.rotation(self.frame.origin, self.frame.w,
                              curvilinear_abscissa / self.radius)

    def line_intersections(self, line: Line3D, abs_tol: float = 1e-6):
        """
        Calculates the intersections between the Circle3D and a line 3D.

        :param line: line 3D to verify intersections
        :param abs_tol: Tolerance.
        :return: list of points intersecting Circle
        """
        circle3d_line_intersections = volmdlr_intersections.circle_3d_line_intersections(self, line, abs_tol)
        return circle3d_line_intersections

    def linesegment_intersections(self, linesegment: 'volmdlr.edges.LineSegment3D', abs_tol: float = 1e-6):
        """
        Calculates the intersections between the Circle3D and a line segment 3D.

        :param linesegment: line segment 3D to verify intersections
        :param abs_tol: tolerance to be considered while validating an intersection.
        :return: list of points intersecting Circle
        """
        intersections = []
        circle3d_line_intersections = volmdlr_intersections.circle_3d_line_intersections(self, linesegment.line)
        for intersection in circle3d_line_intersections:
            if linesegment.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def circle_intersections(self, other_circle, abs_tol: float = 1e-6):
        """
        Calculates the intersections between two Circle3D.

        :param other_circle: Circle 3D to verify intersections.
        :param abs_tol: tolerance.
        :return: list of points intersecting Circle
        """
        plane_intersections = volmdlr_intersections.get_two_planes_intersections(self.frame, other_circle.frame)
        if not plane_intersections:
            return []
        plane_intersections = Line3D(plane_intersections[0], plane_intersections[1])
        circle3d_line_intersections1 = volmdlr_intersections.circle_3d_line_intersections(self, plane_intersections)
        circle3d_line_intersections2 = volmdlr_intersections.circle_3d_line_intersections(other_circle,
                                                                                          plane_intersections)
        intersections = []
        for intersection in circle3d_line_intersections1 + circle3d_line_intersections2:
            if volmdlr.core.point_in_list(intersection, intersections):
                continue
            if self.point_belongs(intersection, abs_tol) and other_circle.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def ellipse_intersections(self, ellipse, abs_tol: float = 1e-6):
        """
        Calculates the intersections between two Circle3D.

        :param ellipse: Ellipse 3D to verify intersections.
        :param abs_tol: tolerance.
        :return: list of points intersecting Circle
        """
        intersections = []
        if self.frame.w.is_colinear_to(ellipse.frame.w) and \
                math.isclose(self.frame.w.dot(ellipse.frame.origin - self.frame.origin), 0, abs_tol=1e-6):
            ellipse2d = ellipse.to_2d(self.frame.origin, self.frame.u, self.frame.v)
            circle2d = self.to_2d(self.frame.origin, self.frame.u, self.frame.v)
            intersections_2d = circle2d.ellipse_intersections(ellipse2d)
            for intersection in intersections_2d:
                intersections.append(intersection.to_3d(self.frame.origin, self.frame.u, self.frame.v))
            return intersections

        plane_intersections = volmdlr_intersections.get_two_planes_intersections(self.frame, ellipse.frame)
        if not plane_intersections:
            return []
        plane_intersections = Line3D(plane_intersections[0], plane_intersections[1])
        circle3d_line_intersections = volmdlr_intersections.circle_3d_line_intersections(self, plane_intersections)
        ellipse3d_line_intersections = volmdlr_intersections.conic3d_line_intersections(
            ellipse, plane_intersections)
        for intersection in circle3d_line_intersections + ellipse3d_line_intersections:
            if volmdlr.core.point_in_list(intersection, intersections):
                continue
            if self.point_belongs(intersection, abs_tol) and ellipse.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Circle3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding Circle3D object.
        :rtype: :class:`volmdlr.wires.Circle3D`
        """
        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        center = object_dict[arguments[1]].origin
        radius = float(arguments[2]) * length_conversion_factor
        normal = object_dict[arguments[1]].w
        normal = normal.unit_vector()
        return cls.from_center_normal(center, normal, radius, arguments[0][1:-1])

    def to_step(self, current_id, *args, **kwargs):
        """
        Exports the circle 3d to STEP.

        """
        content, frame_id = self.frame.to_step(current_id)
        curve_id = frame_id + 1
        content += f"#{curve_id} = CIRCLE('{self.name}',#{frame_id},{self.radius * 1000});\n"
        current_id = curve_id
        # if surface_id:
        #     content += f"#{curve_id + 1} = SURFACE_CURVE('',#{curve_id},(#{surface_id}),.PCURVE_S1.);\n"
        #     curve_id += 1

        # point1 = self.frame.origin + self.frame.u * self.radius
        # point3 = self.frame.origin - self.frame.u * self.radius
        #
        # p1_content, p1_id = point1.to_step(curve_id + 1, vertex=True)
        # p3_content, p3_id = point3.to_step(p1_id + 1, vertex=True)
        # content += p1_content + p3_content
        #
        # arc1_id = p3_id + 1
        # content += f"#{arc1_id} = EDGE_CURVE('{self.name}',#{p1_id},#{p3_id},#{curve_id},.T.);\n"
        # oriented_edge1_id = arc1_id + 1
        # content += f"#{oriented_edge1_id} = ORIENTED_EDGE('',*,*,#{arc1_id},.T.);\n"
        #
        # arc2_id = oriented_edge1_id + 1
        # content += f"#{arc2_id} = EDGE_CURVE('{self.name}',#{p3_id},#{p1_id},#{curve_id},.T.);\n"
        # oriented_edge2_id = arc2_id + 1
        # content += f"#{oriented_edge2_id} = ORIENTED_EDGE('',*,*,#{arc2_id},.T.);\n"
        #
        # current_id = oriented_edge2_id + 1
        # content += f"#{current_id} = EDGE_LOOP('{self.name}',(#{oriented_edge1_id},#{oriented_edge2_id}));\n"

        return content, current_id

    @property
    def bounding_box(self):
        """Bounding box for Arc 3D."""
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    def _bounding_box(self):
        """
        Computes the bounding box.

        """
        points = [self.frame.origin + self.radius * v
                  for v in [self.frame.u, -self.frame.u,
                            self.frame.v, -self.frame.v]]
        return core.BoundingBox.from_points(points)

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Circle3D into an Circle2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Circle2D.
        """
        center = self.center.to_2d(plane_origin, x, y)
        return Circle2D(center, self.radius)

    @classmethod
    def from_center_normal(cls, center: volmdlr.Point3D,
                           normal: volmdlr.Vector3D,
                           radius: float,
                           name: str = ''):
        """Creates a Circle 3D from a center point and a normal vector, along with is radius."""
        u = normal.deterministic_unit_normal_vector()
        v = normal.cross(u)
        return cls(volmdlr.Frame3D(center, u, v, normal), radius, name)

    @classmethod
    def from_3_points(cls, point1, point2, point3, name: str = ''):
        """
        Create a Circle3D object from three points.

        This class method constructs a Circle3D object given three 3D points (Point3D objects).
        The three points are used to uniquely define a circle in 3D space.

        :param (Point3D) point1: The first point on the circumference of the circle.
        :param (Point3D) point2: The second point on the circumference of the circle.
        :param (Point3D) point3: The third point on the circumference of the circle.

        return: A Circle3D object that represents the circle uniquely defined by the three input points.

        :raise ZeroDivisionError: If the three input points are not distinct, a ZeroDivisionError is raised.
        :raise ZeroDivisionError: If the start, end, and interior points of the arc are not distinct,
        a ZeroDivisionError is raised.
        """
        # The implementation details are not described in the docstring as they are quite involved.
        # The method calculates the center, radius, and frame of the circle from the three input points in 3D space.
        # The frame represents the orientation of the circle in 3D space.
        # The method uses various geometric calculations to find these properties.
        vector_u1 = point2 - point1
        vector_u2 = point2 - point3
        try:
            vector_u1 = vector_u1.unit_vector()
            vector_u2 = vector_u2.unit_vector()
        except ZeroDivisionError as exc:
            raise ZeroDivisionError('the 3 points must be distincts') from exc

        normal = vector_u2.cross(vector_u1)
        normal = normal.unit_vector()

        if vector_u1.is_close(vector_u2):
            vector_u2 = normal.cross(vector_u1)
            vector_u2 = vector_u2.unit_vector()

        point11 = 0.5 * (point1 + point2)  # Mid-point of segment s,m
        point21 = 0.5 * (point2 + point3)  # Mid-point of segment s,m

        line1 = Line3D(point11, point11 + normal.cross(vector_u1))
        line2 = Line3D(point21, point21 + normal.cross(vector_u2))

        try:
            center, _ = line1.minimum_distance_points(line2)
        except ZeroDivisionError as exc:
            raise ZeroDivisionError('Start, end and interior points  of an arc must be distincts') from exc

        return cls(frame=volmdlr.Frame3D(center, vector_u1, normal.cross(vector_u1), normal),
                   radius=(center - point1).norm(), name=name)

    def extrusion(self, extrusion_vector):
        """
        Returns the cylindrical face generated by extrusion of the circle.
        """
        if self.normal.is_colinear_to(extrusion_vector):
            u = self.normal.deterministic_unit_normal_vector()
            v = self.normal.cross(u)
            w = extrusion_vector.copy()
            w = w.unit_vector()
            cylinder = volmdlr.surfaces.CylindricalSurface3D(
                volmdlr.Frame3D(self.center, u, v, w), self.radius)
            return [volmdlr.faces.CylindricalFace3D.from_surface_rectangular_cut(cylinder, 0, volmdlr.TWO_PI,
                                                                                 0, extrusion_vector.norm())]
        raise NotImplementedError(
            f'Extrusion along vector not colinar to normal for circle not '
            f'handled yet: dot={self.normal.dot(extrusion_vector)}')

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        """
        Return the Toroidal face generated by the revolution of the circle.
        """
        line3d = Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)
        u = self.center - tore_center
        u = u.unit_vector()
        v = axis.cross(u)
        if not math.isclose(self.normal.dot(u), 0., abs_tol=1e-9):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        tore_radius = tore_center.point_distance(self.center)
        surface = volmdlr.surfaces.ToroidalSurface3D(
            volmdlr.Frame3D(tore_center, u, v, axis),
            tore_radius, self.radius)
        return [volmdlr.faces.ToroidalFace3D.from_surface_rectangular_cut(surface, 0, angle, 0, volmdlr.TWO_PI)]

    def point_belongs(self, point: volmdlr.Point3D, abs_tol: float = 1e-6):
        """
        Returns if given point belongs to the Circle3D.
        """
        distance = point.point_distance(self.center)
        vec = volmdlr.Vector3D(*point - self.center)
        dot = self.normal.dot(vec)
        if math.isclose(distance, self.radius, abs_tol=abs_tol) \
                and math.isclose(dot, 0, abs_tol=abs_tol):
            return True
        return False

    def reverse(self):
        """
        Reverses the direction of the circle.

        """
        frame = volmdlr.Frame3D(self.center, self.frame.u, -self.frame.v, self.frame.u.cross(-self.frame.v))
        return Circle3D(frame, self.radius)

    def split(self, split_start, split_end):
        """
        Splits a circle into two arcs, at two given points.

        :param split_start: split point 1.
        :param split_end:  split point 2.
        :return: A list with two split arc 3D.
        """
        return [volmdlr.edges.Arc3D(self, split_start, split_end),
                volmdlr.edges.Arc3D(self, split_end, split_start)]

    def sweep(self, *args):
        """
        Circle 3D is used as path for sweeping given section through it.

        :return:
        """
        _, section_contour = args
        new_faces = []
        for contour_primitive in section_contour.primitives:
            new_faces.extend(contour_primitive.revolution(
                self.center, self.normal, volmdlr.TWO_PI))
        return new_faces

    def distance_linesegment(self, linesegment3d, return_points=False):
        """
        Gets the minimum distance between an Arc 3D and Line Segment 3D.

        :param linesegment3d: other line segment 3d.
        :param return_points: boolean to decide weather to return the corresponding minimal distance points or not.
        :return: minimum distance / minimal distance with corresponding points.
        """
        point1, point2 = vm_common_operations.minimum_distance_points_circle3d_linesegment3d(self, linesegment3d)
        if return_points:
            return point1.point_distance(point2), point1, point2
        return point1.point_distance(point2)

    def get_arc_point_angle(self, point):
        """Returns the angle of point on the circle."""
        local_start_point = self.frame.global_to_local_coordinates(point)
        u1, u2 = local_start_point.x / self.radius, local_start_point.y / self.radius
        point_angle = volmdlr.geometry.sin_cos_angle(u1, u2)
        return point_angle


class Ellipse2D(ClosedCurve):
    """
    Defines an Ellipse in two-dimensions.

    Ellipse2D defined by a major axis (A), minor axis (B), a center and a frame 2d where its u-component
    represents the direction of the major axis.

    :param major_axis: ellipse's major axis (A)
    :type major_axis: float
    :param minor_axis: ellipse's minor axis (B)
    :type minor_axis: float
    :param frame: ellipse's local frame.
    :type frame: volmdlr.Frame2D.

    :Example:
    >>> ellipse2d = Ellipse2D(4, 2, volmdlr.OXY)
    """

    def __init__(self, major_axis, minor_axis, frame, name=''):
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.center = frame.origin
        self.major_dir = frame.u
        self.minor_dir = frame.v
        # self.frame = volmdlr.Frame2D(self.center, self.major_dir, self.minor_dir)
        self.frame = frame
        if math.isclose(frame.u.cross(frame.v), 1.0, abs_tol=1e-6):
            self.angle_start = 0.0
            self.angle_end = volmdlr.TWO_PI
            self.is_trigo = True
        elif math.isclose(frame.u.cross(frame.v), -1.0, abs_tol=1e-6):
            self.angle_start = volmdlr.TWO_PI
            self.angle_end = 0.0
            self.is_trigo = False
        self.theta = geometry.clockwise_angle(self.major_dir, volmdlr.X2D)
        if self.theta == math.pi * 2:
            self.theta = 0.0
        self._bounding_rectangle = None
        ClosedCurve.__init__(self, name=name)

    def __hash__(self):
        return hash((self.center, self.major_dir, self.major_axis, self.minor_axis))

    def __getitem__(self, key):
        if key == 0:
            return self.major_axis
        if key == 1:
            return self.minor_axis
        if key == 2:
            return self.frame
        raise IndexError

    @property
    def bounding_rectangle(self):
        """
        Gets the bounding rectangle of the ellipse 2d.

        :return: a Bounding Rectangle object.
        """
        if not self._bounding_rectangle:
            self._bounding_rectangle = self.get_bounding_rectangle()
        return self._bounding_rectangle

    def get_bounding_rectangle(self):
        """
        Calculates the bounding rectangle of the ellipse 2d.

        :return: a Bounding Rectangle object.
        """
        point1 = self.center - self.major_dir * self.major_axis
        point2 = self.center + self.major_dir * self.major_axis
        point3 = self.center - self.minor_dir * self.minor_axis
        point4 = self.center + self.minor_dir * self.minor_axis
        x_components = [point1.x, point2.x, point3.x, point4.x]
        y_components = [point1.y, point2.y, point3.y, point4.y]
        return volmdlr.core.BoundingRectangle(min(x_components), max(x_components),
                                              min(y_components), max(y_components))

    def area(self):
        """
        Calculates the ellipse's area.

        :return: ellipse's area, float.
        """
        return math.pi * self.major_axis * self.minor_axis

    def length(self):
        """
        Calculates the ellipse's length.

        :return: ellipse's length.
        """
        mid_point = self.center - self.major_axis * self.major_dir
        if self.theta != 0.0:
            mid_point = self.center - volmdlr.Point2D(self.major_axis, 0)
            mid_point = mid_point.rotation(self.center, self.theta)
        length = 2 * self.abscissa(mid_point)
        return length

    def to_3d(self, plane_origin, x, y):
        """
        Transforms a Ellipse2D into an Ellipse3D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Ellipse3D.
        """
        center3d = self.frame.origin.to_3d(plane_origin, x, y)
        major_dir_pointd2d = self.center + self.major_axis * self.major_dir
        major_dir_point = major_dir_pointd2d.to_3d(plane_origin, x, y)
        u_vector = major_dir_point - center3d
        u_vector = u_vector.unit_vector()
        minor_dir_point2d = self.center + self.minor_axis * self.minor_dir
        minor_dir_point = minor_dir_point2d.to_3d(plane_origin, x, y)
        v_vector = minor_dir_point - center3d
        v_vector = v_vector.unit_vector()
        w_vector = u_vector.cross(v_vector)
        frame3d = volmdlr.Frame3D(center3d, u_vector, v_vector, w_vector)
        return Ellipse3D(self.major_axis, self.minor_axis, frame3d)

    def point_over_ellipse(self, point, abs_tol=1e-6):
        """
        Verifies if a point is on the ellipse.

        :param point: point to be verified.
         :param abs_tol: tolerance.
        :return: True or False.
        """
        return math.isclose(
            round(((point.x - self.center.x) * math.cos(self.theta) +
                   (point.y - self.center.y) * math.sin(self.theta)) ** 2 / self.major_axis ** 2 +
                  ((point.x - self.center.x) * math.sin(self.theta) -
                   (point.y - self.center.y) * math.cos(self.theta)) ** 2 / self.minor_axis ** 2, 2), 1.0,
            abs_tol=abs_tol)

    def point_over_contour(self, point, abs_tol=1e-6):
        """
        Verifies if a point is on the ellipse.

        :param point: point to be verified.
        :param abs_tol: tolerance.
        :return: True or False.
        """
        return self.point_over_ellipse(point, abs_tol)

    def line_intersections(self, line: 'Line2D'):
        """
        Calculates the intersections between a line and an ellipse.

        :param line: line to calculate intersections
        :return: list of points intersections, if there are any
        """
        intersections = volmdlr_intersections.ellipse2d_line_intersections(self, line)
        return intersections

    def linesegment_intersections(self, linesegment: 'volmdlr.edges.LineSegment2D', abs_tol: float = 1e-6):
        """
        Calculates the intersections between a line segment and an ellipse.

        :param linesegment: line segment to calculate intersections.
        :param abs_tol: tolerance to be considered while validating an intersection.
        :return: list of points intersections, if there are any.
        """
        line_intersections = self.line_intersections(linesegment.line)
        intersections = []
        for intersection in line_intersections:
            if linesegment.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def ellipse_intersections(self, ellipse2d, abs_tol: float = 1e-7):
        """
        Gets the intersections between two Ellipse 2D.

        :param ellipse2d: The other ellipse.
        :param abs_tol: Tolerance.
        :return:
        """
        if self.bounding_rectangle.distance_to_b_rectangle(ellipse2d.bounding_rectangle) > abs_tol:
            return []
        intersections = volmdlr_intersections.get_bsplinecurve_intersections(ellipse2d, self, abs_tol)
        return intersections

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Calculates the discretized points for the ellipse.

        :param number_points: number of point to have in the discretized points.
        :param angle_resolution: the angle resolution to be used to discretize points.
        :return: discretized points.
        """
        if number_points:
            angle_resolution = number_points
        discretization_points = [self.frame.local_to_global_coordinates(
            volmdlr.Point2D(self.major_axis * math.cos(theta), self.minor_axis * math.sin(theta)))
            for theta in npy.linspace(self.angle_start, self.angle_end, angle_resolution)]
        return discretization_points

    def abscissa(self, point: volmdlr.Point2D, tol: float = 1e-6):
        """
        Calculates the abscissa for a given point.

        :param point: point to calculate the abscissa.
        :param tol: tolerance.
        :return: the corresponding abscissa, 0 < abscissa < ellipse's length.
        """
        if self.point_over_ellipse(point, tol):
            angle_abscissa = self.point_angle_with_major_dir(point)

            def arc_length(theta):
                return math.sqrt((self.major_axis ** 2) * math.sin(theta) ** 2 +
                                 (self.minor_axis ** 2) * math.cos(theta) ** 2)

            res, _ = scipy_integrate.quad(arc_length, 0, angle_abscissa)
            return res
        raise ValueError(f'point {point} does not belong to ellipse')

    def point_at_abscissa(self, abscissa):
        """Get a point at given abscissa."""
        if math.isclose(abscissa, 0.0, abs_tol=1e-6) or math.isclose(abscissa, self.length(), abs_tol=1e-6):
            return self.center + self.major_axis * self.major_dir
        discretized_points = self.discretization_points(number_points=100)
        aproximation_abscissa = 0
        aproximation_point = None
        for point1, point2 in zip(discretized_points[:-1], discretized_points[1:]):
            dist1 = point1.point_distance(point2)
            if aproximation_abscissa + dist1 > abscissa:
                aproximation_point = point1
                break
            aproximation_abscissa += dist1
        initial_point = self.frame.global_to_local_coordinates(aproximation_point)
        u1, u2 = initial_point.x / self.major_axis, initial_point.y / self.minor_axis
        initial_angle = geometry.sin_cos_angle(u1, u2)
        angle_start = 0
        abscissa_angle = vm_common_operations.ellipse_abscissa_angle_integration(
            self, abscissa, angle_start, initial_angle)
        return self.frame.local_to_global_coordinates(
            volmdlr.Point2D(self.major_axis * math.cos(abscissa_angle),
                            self.minor_axis * math.sin(abscissa_angle)))

    def point_distance(self, point):
        """
        Calculates the distance between an Ellipse 2d and point 2d.

        :param point: Other point to calculate distance.
        :type point: volmdlr.Point3D.
        :return: The distance between ellipse and point
        :rtype: float.
        """
        start = self.point_at_abscissa(0.0)
        return vm_common_operations.get_point_distance_to_edge(self, point, start, start)

    def point_angle_with_major_dir(self, point2d):
        """
        Given a point in the ellipse, calculates it angle with the major direction vector.

        """
        initial_point = self.frame.global_to_local_coordinates(point2d)
        u1, u2 = initial_point.x / self.major_axis, initial_point.y / self.minor_axis
        angle_abscissa = geometry.sin_cos_angle(u1, u2)
        return angle_abscissa

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Matplotlib plot for an ellipse.

        """
        if ax is None:
            _, ax = plt.subplots()
        ax = vm_common_operations.plot_from_discretization_points(ax, edge_style, self,
                                                                  number_points=100, close_plot=True)
        if edge_style.equal_aspect:
            ax.set_aspect('equal')
        return ax

    def rotation(self, center, angle: float):
        """
        Rotation of ellipse around a center and an angle.

        :param center: center of the rotation.
        :param angle: angle to rotated of.
        :return: a rotated new ellipse.
        """
        rotated_center = self.center.rotation(center, angle)
        point_major_dir = self.center + self.major_dir * self.major_axis
        rotated_major_dir_point = point_major_dir.rotation(center, angle)
        major_dir = rotated_major_dir_point - rotated_center
        major_dir = major_dir.unit_vector()
        minor_dir = major_dir.normal_vector()
        if not self.is_trigo:
            minor_dir = -minor_dir
        new_frame = volmdlr.Frame2D(rotated_center, major_dir, minor_dir)
        return Ellipse2D(self.major_axis, self.minor_axis, new_frame)

    def translation(self, offset: volmdlr.Vector2D):
        """
        Translation of ellipse from an offset vector.

        :param offset: corresponding translation vector.
        :return: translated new ellipse 2d.
        """
        return Ellipse2D(self.major_axis, self.minor_axis, self.frame.translation(offset))

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new Ellipse2D.

        side = 'old' or 'new'.
        """
        return Ellipse2D(self.major_axis, self.minor_axis, self.frame.frame_mapping(frame, side))

    def reverse(self):
        """
        Reverses the direction of the Ellipse.

        """
        frame = volmdlr.Frame2D(self.center, self.frame.u, -self.frame.v)
        return Ellipse2D(self.major_axis, self.minor_axis, frame)


class Ellipse3D(ClosedCurve):
    """
    Defines a 3D ellipse.

    :param major_axis: Largest radius of the ellipse
    :type major_axis: float
    :param minor_axis: The Smallest radius of the ellipse
    :type minor_axis: float
    :param frame: frame 3d where the ellipse is located.
    """

    def __init__(self, major_axis: float, minor_axis: float,
                 frame, name: str = ''):
        self.frame = frame
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.center = frame.origin
        self.normal = frame.w
        self.major_dir = frame.u
        self.minor_dir = frame.v
        self._self_2d = None
        ClosedCurve.__init__(self, name=name)

    def __getitem__(self, key):
        if key == 0:
            return self.major_axis
        if key == 1:
            return self.minor_axis
        if key == 2:
            return self.frame
        raise IndexError

    @property
    def self_2d(self):
        """
        Version 2d of the ellipse 3d as a property.
        """
        if not self._self_2d:
            self._self_2d = self.to_2d(self.center, self.frame.u, self.frame.v)
        return self._self_2d

    def point_belongs(self, point, tol: float = 1e-6):
        """
        Verifies if a given point lies on the Ellipse3D.

        :param point: point to be verified.
        :param tol: tolerance.
        :return: True is point lies on the Ellipse, False otherwise
        """
        new_point = self.frame.global_to_local_coordinates(point)
        return math.isclose(new_point.x ** 2 / self.major_axis ** 2 +
                            new_point.y ** 2 / self.minor_axis ** 2, 1.0, abs_tol=tol)

    def length(self):
        """
        Calculates the length of the ellipse.

        Ramanujan's approximation for the perimeter of the ellipse.
        P = π (a + b) [ 1 + (3h) / (10 + √(4 - 3h) ) ], where h = (a - b)**2/(a + b)**2
        :return:
        """
        perimeter_formular_h = (self.major_axis - self.minor_axis) ** 2 / (self.major_axis + self.minor_axis) ** 2
        return math.pi * (self.major_axis + self.minor_axis) * \
            (1 + (3 * perimeter_formular_h / (10 + math.sqrt(4 - 3 * perimeter_formular_h))))

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Discretize a Contour to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned.
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc.
        :return: a list of sampled points.
        """
        if number_points:
            angle_resolution = number_points
        discretization_points_3d = [
                                       self.center + self.major_axis * math.cos(
                                           theta) * self.major_dir
                                       + self.minor_axis * math.sin(
                                           theta) * self.minor_dir for theta in
                                       npy.linspace(0, volmdlr.TWO_PI, angle_resolution)]

        return discretization_points_3d

    def to_2d(self, plane_origin, x, y):
        """
        Transforms an Ellipse 3D into an Ellipse 2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Ellipse2D.
        """
        center = self.center.to_2d(plane_origin, x, y)
        major_dir_point3d = self.center + self.major_axis * self.major_dir
        major_dir_point2d = major_dir_point3d.to_2d(plane_origin, x, y)
        major_dir_2d = major_dir_point2d - center
        major_dir_2d = major_dir_2d.unit_vector()
        minor_dir_point3d = self.center + self.minor_axis * self.minor_dir
        minor_dir_point2d = minor_dir_point3d.to_2d(plane_origin, x, y)
        minor_dir_2d = minor_dir_point2d - center
        minor_dir_2d = minor_dir_2d.unit_vector()
        # major_dir_2d = self.major_dir.to_2d()
        # _d2 = self.minor_dir.to_2d(plane_origin, x, y)
        return Ellipse2D(self.major_axis, self.minor_axis, volmdlr.Frame2D(center, major_dir_2d, minor_dir_2d))

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :param tol: tolerance.
        :return: abscissa
        """
        if not self.point_belongs(point, tol):
            raise ValueError('Point is not on ellipse.')
        point2d = point.to_2d(self.center, self.frame.u, self.frame.v)
        return self.self_2d.abscissa(point2d)

    def point_at_abscissa(self, abscissa: float):
        """
        Calculates the 3D point on the curve at a given fraction of its length (abscissa).

        :param abscissa: The fraction of the curve's length at which to calculate the point.
        :type abscissa: (float)
        Returns: The calculated 3D point on the curve.
        :rtype: Point3D.
        """
        point2d = self.self_2d.point_at_abscissa(abscissa)
        return point2d.to_3d(self.center, self.frame.u, self.frame.v)

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D, same_sense: bool = True):
        """
        Trims an ellipse between two points.

        :param point1: point1 used to trim ellipse.
        :param point2: point2 used to trim ellipse.
        :param same_sense: indicates whether the curve direction agrees with (True) or is in the opposite
               direction (False) to the edge direction. By default, it's assumed True
        :return: arc of ellipse between these two points.
        """
        ellipse = self
        if not same_sense:
            ellipse = self.reverse()
        if point1.is_close(point2):
            return volmdlr.edges.FullArcEllipse3D(ellipse, point1, self.name)
        return volmdlr.edges.ArcEllipse3D(ellipse, point1, point2)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Ellipse3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated Ellipse3D.
        """
        return Ellipse3D(self.major_axis, self.minor_axis, self.frame.rotation(center, axis, angle), self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Ellipse 3D translation.

        :param offset: translation vector.
        :return: A new translated Ellipse 3D.
        """
        return Ellipse3D(self.major_axis, self.minor_axis, self.frame.translation(offset), self.name)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Ellipse3D.

        side = 'old' or 'new'.
        """

        return Ellipse3D(self.major_axis, self.minor_axis, self.frame.frame_mapping(frame, side))

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plots an ellipse using Matplotlib."""
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')

        return vm_common_operations.plot_from_discretization_points(ax, edge_style, self, close_plot=True,
                                                                    number_points=100)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Ellipse3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding Ellipse3D object.
        :rtype: :class:`volmdlr.wires.Ellipse3D`
        """
        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        center = object_dict[arguments[1]].origin
        normal = object_dict[arguments[1]].w
        major_dir = object_dict[arguments[1]].u
        major_axis = float(arguments[2]) * length_conversion_factor
        minor_axis = float(arguments[3]) * length_conversion_factor
        return cls(major_axis, minor_axis, volmdlr.Frame3D(center, major_dir, normal.cross(major_dir), normal),
                   arguments[0][1:-1])

    def reverse(self):
        """
        Reverses the direction of the Ellipse.

        """
        frame = volmdlr.Frame3D(self.center, self.frame.u, -self.frame.v,
                                self.frame.u.cross(-self.frame.v))
        return Ellipse3D(self.major_axis, self.minor_axis, frame)

    def line_intersections(self, line):
        """
        Gets intersections between an Ellipse 3D and a Line3D.

        :param line: Other Line 3D.
        :return: A list of points, containing all intersections between the Line 3D and the Ellipse3D.
        """
        return volmdlr_intersections.conic3d_line_intersections(self, line)

    def linesegment_intersections(self, linesegment, abs_tol: float = 1e-6):
        """
        Gets intersections between an Ellipse 3D and a Line3D.

        :param linesegment: Other Line 3D.
        :param abs_tol: tolerance.
        :return: A list of points, containing all intersections between the Line 3D and the Ellipse3D.
        """
        ellipse3d_line_intersections = self.line_intersections(linesegment.line)
        intersections = []
        for intersection in ellipse3d_line_intersections:
            if linesegment.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def ellipse_intersections(self, ellipse, abs_tol: float = 1e-6):
        """
        Gets intersections between an Ellipse 3D and a Line3D.

        :param ellipse: Other Ellipse 3D.
        :param abs_tol: tolerance.
        :return: A list of points, containing all intersections between the two Ellipse3D.
        """
        intersections = []
        if self.frame.w.is_colinear_to(ellipse.frame.w) and\
                math.isclose(self.frame.w.dot(ellipse.frame.origin - self.frame.origin), 0, abs_tol=1e-6):
            ellipse2d = ellipse.to_2d(self.frame.origin, self.frame.u, self.frame.v)
            self_ellipse2d = self.to_2d(self.frame.origin, self.frame.u, self.frame.v)
            intersections_2d = self_ellipse2d.ellipse_intersections(ellipse2d)
            for intersection in intersections_2d:
                intersections.append(intersection.to_3d(self.frame.origin, self.frame.u, self.frame.v))
            return intersections

        plane_intersections = volmdlr_intersections.get_two_planes_intersections(self.frame, ellipse.frame)
        if not plane_intersections:
            return []
        plane_intersections = Line3D(plane_intersections[0], plane_intersections[1])
        self_ellipse3d_line_intersections = volmdlr_intersections.conic3d_line_intersections(self,
                                                                                             plane_intersections)
        ellipse3d_line_intersections = volmdlr_intersections.conic3d_line_intersections(ellipse, plane_intersections)
        for intersection in self_ellipse3d_line_intersections + ellipse3d_line_intersections:
            if volmdlr.core.point_in_list(intersection, intersections):
                continue
            if self.point_belongs(intersection, abs_tol) and ellipse.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections


class HyperbolaMixin(Curve):
    """
    Abstract class for a Hyperbola.
    """
    def __init__(self, frame: Union[volmdlr.Frame2D, volmdlr.Frame3D], semi_major_axis: float,
                 semi_minor_axis: float, name: str = ''):
        self.frame = frame
        self.semi_major_axis = semi_major_axis
        self.semi_minor_axis = semi_minor_axis
        Curve.__init__(self, name=name)

    def __getitem__(self, key):
        if key == 0:
            return self.frame
        if key == 1:
            return self.semi_major_axis
        if key == 2:
            return self.semi_minor_axis
        raise IndexError

    def get_x(self, y):
        """
        For given y component, get the corresponding hyperbola x component, in local coordinates.

        :param y: y component.
        :return: x component.
        """
        x_positive = npy.sqrt(((y ** 2) / (self.semi_minor_axis ** 2) + 1)*(self.semi_major_axis ** 2))
        return x_positive

    def linesegment_intersections(self, linesegment, abs_tol: float = 1e-6):
        """
        Calculates the intersections between a Hyperbola and a Line Segment.

        :param linesegment: the Line Segment.
        :param abs_tol: tolerance.
        :return:a list containing all intersections between the two objects, if any exists.
        """
        line_intersections = self.line_intersections(linesegment.line)
        intersections = []
        for intersection in line_intersections:
            if linesegment.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def trim(self, point1, point2):
        """
        Trims a hyperbola between two points.

        :param point1: point 1 used to trim circle.
        :param point2: point2 used to trim circle.
        """
        _bspline_class = getattr(volmdlr.edges, 'BSplineCurve'+self.__class__.__name__[-2:])
        local_split_start = self.frame.global_to_local_coordinates(point1)
        local_split_end = self.frame.global_to_local_coordinates(point2)
        max_y = max(local_split_start.y, local_split_end.y)
        min_y = min(local_split_start.y, local_split_end.y)
        hyperbola_points = self.get_points(min_y, max_y, 3)
        if not hyperbola_points[0].is_close(point1):
            hyperbola_points = hyperbola_points[::-1]
        point, weight1 = hyperbola_parabola_control_point_and_weight(
            hyperbola_points[0], self.tangent(hyperbola_points[0]),
            hyperbola_points[2], self.tangent(hyperbola_points[2]),
            hyperbola_points[1])
        knotvector = generate_knot_vector(2, 3)
        knot_multiplicity = [1] * len(knotvector)

        bspline = _bspline_class(2, [point1, point, point2], knot_multiplicity, knotvector, [1, weight1, 1])
        return bspline


class Hyperbola2D(HyperbolaMixin):
    """
    Class for Hyperbola 2D.

    :param frame: 2D frame where the hyperbola is at.
    :param semi_major_axis: hyperbola's semi major axis.
    :param semi_minor_axis: hyperbola's semi minor axis.
    """
    def __init__(self, frame: volmdlr.Frame2D, semi_major_axis, semi_minor_axis, name: str = ''):
        self.frame = frame
        self.semi_major_axis = semi_major_axis
        self.semi_minor_axis = semi_minor_axis
        HyperbolaMixin.__init__(self, frame, semi_major_axis, semi_minor_axis, name=name)
  
    def __eq__(self, other):
        if self.frame != other.frame:
            return False
        if not math.isclose(self.semi_major_axis, other.semi_major_axis, abs_tol=1e-6) or\
                not math.isclose(self.semi_minor_axis, other.semi_minor_axis, abs_tol=1e-6):
            return False
        return True

    def get_points(self, min_y: float = None, max_y: float = None, number_points: int = 30):
        """
        Gets hyperbola positive branch points.

        :param number_points: number of point to be generated.
        :param min_y: y-minimal value.
        :param max_y: y-maximal value.
        :return: a List of 2D points for the hyperbola positive branch.
        """
        if not min_y and not max_y:
            min_y, max_y = -self.semi_major_axis * 5, self.semi_major_axis * 5
        y_vals = npy.linspace(min_y, max_y, number_points)
        x_positive_vals = self.get_x(y_vals)
        points_positive_branch = []
        for i, y in enumerate(y_vals):
            points_positive_branch.append(volmdlr.Point2D(x_positive_vals[i], y))
        points_positive_branch = [self.frame.local_to_global_coordinates(point) for point in points_positive_branch]
        return points_positive_branch

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Verifies if point belongs to the Hyperbola.

        :param point: other point.
        :param abs_tol: tolerance.
        :return: True if point belongs, and False otherwise.
        """
        local_point = self.frame.global_to_local_coordinates(point)
        if math.isclose(
                local_point.x ** 2 / self.semi_major_axis ** 2 - local_point.y**2 / self.semi_minor_axis ** 1,
                1, abs_tol=abs_tol):
            return True
        return False

    def get_dx_dy(self, point):
        """
        Gets the dx/dy at a given point of the hyperbola 2d.

        :param point: the other point.
        :return: the dx/dy slope at given point.
        """
        return (self.semi_major_axis ** 2 * point.y) / (self.semi_minor_axis ** 2 * math.sqrt(
            self.semi_major_axis ** 2 * point.y ** 2 / self.semi_minor_axis ** 2 + self.semi_major_axis ** 2))

    def tangent(self, point):
        """
        Calculates the tangent vector to a hyperbola at a given point.

        :param point: The point at which the tangent vector is to be calculated.
        :type point: volmdlr.Point2D.
        :return: The tangent vector to the hyperbola at the given point.
        :rtype: volmdlr.Vector2D.
        """
        # Convert the point to local coordinates within the hyperbola's frame
        point_at_local_coord = self.frame.global_to_local_coordinates(point)

        # Calculate the slope of the tangent line at the given abscissa
        dx_dy = self.get_dx_dy(point)

        # Construct the second point on the tangent line still on hyperbola's frame.
        tangent_second_point = point_at_local_coord + volmdlr.Point2D(dx_dy, 1)

        # Convert the second point back to global coordinates
        global_coord_second_point = self.frame.local_to_global_coordinates(tangent_second_point)

        tangent_vector = global_coord_second_point - point
        tangent_vector = tangent_vector.to_vector()

        return tangent_vector

    def line_intersections(self, line: Line2D):
        """
        Calculates the intersections between a Hyperbola and an infinite Line in 2D.

        :param line: the infinite 2d line.
        :return:a list containing all intersections between the two objects, if any exists.
        """
        line_to_local_coodinates = line.frame_mapping(self.frame, 'new')
        m = line_to_local_coodinates.get_slope()
        c = line_to_local_coodinates.get_y_intersection()
        a_quad_equation = (self.semi_major_axis**2) * (m**2) - self.semi_minor_axis**2
        b_quad_equation = 2*(self.semi_major_axis**2)*m*c
        c_quad_equation = self.semi_major_axis**2 * (self.semi_minor_axis**2 + c**2)
        if c**2 < (self.semi_major_axis**2)*(m**2) - self.semi_minor_axis**2:
            return []
        if a_quad_equation == 0.0:
            return []
        delta = math.sqrt(b_quad_equation**2 - 4*a_quad_equation*c_quad_equation)
        x1 = (-b_quad_equation + delta) / (2*a_quad_equation)
        x2 = (-b_quad_equation - delta) / (2*a_quad_equation)
        y1 = m * x1 + c
        y2 = m * x2 + c
        intersections = []
        if x1 > 0:
            intersections.append(volmdlr.Point2D(x1, y1))
        if x2 > 0:
            intersections.append(volmdlr.Point2D(x2, y2))
        intersections = [self.frame.local_to_global_coordinates(point) for point in intersections]
        return intersections

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new Hyperbola2D.

        side = 'old' or 'new'.
        """
        return Hyperbola2D(self.frame.frame_mapping(frame, side), self.semi_major_axis, self.semi_minor_axis)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Matplotlib plot for a hyperbola in 2D.

        :param ax: Matplotlib 2D axes.
        :param edge_style: the Styles to be applied to the Hyperbola.
        :return: Matplotlib 2D axes.
        """
        if ax is None:
            _, ax = plt.subplots()
        points_positive_branch = self.get_points()
        components_positive_branch = vm_common_operations.plot_components_from_points(points_positive_branch)
        ax.plot(*components_positive_branch, color=edge_style.color, alpha=edge_style.alpha)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True)
        ax.axis('equal')
        return ax


class Hyperbola3D(HyperbolaMixin):
    """
    Class for Hyperbola 3D.

    :param frame: 3D frame where the hyperbola is at.
    :param semi_major_axis: hyperbola's semi major axis.
    :param semi_minor_axis: hyperbola's semi minor axis.
    """
    def __init__(self, frame: volmdlr.Frame3D, semi_major_axis, semi_minor_axis, name: str = ''):
        self._self_2d = None
        HyperbolaMixin.__init__(self, frame, semi_major_axis, semi_minor_axis, name=name)

    @property
    def self_2d(self):
        """Version 2d of the ellipse 3d as a property."""
        if not self._self_2d:
            self._self_2d = self.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        return self._self_2d

    def get_points(self, min_y: float = None, max_y: float = None, number_points: int = 30):
        """
        Gets hyperbola positive branch points.

        :param number_points: number of point to be generated.
        :param min_y: y-minimal value.
        :param max_y: y-maximal value.
        :return: a List of 3D points for the hyperbola positive branch.
        """
        if not min_y and not max_y:
            min_y, max_y = -self.semi_major_axis * 5, self.semi_major_axis * 5
        y_vals = npy.linspace(min_y, max_y, number_points)
        x_positive_vals = self.get_x(y_vals)
        points_positive_branch = []
        for i, y in enumerate(y_vals):
            points_positive_branch.append(volmdlr.Point3D(x_positive_vals[i], y, 0))
        points_positive_branch = [self.frame.local_to_global_coordinates(point) for point in points_positive_branch]
        return points_positive_branch

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Hyperbola 3D into a Hyperbola 2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Hyperbola2D.
        """
        origin = self.frame.origin.to_2d(plane_origin, x, y)
        u_point = self.frame.origin + self.frame.u
        v_point = self.frame.origin + self.frame.v
        u_point2d = u_point.to_2d(plane_origin, x, y)
        v_point2d = v_point.to_2d(plane_origin, x, y)
        u_vector = (u_point2d - origin).to_vector().unit_vector()
        v_vector = (v_point2d - origin).to_vector().unit_vector()
        frame = volmdlr.Frame2D(origin, u_vector, v_vector)
        return Hyperbola2D(frame, self.semi_major_axis, self.semi_minor_axis)

    def point_belongs(self, point, tol: float = 1e-6):
        """
        Verifies if a given point lies on the Hyperbola 3D.

        :param point: point to be verified.
        :param tol: tolerance.
        :return: True is point lies on the Hyperbola 3D, False otherwise
        """
        new_point = self.frame.global_to_local_coordinates(point)
        return math.isclose(new_point.x ** 2 / self.semi_major_axis ** 2 -
                            new_point.y ** 2 / self.semi_minor_axis ** 2, 1.0, abs_tol=tol)

    def tangent(self, point):
        """
        Calculates the tangent vector to a hyperbola at a given point.

        :param point: The point at which the tangent vector is to be calculated.
        :type point: volmdlr.Point3D.
        :return: The tangent vector to the hyperbola at the given point.
        :rtype: volmdlr.Vector3D.
        """
        point_2d = point.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        tangent_2d = self.self_2d.tangent(point_2d)
        point_tangent_2d = point_2d + tangent_2d
        point_tangent_3d = point_tangent_2d.to_3d(self.frame.origin, self.frame.u, self.frame.v)
        return (point_tangent_3d - point).to_vector()

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Hyperbola3D.

        side = 'old' or 'new'.
        """

        return Hyperbola3D(self.frame.frame_mapping(frame, side), self.semi_major_axis, self.semi_minor_axis)

    def line_intersections(self, line):
        """
        Gets intersections between a Hyperbola 3D and a Line 3D.

        :param line: Other Line 3D.
        :return: A list of points, containing all intersections between the Line 3D and the Hyperbola3D.
        """
        return volmdlr_intersections.conic3d_line_intersections(self, line)

    def circle_intersections(self, circle, abs_tol: float = 1e-6):
        """
        Gets intersections between a Hyperbola and Circle 3D.

        :param circle: Other Circle 3D.
        :param abs_tol: tolerance.
        :return: A list of points, containing all intersections between the Line 3D and the Parabola3D.
        """
        return volmdlr_intersections.conic_intersections(self, circle, abs_tol)

    def ellipse_intersections(self, ellipse, abs_tol: float = 1e-6):
        """
        Gets intersections between a Hyperbola and Circle 3D.

        :param ellipse: Other Circle 3D.
        :param abs_tol: tolerance.
        :return: A list of points, containing all intersections between the Line 3D and the Parabola3D.
        """
        return volmdlr_intersections.conic_intersections(self, ellipse, abs_tol)

    def sort_points_along_curve(self, points: List[Union[volmdlr.Point2D, volmdlr.Point3D]]):
        """
        Sort point along a curve.

        :param points: list of points to be sorted.
        :return: sorted points.
        """
        points_ = [self.frame.global_to_local_coordinates(point) for point in points]
        localy_sorted = sorted(points_, key=lambda ip: ip.y)
        sorted_points = [self.frame.local_to_global_coordinates(point) for point in localy_sorted]
        return sorted_points

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Matplotlib plot for a hyperbola in 3D.

        :param ax: Matplotlib 3D axes.
        :param edge_style: the Styles to be applied to the Hyperbola.
        :return: Matplotlib 3D axes.
        """
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')
        points_positive_branch = self.get_points()
        components_positive_branch = vm_common_operations.plot_components_from_points(points_positive_branch)
        ax.plot(*components_positive_branch, color=edge_style.color, alpha=edge_style.alpha)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.grid(True)
        return ax


class ParabolaMixin(Curve):
    """Abstract class for Parabola."""
    def __getitem__(self, key):
        if key == 0:
            return self.frame
        if key == 1:
            return self.focal_length
        raise IndexError

    def _get_y(self, x):
        """
        Evaluate the y-coordinate of the parabola at a given x-coordinate.

        :param x: The x-coordinate of the point.
        :type x: float.
        :return float: The y-coordinate of the point on the parabola.
        """
        return 0.5 * (x ** 2) / (2 * self.focal_length)

    def trim(self, point1, point2):
        """
        Trims a Parabola between two points.

        :param point1: point 1 used to trim circle.
        :param point2: point2 used to trim circle.
        """
        _bspline_class = getattr(volmdlr.edges, 'BezierCurve' + self.__class__.__name__[-2:])
        _line_class = getattr(sys.modules[__name__], 'Line' + self.__class__.__name__[-2:])
        tangent_vector1 = self.tangent(point1)
        tangent_vector2 = self.tangent(point2)
        lineseg1 = _line_class(point1, point1 + tangent_vector1)
        lineseg2 = _line_class(point2, point2 + tangent_vector2)
        line_inters = lineseg1.line_intersections(lineseg2)
        bezier_parabola = _bspline_class(2, [point1, line_inters[0], point2])
        return bezier_parabola


class Parabola2D(ParabolaMixin):
    """
    Class for a Parabola in 2D.

    :param frame: 2D frame where the Parabola will be at. The parabola Opens int the frame.v direction.
    :type frame: Frame2D.
    :param focal_length: the parabola's focal length.
    :type focal_length: float.
    """
    def __init__(self, frame, focal_length: float, name: str = ''):
        self.vertex = frame.origin
        self.focal_length = focal_length
        self.focus = self.vertex + focal_length * frame.v
        self.frame = frame
        self.vrtx_equation_a = 1 / (4 * focal_length)
        ParabolaMixin.__init__(self, name=name)

    def get_points(self, min_x: float = None, max_x: float = None, number_points: int = 30):
        """
        Gets parabola points.

        :param number_points: number of point to be generated.
        :param min_x: x-minimal value.
        :param max_x: x-maximal value.
        :return: a List of 2D points.
        """
        if not min_x and not max_x:
            min_x, max_x = -self.focal_length * 5, self.focal_length * 5
        x_vals = npy.linspace(min_x, max_x, number_points)
        points = []
        for x in x_vals:
            y = self._get_y(x)
            points.append(self.frame.local_to_global_coordinates(volmdlr.Point2D(x, y)))
        return points

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Verifies if point belongs to the Parabola.

        :param point: other point.
        :param abs_tol: tolerance.
        :return: True if point belongs, and False otherwise.
        """
        local_point = self.frame.global_to_local_coordinates(point)
        if math.isclose(local_point.y,
                        self.vrtx_equation_a * local_point.x ** 2, abs_tol=abs_tol):
            return True
        return False

    def line_intersections(self, line: Line2D):
        """
        Gets intersections between a Parabola 2D and a Line 2D.

        :param line: Other Line 2D.
        :return: A list of points, containing all intersections between the Line 2D and the Parabola 2D.
        """
        line_to_local_coodinates = line.frame_mapping(self.frame, 'new')
        m = line_to_local_coodinates.get_slope()
        c = line_to_local_coodinates.get_y_intersection()
        if m**2 > - 4 * self.vrtx_equation_a * c:
            delta = math.sqrt(m**2 - 4 * self.vrtx_equation_a * (-c))
            x1 = (m + delta) / (2 * self.vrtx_equation_a)
            x2 = (m - delta) / (2 * self.vrtx_equation_a)
            y1 = m * x1 + c
            y2 = m * x2 + c
            intersections = [volmdlr.Point2D(x1, y1), volmdlr.Point2D(x2, y2)]
            intersections = [self.frame.local_to_global_coordinates(point) for point in intersections]
            return intersections
        if math.isclose(m**2, - 4 * self.vrtx_equation_a * c, abs_tol=1e-6):
            x = m / (2 * self.vrtx_equation_a)
            return [volmdlr.Point2D(x, m * x + c)]
        return []

    def tangent(self, point):
        """
        Calculates the tangent vector to a parabola at a given point.

        :param point: The point at which the tangent vector is to be calculated.
        :type point: volmdlr.Point2D.
        :return: The tangent vector to the ellipse at the given point.
        :rtype: volmdlr.Vector2D.
        """
        # Convert the point to local coordinates within the parabola's frame
        point_at_local_coord = self.frame.global_to_local_coordinates(point)

        # Calculate the slope of the tangent line at the point
        dy_dx = 2 * self.vrtx_equation_a * point_at_local_coord.x

        # Construct the second point on the tangent line still on parabola's frame.
        tangent_second_point = point_at_local_coord + volmdlr.Point2D(1, dy_dx)

        # Convert the second point back to global coordinates
        global_coord_second_point = self.frame.local_to_global_coordinates(tangent_second_point)

        tangent_vector = global_coord_second_point - point
        tangent_vector = tangent_vector.to_vector()

        return tangent_vector

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Matplotlib plot for a parabola in 2D.

        :param ax: Matplotlib 2D axes.
        :param edge_style: the Styles to be applied to the Parabola.
        :return: Matplotlib 2D axes.
        """
        if ax is None:
            _, ax = plt.subplots()

        points_positive_branch = self.get_points()
        components_positive_branch = vm_common_operations.plot_components_from_points(points_positive_branch)
        ax.plot(*components_positive_branch, color=edge_style.color, alpha=edge_style.alpha)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.grid(True)
        ax.axis('equal')
        return ax


class Parabola3D(ParabolaMixin):
    """
    Class for a Parabola in 3D.

    :param frame: 3D frame where the Parabola will be at. The parabola Opens int the frame.v direction.
    :type frame: Frame3D.
    :param focal_length: the parabola's focal length.
    :type focal_length: float.
    """
    def __init__(self, frame: volmdlr.Frame3D, focal_length: float, name: str = ''):
        self.vertex = frame.origin
        self.focal_length = focal_length
        self.focus = self.vertex + focal_length * frame.v
        self.frame = frame
        self.vrtx_equation_a = 1 / (4 * focal_length)
        ParabolaMixin.__init__(self, name=name)
        self._self_2d = None

    @property
    def self_2d(self):
        """Version 2d of the ellipse 3d as a property."""
        if not self._self_2d:
            self._self_2d = self.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        return self._self_2d

    def get_points(self, min_x: float = None, max_x: float = None, number_points: int = 30):
        """
        Gets parabola 3D points.

        :param number_points: number of point to be generated.
        :param min_x: x-minimal value.
        :param max_x: x-maximal value.
        :return: a List of 3D points.
        """
        if not min_x and not max_x:
            min_x, max_x = -self.focal_length * 5, self.focal_length * 5
        x_vals = npy.linspace(min_x, max_x, number_points)
        points = []
        for x in x_vals:
            y = self._get_y(x)
            points.append(self.frame.local_to_global_coordinates(volmdlr.Point3D(x, y, 0)))
        return points

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Parabola 3D into a Parabola 2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Parabola2D.
        """
        origin = self.frame.origin.to_2d(plane_origin, x, y)
        u_point = self.frame.origin + self.frame.u
        v_point = self.frame.origin + self.frame.v
        u_point2d = u_point.to_2d(plane_origin, x, y)
        v_point2d = v_point.to_2d(plane_origin, x, y)
        u_vector = (u_point2d - origin).to_vector().unit_vector()
        v_vector = (v_point2d - origin).to_vector().unit_vector()
        frame = volmdlr.Frame2D(origin, u_vector, v_vector)
        return Parabola2D(frame, self.focal_length)

    def tangent(self, point):
        """
        Calculates the tangent vector to a parabola at a given point.

        :param point: The point at which the tangent vector is to be calculated.
        :type point: volmdlr.Point3D.
        :return: The tangent vector to the parabola at the given point.
        :rtype: volmdlr.Vector3D.
        """
        point_2d = point.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        tangent_2d = self.self_2d.tangent(point_2d)
        point_tangent_2d = point_2d + tangent_2d
        point_tangent_3d = point_tangent_2d.to_3d(self.frame.origin, self.frame.u, self.frame.v)
        return (point_tangent_3d - point).to_vector()

    def point_belongs(self, point, tol: float = 1e-6):
        """
        Verifies if a given point lies on the Hyperbola 3D.

        :param point: point to be verified.
        :param tol: tolerance.
        :return: True is point lies on the Hyperbola 3D, False otherwise
        """
        new_point = self.frame.global_to_local_coordinates(point)
        return math.isclose(new_point.y, self.vrtx_equation_a * new_point.x**2, abs_tol=tol)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Hyperbola3D.

        side = 'old' or 'new'.
        """

        return Parabola3D(self.frame.frame_mapping(frame, side), self.focal_length)

    def line_intersections(self, line):
        """
        Gets intersections between a Parabola 3D and a Line 3D.

        :param line: Other Line 3D.
        :return: A list of points, containing all intersections between the Line 3D and the Parabola3D.
        """
        return volmdlr_intersections.conic3d_line_intersections(self, line)

    def conic_intersections(self, conic, abs_tol: float = 1e-6):
        """
        Gets intersections between a two conic curves 3D.

        :param conic: Other Line 3D.
        :param abs_tol: tolerance.
        :return: A list of points, containing all intersections between the Line 3D and the Parabola3D.
        """
        if self.frame.w.is_colinear_to(conic.frame.w) and \
                math.isclose(self.frame.w.dot(conic.frame.origin - self.frame.origin), 0, abs_tol=1e-6):
            raise NotImplementedError
        intersections = []
        plane_intersections = volmdlr_intersections.get_two_planes_intersections(self.frame, conic.frame)
        if not plane_intersections:
            return []
        plane_intersections = Line3D(plane_intersections[0], plane_intersections[1])
        self_ellipse3d_line_intersections = volmdlr_intersections.conic3d_line_intersections(self,
                                                                                             plane_intersections)
        ellipse3d_line_intersections = volmdlr_intersections.conic3d_line_intersections(conic, plane_intersections)
        for intersection in self_ellipse3d_line_intersections + ellipse3d_line_intersections:
            if volmdlr.core.point_in_list(intersection, intersections):
                continue
            if self.point_belongs(intersection, abs_tol) and conic.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def sort_points_along_curve(self, points: List[Union[volmdlr.Point2D, volmdlr.Point3D]]):
        """
        Sort point along a curve.

        :param points: list of points to be sorted.
        :return: sorted points.
        """
        points_ = [self.frame.global_to_local_coordinates(point) for point in points]
        localy_sorted = sorted(points_, key=lambda ip: ip.x)
        sorted_points = [self.frame.local_to_global_coordinates(point) for point in localy_sorted]
        return sorted_points

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Matplotlib plot for a parabola in 3D.

        :param ax: Matplotlib 3D axes.
        :param edge_style: the Styles to be applied to the Parabola.
        :return: Matplotlib 3D axes.
        """
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')
        points_positive_branch = self.get_points()
        components_positive_branch = vm_common_operations.plot_components_from_points(points_positive_branch)
        ax.plot(*components_positive_branch, color=edge_style.color, alpha=edge_style.alpha)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.grid(True)
        return ax
