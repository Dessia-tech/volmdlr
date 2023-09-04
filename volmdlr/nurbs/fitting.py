# cython: language_level=3
# distutils: language = c++
# pylint: disable=no-member, used-before-assignment, no-name-in-module, import-error, undefined-variable
"""
Provides curve and surface fitting functions.

This module provides functions for working with Non-Uniform Rational B-Splines (NURBS) geometries.
Some functions are based on the geomdl open-source library (MIT license),
originally created by the geomdl authors: https://github.com/orbingol/NURBS-Python.

We also use as main reference Les Piegl and Wayne Tiller. The NURBS Book. Springer Science & Business Media, 1997

For more information on NURBs and their use in geometric modeling, please refer to the literature
in the field.
The authors of this module make no warranty or guarantee as to the accuracy or suitability of this
code for any particular purpose.
"""

import cython.cimports.libc.math as math_c
import cython
import numpy as np
from scipy.linalg import lu_factor, lu_solve

from cython.cimports.libcpp.vector import vector

from volmdlr.nurbs import core, helpers


@cython.wraparound(False)
@cython.boundscheck(False)
def interpolate_curve(points: np.ndarray[np.double_t, ndim == 2], degree: cython.int,
                      centripetal: cython.bint):
    """
    Curve interpolation through the data points.

    Please refer to Algorithm A9.1 on The NURBS Book (2nd Edition), pp.369-370 for details.

    Keyword Arguments:
        * ``centripetal``: activates centripetal parametrization method. *Default: False*

    :param points: data points
    :type points: list, tuple
    :param degree: degree of the output parametric curve
    :type degree: int
    :return: interpolated B-Spline curve
    :rtype: BSpline.Curve

    """
    # Number of control points
    num_points: cython.Py_ssize_t = len(points)

    u_k: cython.double[:] = compute_params_curve(points, centripetal)

    # Compute knot vector
    knotvector: vector[cython.double] = compute_knot_vector(degree, num_points, u_k)

    # Do global interpolation
    matrix_a: cython.double[:, :] = _build_coeff_matrix(degree, knotvector, u_k, num_points)
    lu_matrix: np.ndarray[np.double_t, ndim == 2]
    piv: np.ndarray[np.int_t, ndim == 1]
    lu_matrix, piv = lu_factor(matrix_a)
    ctrlpts: np.ndarray[np.double_t, ndim == 2] = lu_solve((lu_matrix, piv), points)

    knots = np.unique(knotvector)
    knot_multiplicities = [core.find_multiplicity(knot, knotvector) for knot in knots]

    return ctrlpts, knots, knot_multiplicities


@cython.wraparound(False)
@cython.boundscheck(False)
def interpolate_surface(points, size_u: cython.int, size_v: cython.int, degree_u: cython.int, degree_v: cython.int,
                        **kwargs):
    """
    Surface interpolation through the data points.

    Please refer to the Algorithm A9.4 on The NURBS Book (2nd Edition), pp.380 for details.

    Keyword Arguments:
        * ``centripetal``: activates centripetal parametrization method. *Default: False*

    :param points: data points
    :type points: list, tuple
    :param size_u: number of data points on the u-direction
    :type size_u: int
    :param size_v: number of data points on the v-direction
    :type size_v: int
    :param degree_u: degree of the output surface for the u-direction
    :type degree_u: int
    :param degree_v: degree of the output surface for the v-direction
    :type degree_v: int
    :return: interpolated B-Spline surface
    :rtype: BSpline.Surface

    """
    # Keyword arguments
    use_centripetal: cython.bint = kwargs.get("centripetal", False)
    u_k: np.ndarray[np.double_t, ndim == 1]
    v_l: np.ndarray[np.double_t, ndim == 1]
    u_k, v_l = compute_params_surface(points, size_u, size_v, use_centripetal)

    # Compute knot vectors
    kv_u: vector[cython.double] = compute_knot_vector(degree_u, size_u, u_k)
    kv_v: vector[cython.double] = compute_knot_vector(degree_v, size_v, v_l)

    j: cython.int
    u: cython.int
    v: cython.int
    dim: cython.size_t = len(points[0])
    # Do global interpolation on the u-direction
    ctrlpts_r: cython.double[:, :] = np.zeros((size_u * size_v, dim), dtype=np.float64)
    temp: cython.double[:, :]
    matrix_a: cython.double[:, :]
    for v in range(size_v):
        pts = np.asarray([points[v + (size_v * u)] for u in range(size_u)], dtype=np.float64)
        matrix_a = _build_coeff_matrix(degree_u, kv_u, u_k, pts.shape[0])
        temp = lu_solve(lu_factor(matrix_a), pts)
        for u in range(size_u):
            for j in range(dim):
                ctrlpts_r[u + (size_u * v)][j] = temp[u][j]

    # Do global interpolation on the v-direction
    ctrlpts: np.ndarray[np.double_t, ndim == 2] = np.zeros((size_u * size_v, dim), dtype=np.float64)
    for u in range(size_u):
        pts = np.asarray([ctrlpts_r[u + (size_u * v)] for v in range(size_v)], dtype=np.float64)
        matrix_a = _build_coeff_matrix(degree_v, kv_v, v_l, pts.shape[0])
        temp = lu_solve(lu_factor(matrix_a), pts)
        for v in range(size_v):
            for j in range(dim):
                ctrlpts[v + (size_v * u)][j] = temp[v][j]

    knots_u = np.unique(kv_u)
    knot_multiplicities_u = [core.find_multiplicity(knot, kv_u) for knot in knots_u]
    knots_v = np.unique(kv_v)
    knot_multiplicities_v = [core.find_multiplicity(knot, kv_v) for knot in knots_v]

    return ctrlpts, knots_u, knot_multiplicities_u, knots_v, knot_multiplicities_v


def approximate_curve(points, degree: cython.int, **kwargs):
    """
    Curve approximation using least squares method with fixed number of control points.

    Please refer to The NURBS Book (2nd Edition), pp.410-413 for details.

    Keyword Arguments:
        * ``centripetal``: activates centripetal parametrization method. *Default: False*
        * ``ctrlpts_size``: number of control points. *Default: len(points) - 1*

    :param points: data points
    :type points: list, tuple
    :param degree: degree of the output parametric curve
    :type degree: int
    :return: approximated B-Spline curve
    :rtype: BSpline.Curve

    """
    # Number of data points
    num_dpts: cython.size_t = len(points)  # corresponds to variable "r" in the algorithm

    # Get keyword arguments
    use_centripetal: cython.bint = kwargs.get("centripetal", False)
    num_cpts: cython.int = kwargs.get("ctrlpts_size", num_dpts - 1)

    # Dimension
    dim: cython.size_t = len(points[0])

    u_k: np.ndarray[np.double_t, ndim == 1] = compute_params_curve(points, use_centripetal)

    # Compute knot vector
    knotvector: vector[cython.double] = compute_knot_vector2(degree, num_dpts, num_cpts, u_k)

    # Compute matrix N
    matrix_n: np.ndarray[np.double_t, ndim == 2] = np.zeros((num_dpts - 2, num_cpts - 2), dtype=np.float64)
    for i in range(1, num_dpts - 1):
        for j in range(1, num_cpts - 1):
            matrix_n[i - 1, j - 1] = core.basis_function_one(degree, knotvector, j, u_k[i])

    matrix_nt = matrix_n.T

    matrix_ntn = np.dot(matrix_nt, matrix_n)

    matrix_l, matrix_u = helpers.lu_decomposition(matrix_ntn)
    # Initialize control points array
    ctrlpts = [[0.0 for _ in range(dim)] for _ in range(num_cpts)]

    # Fix start and end points
    ctrlpts[0] = list(points[0])
    ctrlpts[num_cpts - 1] = list(points[num_dpts - 1])

    # Compute - Eq 9.63
    pt0 = points[0]
    ptm = points[num_dpts - 1]
    r_k = []
    for i in range(1, num_dpts - 1):
        ptk = points[i]
        n0p = core.basis_function_one(degree, knotvector, 0, u_k[i])
        nnp = core.basis_function_one(degree, knotvector, num_cpts - 1, u_k[i])
        elem2 = np.array([c * n0p for c in pt0], dtype=np.float64)
        elem3 = np.array([c * nnp for c in ptm], dtype=np.float64)
        r_k.append(ptk - elem2 - elem3)

    # Compute - Eq 9.67
    vector_r = np.zeros((num_cpts - 2, dim), dtype=np.float64)
    for i in range(1, num_cpts - 1):
        ru_tmp = []
        for idx, point in enumerate(r_k):
            ru_tmp.append(point * core.basis_function_one(degree, knotvector, i, u_k[idx + 1]))
        vector_r[i - 1] = np.sum(ru_tmp, axis=0)

    # Compute control points
    for i in range(dim):
        b = [pt[i] for pt in vector_r]
        y = helpers.forward_substitution(matrix_l, b)
        x = helpers.backward_substitution(matrix_u, y)
        for j in range(1, num_cpts - 1):
            ctrlpts[j][i] = x[j - 1]

    knots = np.unique(knotvector)
    knot_multiplicities = [core.find_multiplicity(knot, knotvector) for knot in knots]

    return ctrlpts, knots.tolist(), knot_multiplicities


@cython.wraparound(False)
@cython.boundscheck(False)
def approximate_surface(points, size_u: cython.int, size_v: cython.int, degree_u: cython.int, degree_v: cython.int,
                        **kwargs):
    """
    Surface approximation using least squares method with fixed number of control points.

    This algorithm interpolates the corner control points and approximates the remaining control points.
    Please refer to Algorithm A9.7 of The NURBS Book (2nd Edition), pp.422-423 for details.

    Keyword Arguments:
        * ``centripetal``: activates centripetal parametrization method. *Default: False*
        * ``ctrlpts_size_u``: number of control points on the u-direction. *Default: size_u - 1*
        * ``ctrlpts_size_v``: number of control points on the v-direction. *Default: size_v - 1*

    :param points: data points
    :type points: list, tuple
    :param size_u: number of data points on the u-direction, :math:`r`
    :type size_u: int
    :param size_v: number of data points on the v-direction, :math:`s`
    :type size_v: int
    :param degree_u: degree of the output surface for the u-direction
    :type degree_u: int
    :param degree_v: degree of the output surface for the v-direction
    :type degree_v: int
    :return: approximated B-Spline surface
    :rtype: BSpline.Surface

    """
    # Keyword arguments
    use_centripetal: cython.bint = kwargs.get("centripetal", False)
    num_cpts_u: cython.int = kwargs.get("ctrlpts_size_u", size_u - 1)
    num_cpts_v: cython.int = kwargs.get("ctrlpts_size_v", size_v - 1)

    # Dimension
    dim: cython.size_t = len(points[0])
    u_k: np.ndarray[np.double_t, ndim == 1]
    v_l: np.ndarray[np.double_t, ndim == 1]
    u_k, v_l = compute_params_surface(points, size_u, size_v, use_centripetal)

    # Compute knot vectors
    kv_u: vector[cython.double] = compute_knot_vector2(degree_u, size_u, num_cpts_u, u_k)
    kv_v: vector[cython.double] = compute_knot_vector2(degree_v, size_v, num_cpts_v, v_l)

    matrix_nu: np.ndarray[np.double_t, ndim == 2] = np.zeros((size_u - 2, num_cpts_u - 2), dtype=np.float64)
    i: cython.int
    j: cython.int
    for i in range(1, size_u - 1):
        for j in range(1, num_cpts_u - 1):
            matrix_nu[i - 1, j - 1] = core.basis_function_one(degree_u, kv_u, j, u_k[i])

    matrix_ntu = matrix_nu.T

    matrix_ntnu = np.dot(matrix_ntu, matrix_nu)

    matrix_ntnul, matrix_ntnuu = helpers.lu_decomposition(matrix_ntnu)

    # Fit u-direction
    ctrlpts_tmp = [[0.0 for _ in range(dim)] for _ in range(num_cpts_u * size_v)]
    for j in range(size_v):
        ctrlpts_tmp[j + (size_v * 0)] = list(points[j + (size_v * 0)])
        ctrlpts_tmp[j + (size_v * (num_cpts_u - 1))] = list(points[j + (size_v * (size_u - 1))])
        # Compute Eq 9.63
        pt0 = points[j + (size_v * 0)]
        ptm = points[j + (size_v * (size_u - 1))]
        rku = []
        for i in range(1, size_u - 1):
            ptk = points[j + (size_v * i)]
            n0p = core.basis_function_one(degree_u, kv_u, 0, u_k[i])
            nnp = core.basis_function_one(degree_u, kv_u, num_cpts_u - 1, u_k[i])
            elem2 = [c * n0p for c in pt0]
            elem3 = [c * nnp for c in ptm]
            rku.append([a - b - c for a, b, c in zip(ptk, elem2, elem3)])
        # Compute - Eq 9.67
        r_u = [[0.0 for _ in range(dim)] for _ in range(num_cpts_u - 2)]
        for i in range(1, num_cpts_u - 1):
            ru_tmp = []
            for idx, point in enumerate(rku):
                ru_tmp.append([p * core.basis_function_one(degree_u, kv_u, i, u_k[idx + 1]) for p in point])
            for d in range(dim):
                for idx in range(len(ru_tmp)):
                    r_u[i - 1][d] += ru_tmp[idx][d]
        # Get intermediate control points
        for d in range(dim):
            b = [pt[d] for pt in r_u]
            y = helpers.forward_substitution(matrix_ntnul, b)
            x = helpers.backward_substitution(matrix_ntnuu, y)
            for i in range(1, num_cpts_u - 1):
                ctrlpts_tmp[j + (size_v * i)][d] = x[i - 1]

    matrix_nv: np.ndarray[np.double_t, ndim == 1] = np.zeros((size_v - 2, num_cpts_v - 2), dtype=np.float64)
    for i in range(1, size_v - 1):
        for j in range(1, num_cpts_v - 1):
            matrix_nv[i - 1, j - 1] = core.basis_function_one(degree_v, kv_v, j, v_l[i])

    matrix_ntv = matrix_nv.T

    matrix_ntnv = np.dot(matrix_ntv, matrix_nv)

    matrix_ntnvl, matrix_ntnvu = helpers.lu_decomposition(matrix_ntnv)

    # Fit v-direction
    ctrlpts = [[0.0 for _ in range(dim)] for _ in range(num_cpts_u * num_cpts_v)]
    for i in range(num_cpts_u):
        ctrlpts[0 + (num_cpts_v * i)] = list(ctrlpts_tmp[0 + (size_v * i)])
        ctrlpts[num_cpts_v - 1 + (num_cpts_v * i)] = list(ctrlpts_tmp[size_v - 1 + (size_v * i)])
        # Compute - Eq 9.63
        pt0 = ctrlpts_tmp[0 + (size_v * i)]
        ptm = ctrlpts_tmp[size_v - 1 + (size_v * i)]
        rkv = []
        for j in range(1, size_v - 1):
            ptk = ctrlpts_tmp[j + (size_v * i)]
            n0p = core.basis_function_one(degree_v, kv_v, 0, v_l[j])
            nnp = core.basis_function_one(degree_v, kv_v, num_cpts_v - 1, v_l[j])
            elem2 = [c * n0p for c in pt0]
            elem3 = [c * nnp for c in ptm]
            rkv.append([a - b - c for a, b, c in zip(ptk, elem2, elem3)])
        # Compute - Eq 9.67
        r_v = [[0.0 for _ in range(dim)] for _ in range(num_cpts_v - 2)]
        for j in range(1, num_cpts_v - 1):
            rv_tmp = []
            for idx, point in enumerate(rkv):
                rv_tmp.append([p * core.basis_function_one(degree_v, kv_v, j, v_l[idx + 1]) for p in point])
            for d in range(dim):
                for idx in range(len(rv_tmp)):
                    r_v[j - 1][d] += rv_tmp[idx][d]
        # Get intermediate control points
        for d in range(dim):
            b = [pt[d] for pt in r_v]
            y = helpers.forward_substitution(matrix_ntnvl, b)
            x = helpers.backward_substitution(matrix_ntnvu, y)
            for j in range(1, num_cpts_v - 1):
                ctrlpts[j + (num_cpts_v * i)][d] = x[j - 1]

    knots_u = np.unique(kv_u)
    knot_multiplicities_u = [core.find_multiplicity(knot, kv_u) for knot in knots_u]
    knots_v = np.unique(kv_v)
    knot_multiplicities_v = [core.find_multiplicity(knot, kv_v) for knot in knots_v]

    return ctrlpts, knots_u, knot_multiplicities_u, knots_v, knot_multiplicities_v


@cython.cfunc
@cython.cdivision
@cython.wraparound(False)
@cython.boundscheck(False)
def compute_knot_vector(
    degree: cython.int, num_points: cython.size_t, params: cython.double[:]
) -> vector[cython.double]:
    """
    Computes a knot vector from the parameter list using averaging method.

    Please refer to the Equation 9.8 on The NURBS Book (2nd Edition), pp.365 for details.

    :param degree: degree
    :type degree: int
    :param num_points: number of data points
    :type num_points: int
    :param params: list of parameters, :math:`\\overline{u}_{k}`
    :type params: list, tuple
    :return: knot vector
    :rtype: list
    """
    # Start knot vector
    knotvector: vector[cython.double] = vector[cython.double](degree + 1, 0.0)
    temp_kv: cython.double
    i: cython.size_t
    j: cython.size_t
    # Use averaging method (Eq 9.8) to compute internal knots in the knot vector
    for i in range(num_points - degree - 1):
        temp_kv = (1.0 / float(degree)) * sum(params[j] for j in range(i + 1, i + degree + 1))
        knotvector.push_back(temp_kv)

    # End knot vector
    for _ in range(degree + 1):
        knotvector.push_back(1.0)

    return knotvector


def compute_knot_vector2(degree, num_dpts, num_cpts, params):
    """
    Computes a knot vector ensuring that every knot span has at least one :math:`\\overline{u}_{k}`.

    Please refer to the Equations 9.68 and 9.69 on The NURBS Book (2nd Edition), p.412 for details.

    :param degree: degree
    :type degree: int
    :param num_dpts: number of data points
    :type num_dpts: int
    :param num_cpts: number of control points
    :type num_cpts: int
    :param params: list of parameters, :math:`\\overline{u}_{k}`
    :type params: list, tuple
    :return: knot vector
    :rtype: list
    """
    # Start knot vector
    knotvector = [0.0 for _ in range(degree + 1)]

    # Compute "d" value - Eq 9.68
    d = float(num_dpts) / float(num_cpts - degree)
    # Find internal knots
    for j in range(1, num_cpts - degree):
        i = int(j * d)
        alpha = (j * d) - i
        temp_kv = ((1.0 - alpha) * params[i - 1]) + (alpha * params[i])
        knotvector.append(temp_kv)

    # End knot vector
    knotvector += [1.0 for _ in range(degree + 1)]

    return knotvector


@cython.cfunc
def compute_params_curve(points: np.ndarray[np.double_t, ndim == 2], centripetal: cython.bint = False):
    """
    Computes ū_k for curves.

    Please refer to the Equations 9.4 and 9.5 for chord length parametrization, and Equation 9.6 for centripetal method
    on The NURBS Book (2nd Edition), pp.364-365.

    :param points: data points
    :type points: np.ndarray
    :param centripetal: activates centripetal parametrization method
    :type centripetal: bool
    :return: parameter array, ū_k
    :rtype: np.array
    """
    # Length of the points array
    num_points: cython.Py_ssize_t = points.shape[0]

    # Calculate chord lengths
    cds: np.ndarray[np.double_t, ndim == 1] = np.zeros(num_points + 1, dtype=np.double)
    cds[num_points] = 1.0
    i: cython.Py_ssize_t
    for i in range(1, num_points):
        distance: cython.double = np.linalg.norm(points[i] - points[i - 1])
        cds[i] = math_c.sqrt(distance) if centripetal else distance

    # Find the total chord length
    d: cython.double = np.sum(cds[1:num_points])

    # Divide individual chord lengths by the total chord length
    u_k: np.ndarray[np.double_t, ndim == 1] = np.zeros(num_points, dtype=np.double)
    for i in range(num_points):
        u_k[i] = np.sum(cds[0 : i + 1]) / d

    return u_k


@cython.cfunc
def compute_params_surface(points: np.ndarray[np.double_t, ndim == 2], size_u: cython.int, size_v: cython.int,
                           centripetal: cython.bint = False) -> tuple:
    """
    Computes :math:`\\overline{u}_{k}` and :math:`\\overline{u}_{l}` for surfaces.

    The data points array has a row size of ``size_v`` and column size of ``size_u`` and it is 1-dimensional. Please
    refer to The NURBS Book (2nd Edition), pp.366-367 for details on how to compute :math:`\\overline{u}_{k}` and
    :math:`\\overline{u}_{l}` arrays for global surface interpolation.

    Please note that this function is not a direct implementation of Algorithm A9.3 which can be found on The NURBS Book
    (2nd Edition), pp.377-378. However, the output is the same.

    :param points: data points
    :type points: list, tuple
    :param size_u: number of points on the u-direction
    :type size_u: int
    :param size_v: number of points on the v-direction
    :type size_v: int
    :param centripetal: activates centripetal parametrization method
    :type centripetal: bool
    :return: :math:`\\overline{u}_{k}` and :math:`\\overline{u}_{l}` parameter arrays as a tuple
    :rtype: tuple
    """
    u_k: np.ndarray[np.double_t, ndim == 1] = np.zeros(size_u, dtype=np.double)

    # Compute for each curve on the v-direction
    uk_temp: np.ndarray[np.double_t, ndim == 1] = np.zeros(size_u * size_v, dtype=np.double)
    pts_u: np.ndarray[np.double_t, ndim == 2]
    temp: cython.double[:]
    u: cython.int
    v: cython.int
    for v in range(size_v):
        pts_u = np.asarray([points[v + (size_v * u)] for u in range(size_u)], dtype=np.float64)
        temp = compute_params_curve(pts_u, centripetal)
        for u in range(size_u):
            uk_temp[u + (size_u * v)] = temp[u]

    # Do averaging on the u-direction
    for u in range(size_u):
        knots_v = [uk_temp[u + (size_u * v)] for v in range(size_v)]
        u_k[u] = sum(knots_v) / size_v

    v_l: np.ndarray[np.double_t, ndim == 1] = np.zeros(size_v, dtype=np.double)

    # Compute for each curve on the u-direction
    vl_temp: np.ndarray[np.double_t, ndim == 1] = np.zeros(size_u * size_v, dtype=np.double)
    pts_u: np.ndarray[np.double_t, ndim == 2]
    for u in range(size_u):
        pts_v = np.asarray([points[v + (size_v * u)] for v in range(size_v)], dtype=np.float64)
        temp = compute_params_curve(pts_v, centripetal)
        for v in range(size_v):
            vl_temp[v + (size_v * u)] = temp[v]

    # Do averaging on the v-direction
    for v in range(size_v):
        knots_u = [vl_temp[v + (size_v * u)] for u in range(size_u)]
        v_l[v] = sum(knots_u) / size_u

    return u_k, v_l


@cython.cfunc
def _build_coeff_matrix(
    degree: cython.int, knotvector: vector[cython.double], params: cython.double[:], num_points: cython.size_t
) -> cython.double[:, :]:
    """
    Builds the coefficient matrix for global interpolation.

    This function only uses data points to build the coefficient matrix. Please refer to The NURBS Book (2nd Edition),
    pp364-370 for details.

    :param degree: degree
    :type degree: int
    :param knotvector: knot vector
    :type knotvector: list, tuple
    :param params: list of parameters
    :type params: list, tuple
    :param points: data points
    :type points: list, tuple
    :return: coefficient matrix
    :rtype: list
    """

    # Set up coefficient matrix
    matrix_a: np.ndarray[np.double_t, ndim == 2] = core.build_coeff_matrix(degree, knotvector, params, num_points)

    return matrix_a
