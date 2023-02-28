from math import sqrt
import numpy as np
from geomdl import BSpline
from scipy.linalg import lu_factor, lu_solve
from volmdlr.bspline_compiled import find_span_linear, basis_function


def interpolate_curve(list points, int degree, bint centripetal: bool = False):
    """ Curve interpolation through the data points.

    Please refer to Algorithm A9.1 on The NURBS Book (2nd Edition), pp.369-370 for details.

    Keyword Arguments:
        * ``centripetal``:

    :param points: data points
    :type points: list, tuple
    :param degree: degree of the output parametric curve
    :type degree: int
    :param centripetal: activates centripetal parametrization method. *Default: False*
    :type centripetal: bool
    :return: interpolated B-Spline curve
    :rtype: BSpline.Curve
    """
    # Number of control points
    cdef int num_points = len(points)
    cdef list uk, kv, matrix_a
    # Get uk
    uk = compute_params_curve(points, centripetal)

    # Compute knot vector
    kv = compute_knot_vector(degree, num_points, uk)

    matrix_a = _build_coeff_matrix(degree, kv, uk, points)
    # Compute the LU decomposition of the coefficient matrix A
    LU, piv = lu_factor(np.asarray(matrix_a))
    x = lu_solve((LU, piv), np.asarray(points))
    ctrlpts = x.tolist()
    # Generate B-spline curve
    curve = BSpline.Curve()
    curve.degree = degree
    curve.ctrlpts = ctrlpts
    curve.knotvector = kv

    return curve


cdef compute_params_curve(list points, bint centripetal: bool = False):
    """
    Computes :math:`\\overline{u}_{k}` for curves.

    Please refer to the Equations 9.4 and 9.5 for chord length parametrization, and Equation 9.6 for centripetal method
    on The NURBS Book (2nd Edition), pp.364-365.

    :param points: data points
    :type points: list, tuple
    :param centripetal: activates centripetal parametrization method
    :type centripetal: bool
    :return: parameter array, :math:`\\overline{u}_{k}`
    :rtype: list
    """
    # Length of the points array
    cdef int num_points = len(points)

    # Allocate memory for chord lengths
    cdef list cds = [0.0] * (num_points + 1)

    # Calculate chord lengths
    cds[-1] = 1.0
    cdef int i
    cdef double distance
    for i in range(1, num_points):
        p1 = np.array(points[i])
        p2 = np.array(points[i - 1])
        distance = np.linalg.norm(p1 - p2)
        cds[i] = sqrt(distance) if centripetal else distance

    # Find the total chord length
    cdef double d = np.sum(cds[1:-1])

    # Allocate memory for parameter array
    cdef list uk = [0] * (num_points)

    # Divide individual chord lengths by the total chord length
    cdef double s = 0.0
    for i in range(num_points):
        s += cds[i]
        uk[i] = s / d

    return uk


cdef compute_knot_vector(int degree, int num_points, list params):
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
    cdef list kv = [0.0] * (degree + 1)

    # Use averaging method (Eqn 9.8) to compute internal knots in the knot vector
    cdef int i, j
    cdef double temp_kv
    for i in range(num_points - degree - 1):
        temp_kv = 0.0
        for j in range(i + 1, i + degree + 1):
            temp_kv += params[j]
        kv.append((1.0 / degree) * temp_kv)

    # End knot vector
    kv += [1.0] * (degree + 1)

    return kv


cdef _build_coeff_matrix(int degree, list knotvector, list params, list points):
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
    # Number of data points
    cdef int num_points = len(points)

    # Set up coefficient matrix
    cdef list matrix_a = [[0.0 for _ in range(num_points)] for _ in range(num_points)]
    cdef int i, span
    for i in range(num_points):
        span = find_span_linear(degree, knotvector, num_points, params[i])
        matrix_a[i][span-degree:span+1] = basis_function(degree, knotvector, span, params[i])

    # Return coefficient matrix
    return matrix_a
