# cython: language_level=3
"""
Helpers.
"""
from copy import deepcopy
from functools import lru_cache
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from cython cimport cdivision
import cython.cimports.libc.math as math_c


@cdivision(True)
cpdef double round_c(double num, int digits=0):
    cdef double multiplier = math_c.pow(10.0, digits)
    return float(math_c.round(num * multiplier)) / multiplier


def linspace(double start, double stop, int num, int decimals=18):
    """Returns a list of evenly spaced numbers over a specified interval.

    Inspired from Numpy's linspace function: https://github.com/numpy/numpy/blob/master/numpy/core/function_base.py

    :param start: starting value
    :type start: float
    :param stop: end value
    :type stop: float
    :param num: number of samples to generate
    :type num: int
    :param decimals: number of significands
    :type decimals: int
    :return: a list of equally spaced numbers
    :rtype: list
    """
    cdef double delta
    cdef int div, x
    cdef double step = 0.0
    cdef double *result = <double *>PyMem_Malloc(num * sizeof(double))

    if abs(start - stop) <= 10e-8:
        return [start]

    div = num - 1
    delta = stop - start
    try:
        for x in range(num):
            step = start + (x * delta / div)
            result[x] = round_c(step, decimals)
        return [result[i] for i in range(num)]
    finally:
        PyMem_Free(result)  # Free the dynamically allocated memory


def generate_knot_vector(degree, num_ctrlpts, **kwargs):
    """Generates an equally spaced knot vector.

    It uses the following equality to generate knot vector: :math:`m = n + p + 1`

    where;

    * :math:`p`, degree
    * :math:`n + 1`, number of control points
    * :math:`m + 1`, number of knots

    Keyword Arguments:

        * ``clamped``: Flag to choose from clamped or unclamped knot vector options. *Default: True*

    :param degree: degree
    :type degree: int
    :param num_ctrlpts: number of control points
    :type num_ctrlpts: int
    :return: knot vector
    :rtype: list
    """
    if degree == 0 or num_ctrlpts == 0:
        raise ValueError("Input values should be different than zero.")

    # Get keyword arguments
    clamped = kwargs.get("clamped", True)

    # Number of repetitions at the start and end of the array
    num_repeat = degree

    # Number of knots in the middle
    num_segments = num_ctrlpts - (degree + 1)

    if not clamped:
        # No repetitions at the start and end
        num_repeat = 0
        # Should conform the rule: m = n + p + 1
        num_segments = degree + num_ctrlpts - 1

    # First knots
    knot_vector = [0.0 for _ in range(0, num_repeat)]

    # Middle knots
    knot_vector += linspace(0.0, 1.0, num_segments + 2)

    # Last knots
    knot_vector += [1.0 for _ in range(0, num_repeat)]

    # Return auto-generated knot vector
    return knot_vector


def standardize_knot_vector(knot_vector):
    """
    Standardize a knot vector to range from 0 to 1.
    """
    decimals = 18
    first_knot = float(knot_vector[0])
    last_knot = float(knot_vector[-1])
    denominator = last_knot - first_knot

    knot_vector_out = [
        float(("{:." + str(decimals) + "f}").format((float(kv) - first_knot) / denominator)) for kv in knot_vector
    ]

    return knot_vector_out


def insert_knots_and_mutiplicity(knots, knot_mutiplicities, knot_to_add, num):
    """
    Compute knot-elements and multiplicities based on the global knot vector.

    """
    new_knots = []
    new_knot_mutiplicities = []
    i = 0
    for i, knot in enumerate(knots):
        if knot > knot_to_add:
            new_knots.extend([knot_to_add])
            new_knot_mutiplicities.append(num)
            new_knots.extend(knots[i:])
            new_knot_mutiplicities.extend(knot_mutiplicities[i:])
            break
        new_knots.append(knot)
        new_knot_mutiplicities.append(knot_mutiplicities[i])
    return new_knots, new_knot_mutiplicities, i


def doolittle(matrix_a):
    """Doolittle's Method for LU-factorization.

    :param matrix_a: Input matrix (must be a square matrix)
    :type matrix_a: list, tuple
    :return: a tuple containing matrices (L,U)
    :rtype: tuple
    """
    # Initialize L and U matrices
    matrix_u = [[0.0 for _ in range(len(matrix_a))] for _ in range(len(matrix_a))]
    matrix_l = [[0.0 for _ in range(len(matrix_a))] for _ in range(len(matrix_a))]

    # Doolittle Method
    for i in range(0, len(matrix_a)):
        for k in range(i, len(matrix_a)):
            # Upper triangular (U) matrix
            matrix_u[i][k] = float(matrix_a[i][k] - sum([matrix_l[i][j] * matrix_u[j][k] for j in range(0, i)]))
            # Lower triangular (L) matrix
            if i == k:
                matrix_l[i][i] = 1.0
            else:
                matrix_l[k][i] = float(matrix_a[k][i] - sum([matrix_l[k][j] * matrix_u[j][i] for j in range(0, i)]))
                # Handle zero division error
                try:
                    matrix_l[k][i] /= float(matrix_u[i][i])
                except ZeroDivisionError:
                    matrix_l[k][i] = 0.0

    return matrix_l, matrix_u

def lu_decomposition(matrix_a):
    """LU-Factorization method using Doolittle's Method for solution of linear systems.

    Decomposes the matrix :math:`A` such that :math:`A = LU`.

    The input matrix is represented by a list or a tuple. The input matrix is **2-dimensional**, i.e. list of lists of
    integers and/or floats.

    :param matrix_a: Input matrix (must be a square matrix)
    :type matrix_a: list, tuple
    :return: a tuple containing matrices L and U
    :rtype: tuple
    """
    # Check if the 2-dimensional input matrix is a square matrix
    q = len(matrix_a)
    for idx, m_a in enumerate(matrix_a):
        if len(m_a) != q:
            raise ValueError(
                "The input must be a square matrix. " + "Row " + str(idx + 1) + " has a size of " + str(len(m_a)) + "."
            )

    # Return L and U matrices
    return doolittle(matrix_a)


def forward_substitution(matrix_l, matrix_b):
    """Forward substitution method for the solution of linear systems.

    Solves the equation :math:`Ly = b` using forward substitution method
    where :math:`L` is a lower triangular matrix and :math:`b` is a column matrix.

    :param matrix_l: L, lower triangular matrix
    :type matrix_l: list, tuple
    :param matrix_b: b, column matrix
    :type matrix_b: list, tuple
    :return: y, column matrix
    :rtype: list
    """
    q = len(matrix_b)
    matrix_y = [0.0 for _ in range(q)]
    matrix_y[0] = float(matrix_b[0]) / float(matrix_l[0][0])
    for i in range(1, q):
        matrix_y[i] = float(matrix_b[i]) - sum([matrix_l[i][j] * matrix_y[j] for j in range(0, i)])
        matrix_y[i] /= float(matrix_l[i][i])
    return matrix_y


def backward_substitution(matrix_u, matrix_y):
    """Backward substitution method for the solution of linear systems.

    Solves the equation :math:`Ux = y` using backward substitution method
    where :math:`U` is a upper triangular matrix and :math:`y` is a column matrix.

    :param matrix_u: U, upper triangular matrix
    :type matrix_u: list, tuple
    :param matrix_y: y, column matrix
    :type matrix_y: list, tuple
    :return: x, column matrix
    :rtype: list
    """
    q = len(matrix_y)
    matrix_x = [0.0 for _ in range(q)]
    matrix_x[q - 1] = float(matrix_y[q - 1]) / float(matrix_u[q - 1][q - 1])
    for i in range(q - 2, -1, -1):
        matrix_x[i] = float(matrix_y[i]) - sum([matrix_u[i][j] * matrix_x[j] for j in range(i, q)])
        matrix_x[i] /= float(matrix_u[i][i])
    return matrix_x
