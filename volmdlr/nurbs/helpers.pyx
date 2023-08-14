# cython: language_level=3
"""
Helpers.
"""
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