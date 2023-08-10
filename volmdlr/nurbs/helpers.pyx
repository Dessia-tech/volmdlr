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