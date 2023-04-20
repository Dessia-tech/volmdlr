"""
This module contains functions for working with Non-Uniform Rational B-Splines (NURBs).
Some of the functions in this module are either partially or totally copied from the geomdl open
source library, which is distributed under the MIT license.
The original source code is credited to the geomdl authors, and their work can be found at
https://github.com/orbingol/NURBS-Python.

Note that this module is also distributed under the MIT license.
Any modifications made to the original geomdl code will be clearly marked in the source code and
will be subject to the same license terms.

For more information on NURBs and their use in geometric modeling, please refer to the literature
in the field.
The authors of this module make no warranty or guarantee as to the accuracy or suitability of this
code for any particular purpose.
"""
from volmdlr.bspline_compiled import find_spans, basis_functions

cdef linspace(double start, double stop, int num, int decimals=18):
    """
    Returns a list of evenly spaced numbers over a specified interval.

    Inspired from Numpy's linspace function:
    https://github.com/numpy/numpy/blob/master/numpy/core/function_base.py

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
    cdef double delta, value
    cdef int div, i
    cdef list result = []

    start = float(start)
    stop = float(stop)
    if abs(start - stop) <= 10e-8:
        return [start]
    num = int(num)
    if num > 1:
        div = num - 1
        delta = stop - start
        for i in range(num):
            value = start + (i * delta / div)
            result.append(float(("{:." + str(decimals) + "f}").format(value)))
        return result
    return [float(("{:." + str(decimals) + "f}").format(start))]


cdef evaluate(dict datadict, span_func, tuple start, tuple stop):
    """
    Evaluates the surface.

    Keyword Arguments:
        * ``start``: starting parametric position for evaluation
        * ``stop``: ending parametric position for evaluation

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :param span_func: Function for searching the right spans for evaluation.
    :type span_func:
    :return: evaluated points
    :rtype: list
    """
    # Geometry data from datadict
    cdef tuple sample_size = datadict["sample_size"]
    cdef int[2] degree = datadict["degree"]
    cdef tuple knotvector = datadict["knotvector"]
    cdef tuple ctrlpts = datadict["control_points"]
    cdef tuple size = datadict["size"]
    cdef int dimension = datadict["dimension"] + 1 if datadict["rational"] else datadict["dimension"]
    cdef int pdimension = datadict["pdimension"]
    cdef int precision = datadict["precision"]

    # Algorithm A3.5
    cdef list spans = [[] for _ in range(pdimension)]
    cdef list basis = [[] for _ in range(pdimension)]
    for idx in range(pdimension):
        knots = linspace(start[idx], stop[idx], sample_size[idx], decimals=precision)
        spans[idx] = find_spans(degree[idx], knotvector[idx], size[idx], knots, span_func)
        basis[idx] = basis_functions(degree[idx], knotvector[idx], spans[idx], knots)

    cdef list eval_points = []
    cdef int i, j, k, m
    cdef int idx_u, idx_v
    cdef list spt, temp
    for i in range(len(spans[0])):
        idx_u = spans[0][i] - degree[0]
        for j in range(len(spans[1])):
            idx_v = spans[1][j] - degree[1]
            spt = [0.0 for _ in range(dimension)]
            for k in range(0, degree[0] + 1):
                temp = [0.0 for _ in range(dimension)]
                for m in range(0, degree[1] + 1):
                    temp[:] = [tmp + (basis[1][j][m] * cp) for tmp, cp in
                               zip(temp, ctrlpts[idx_v + m + (size[1] * (idx_u + k))])]
                spt[:] = [pt + (basis[0][i][k] * tmp) for pt, tmp in zip(spt, temp)]

            eval_points.append(spt)
    return eval_points


cdef evaluate_rational(dict datadict, span_func, tuple start, tuple stop):
    """
    Evaluates the rational surface.

    Keyword Arguments:
        * ``start``: starting parametric position for evaluation
        * ``stop``: ending parametric position for evaluation

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :param span_func: Function for searching the right spans for evaluation.
    :type span_func:
    :return: evaluated points
    :rtype: list
    """
    cdef int dimension = datadict["dimension"] + 1 if datadict["rational"] else datadict["dimension"]
    cdef list cptw = evaluate(datadict, span_func, start=start, stop=stop)
    cdef list eval_points = []
    for pt in cptw:
        cpt = [float(c / pt[-1]) for c in pt[0:(dimension - 1)]]
        eval_points.append(cpt)
    return eval_points


def evaluate_single(tuple param, dict data, bint rational, span_func):
    """ Evaluates the surface at the input (u, v) parameter pair.

    :param param: parameter pair (u, v)
    :type param: list, tuple
    :return: evaluated surface point at the given parameter pair
    :rtype: list
    """
    cdef list pt
    if rational:
        pt = evaluate_rational(datadict=data, span_func=span_func, start=param, stop=param)
        return pt[0]
    # Evaluate the surface point
    pt = evaluate(data, span_func=span_func, start=param, stop=param)
    return pt[0]
