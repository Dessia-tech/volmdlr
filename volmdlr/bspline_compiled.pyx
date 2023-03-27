#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
"""

Cython functions for bspline

"""
from functools import lru_cache
from math import factorial

import cython


@lru_cache(maxsize=10000)
def binomial_coefficient(int k, int i):
    """
    Computes the binomial coefficient (denoted by *k choose i*).

    Please see the following website for details: http://mathworld.wolfram.com/BinomialCoefficient.html

    :param k: size of the set of distinct elements
    :type k: int
    :param i: size of the subsets
    :type i: int
    :return: combination of *k* and *i*
    :rtype: float
    """
    # Special case
    if i > k:
        return float(0)
    # Compute binomial coefficient
    cdef double k_fact = factorial(k)
    cdef double i_fact = factorial(i)
    cdef double k_i_fact = factorial(k - i)
    return k_fact / (k_i_fact * i_fact)


@cython.boundscheck(False)
@cython.wraparound(False)
def find_span_linear(int degree, list knot_vector, int num_ctrlpts, double knot):
    """ Finds the span of a single knot over the knot vector using linear search.

    Alternative implementation for the Algorithm A2.1 from The NURBS Book by Piegl & Tiller.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param num_ctrlpts: number of control points, :math:`n + 1`
    :type num_ctrlpts: int
    :param knot: knot or parameter, :math:`u`
    :type knot: float
    :return: knot span
    :rtype: int
    """
    cdef int span = degree + 1  # Knot span index starts from zero
    while span < num_ctrlpts and knot_vector[span] <= knot:
        span += 1

    return span - 1


def find_spans(int degree, list knot_vector, int num_ctrlpts, list knots, func=find_span_linear):
    """
    Finds spans of a list of knots over the knot vector.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param num_ctrlpts: number of control points, :math:`n + 1`
    :type num_ctrlpts: int
    :param knots: list of knots or parameters
    :type knots: list, tuple
    :param func: function for span finding, e.g. linear or binary search
    :return: list of spans
    :rtype: list
    """
    cdef int i
    cdef double knot
    cdef int num_knots = len(knots)
    cdef list spans = [0.0]*num_knots

    for i in range(num_knots):
        knot = knots[i]
        spans[i] = func(degree, knot_vector, num_ctrlpts, knot)

    return spans


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef basis_function(int degree, list knot_vector, int span, double knot):

    cdef list left = [0.0] * (degree + 1)
    cdef list right = [0.0] * (degree + 1)
    cdef list N = [1.0] * (degree + 1)

    cdef int j, r
    cdef double temp, saved
    for j in range(1, degree + 1):
        left[j] = knot - knot_vector[span + 1 - j]
        right[j] = knot_vector[span + j] - knot
        saved = 0.0
        for r in range(0, j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        N[j] = saved

    return N


def basis_functions(degree, knot_vector, spans, knots):
    """ Computes the non-vanishing basis functions for a list of parameters.

    Wrapper for :func:`.helpers.basis_function` to process multiple span
    and knot values. Uses recurrence to compute the basis functions, also
    known as Cox - de Boor recursion formula.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param spans: list of knot spans
    :type spans:  list, tuple
    :param knots: list of knots or parameters
    :type knots: list, tuple
    :return: basis functions
    :rtype: list
    """
    basis = []
    for span, knot in zip(spans, knots):
        basis.append(basis_function(degree, knot_vector, span, knot))
    return basis


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef basis_function_ders(int degree, list knot_vector, int span, double knot, int order):
    """
    Computes derivatives of the basis functions for a single parameter.

    Implementation of Algorithm A2.3 from The NURBS Book by Piegl & Tiller.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param span: knot span, :math:`i`
    :type span: int
    :param knot: knot or parameter, :math:`u`
    :type knot: float
    :param order: order of the derivative
    :type order: int
    :return: derivatives of the basis functions
    :rtype: list
    """
    # Initialize variables
    cdef int j, k, r
    cdef int s1, s2, j1, j2, pk, rk
    cdef double saved, temp, d
    cdef list left = [1.0 for _ in range(degree + 1)]
    cdef list right = [1.0 for _ in range(degree + 1)]
    cdef list ndu = [[1.0 for _ in range(degree + 1)] for _ in range(degree + 1)]  # N[0][0] = 1.0 by definition

    for j in range(1, degree + 1):
        left[j] = knot - knot_vector[span + 1 - j]
        right[j] = knot_vector[span + j] - knot
        saved = 0.0
        r = 0
        for r in range(r, j):
            # Lower triangle
            ndu[j][r] = right[r + 1] + left[j - r]
            temp = ndu[r][j - 1] / ndu[j][r]
            # Upper triangle
            ndu[r][j] = saved + (right[r + 1] * temp)
            saved = left[j - r] * temp
        ndu[j][j] = saved

    # Load the basis functions
    cdef list ders = [[0.0 for _ in range(degree + 1)] for _ in range((min(degree, order) + 1))]
    for j in range(0, degree + 1):
        ders[0][j] = ndu[j][degree]

    # Start calculating derivatives
    cdef list a = [[1.0 for _ in range(degree + 1)] for _ in range(2)]
    # Loop over function index
    for r in range(0, degree + 1):
        # Alternate rows in array a
        s1 = 0
        s2 = 1
        a[0][0] = 1.0
        # Loop to compute k-th derivative
        for k in range(1, order + 1):
            d = 0.0
            rk = r - k
            pk = degree - k
            if r >= k:
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk]
                d = a[s2][0] * ndu[rk][pk]
            if rk >= -1:
                j1 = 1
            else:
                j1 = -rk
            if (r - 1) <= pk:
                j2 = k - 1
            else:
                j2 = degree - r
            for j in range(j1, j2 + 1):
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j]
                d += (a[s2][j] * ndu[rk + j][pk])
            if r <= pk:
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r]
                d += (a[s2][k] * ndu[r][pk])
            ders[k][r] = d

            # Switch rows
            j = s1
            s1 = s2
            s2 = j

    # Multiply through by the the correct factors
    f = float(degree)
    for k in range(1, order + 1):
        for j in range(0, degree + 1):
            ders[k][j] *= f
        f *= (degree - k)

    # Return the basis function derivatives list
    return ders


def derivatives(dict datadict, tuple parpos, int deriv_order=0):
    """
    Evaluates the n-th order derivatives at the input parametric position.

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :param parpos: parametric position where the derivatives will be computed
    :type parpos: list, tuple
    :param deriv_order: derivative order; to get the i-th derivative
    :type deriv_order: int
    :return: evaluated derivatives
    :rtype: list
    """
    # Geometry data from datadict
    cdef int[2] degree = datadict["degree"]
    cdef tuple knotvector = datadict["knotvector"]
    cdef tuple ctrlpts = datadict["control_points"]
    cdef tuple size = datadict["size"]
    cdef int dimension = datadict["dimension"] + 1 if datadict["rational"] else datadict["dimension"]
    cdef int pdimension = datadict["pdimension"]

    cdef int idx, k, li, s, r, i, dd, cu, cv
    # Algorithm A3.6
    cdef int[2] d = [min(degree[0], deriv_order), min(degree[1], deriv_order)]

    cdef list SKL = [[[0.0 for _ in range(dimension)] for _ in range(deriv_order + 1)] for _ in range(deriv_order + 1)]

    # span = [0 for _ in range(pdimension)]
    cdef int[2] span = [0, 1]
    cdef list basisdrv = [[] for _ in range(pdimension)]
    for idx in range(pdimension):
        span[idx] = find_span_linear(degree[idx], knotvector[idx], size[idx], parpos[idx])
        basisdrv[idx] = basis_function_ders(degree[idx], knotvector[idx], span[idx], parpos[idx], d[idx])
    cdef list t = [0.0] * dimension
    cdef list cp = [0.0] * dimension
    cdef list tmp = [0.0] * dimension
    cdef list temp = [[0.0 for _ in range(dimension)] for _ in range(degree[1] + 1)]
    for k in range(0, d[0] + 1):
        temp = [[0.0 for _ in range(dimension)] for _ in range(degree[1] + 1)]
        for s in range(0, degree[1] + 1):
            tmp = temp[s]
            for r in range(0, degree[0] + 1):
                cu = span[0] - degree[0] + r
                cv = span[1] - degree[1] + s
                cp = ctrlpts[cv + (size[1] * cu)]
                for i in range(dimension):
                    t[i] = tmp[i] + (basisdrv[0][k][r] * cp[i])
                temp[s][:] = t

        dd = min(deriv_order, d[1])
        for li in range(0, dd + 1):
            for s in range(0, degree[1] + 1):
                elem = SKL[k][li]
                tmp = temp[s]
                for i in range(dimension):
                    t[i] = elem[i] + (basisdrv[1][li][s] * tmp[i])
                SKL[k][li][:] = t
    return SKL


def rational_derivatives(dict datadict, tuple parpos, int deriv_order=0):
    """ Evaluates the n-th order derivatives at the input parametric position.

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :param parpos: parametric position where the derivatives will be computed
    :type parpos: list, tuple
    :param deriv_order: derivative order; to get the i-th derivative
    :type deriv_order: int
    :return: evaluated derivatives
    :rtype: list
    """
    cdef int i, j, k, li, ii
    cdef int dimension = datadict["dimension"] + 1 if datadict["rational"] else datadict["dimension"]

    # Call the parent function to evaluate A(u) and w(u) derivatives
    cdef list SKLw = derivatives(datadict, parpos, deriv_order)

    # Generate an empty list of derivatives
    cdef list SKL = [[[0.0 for _ in range(dimension)] for _ in range(deriv_order + 1)] for _ in range(deriv_order + 1)]

    cdef list tmp = [0.0]*(dimension-1)
    cdef list drv = [0.0]*(dimension-1)
    cdef list v = [0.0]*(dimension-1)
    cdef list v2 = [0.0]*(dimension-1)
    cdef list res = [0.0]*(dimension-1)
    # Algorithm A4.4
    for k in range(0, deriv_order + 1):
        for li in range(0, deriv_order + 1):
            # Deep copying might seem a little overkill but we also want to avoid same pointer issues too
            # v = copy.deepcopy(SKLw[k][l])
            v = SKLw[k][li]
            for j in range(1, li + 1):
                drv = SKL[k][li - j]
                for ii in range(dimension - 1):
                    tmp[ii] = v[ii] - (binomial_coefficient(li, j) * SKLw[0][j][-1] * drv[ii])
                v[:] = tmp
            for i in range(1, k + 1):
                drv = SKL[k - i][li]
                for ii in range(dimension - 1):
                    tmp[ii] = v[ii] - (binomial_coefficient(k, i) * SKLw[i][0][-1] * drv[ii])
                v[:] = tmp
                v2 = [0.0 for _ in range(dimension - 1)]
                for j in range(1, li + 1):
                    drv = SKL[k - i][li - j]
                    for ii in range(dimension - 1):
                        tmp[ii] = v2[ii] + (binomial_coefficient(li, j) * SKLw[i][j][-1] * drv[ii])
                    v2[:] = tmp
                for ii in range(dimension - 1):
                    tmp[ii] = v[ii] - (binomial_coefficient(k, i) * v2[ii])
                v[:] = tmp

            for i in range(dimension - 1):
                res[i] = v[i] / SKLw[0][0][-1]

            SKL[k][li][:] = res
    # Return S(u,v) derivatives
    return SKL
