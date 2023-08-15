# cython: language_level=3
# distutils: language=c++
"""
This module provides functions for working with Non-Uniform Rational B-Splines (NURBS) geometries.
Some functions are based on the geomdl open-source library (MIT license),
originally created by the geomdl authors: https://github.com/orbingol/NURBS-Python.

We also use as main reference Les Piegl and Wayne Tiller. The NURBS Book. Springer Science & Business Media, 1997

For more information on NURBs and their use in geometric modeling, please refer to the literature
in the field.
The authors of this module make no warranty or guarantee as to the accuracy or suitability of this
code for any particular purpose.
"""

import numpy as np
cimport numpy as cnp
from cython cimport cdivision
from cython cimport boundscheck, wraparound
from libcpp.vector cimport vector
from cpython.mem cimport PyMem_Malloc, PyMem_Free

import volmdlr.nurbs.helpers as helpers


cnp.import_array()


@cdivision(True)
cdef double binomial_coefficient_c(int k, int i):
    """
    Computes the binomial coefficient (denoted by *k choose i*).
    """
    # Special case
    if i > k:
        return 0.0
    if i == 0 or i == k:
        return 1.0
    # Compute binomial coefficient
    cdef int j
    cdef double result = 1.0
    for j in range(i):
        result *= (float(k - j) / float(i - j))
    return result


def find_span_binsearch(int degree, vector[double] knot_vector, int num_ctrlpts, double knot, **kwargs):
    """Finds the span of the knot over the input knot vector using binary search.

    Implementation of Algorithm A2.1 from The NURBS Book by Piegl & Tiller.

    The NURBS Book states that the knot span index always starts from zero, i.e. for a knot vector [0, 0, 1, 1];
    if FindSpan returns 1, then the knot is between the half-open interval [0, 1).

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
    # Get tolerance value
    cdef double tol = kwargs.get("tol", 10e-6)

    # In The NURBS Book; number of knots = m + 1, number of control points = n + 1, p = degree
    # All knot vectors should follow the rule: m = p + n + 1
    cdef int n = num_ctrlpts - 1

    if abs(knot_vector[n + 1] - knot) <= tol:
        return n

    # Set max and min positions of the array to be searched
    cdef int low = degree
    cdef int high = num_ctrlpts

    # The division could return a float value which makes it impossible to use as an array index
    cdef int mid = int((low + high) / 2)
    # Direct int casting would cause numerical errors due to discarding the significand figures (digits after the dot)
    # The round function could return unexpected results, so we add the floating point with some small number
    # This addition would solve the issues caused by the division operation and how Python stores float numbers.
    # E.g. round(13/2) = 6 (expected to see 7)
    mid = int(helpers.round_c(mid + tol))

    # Search for the span
    while (knot < knot_vector[mid]) or (knot >= knot_vector[mid + 1]):
        if knot < knot_vector[mid]:
            high = mid
        else:
            low = mid
        mid = int((low + high) / 2)

    return mid


@boundscheck(False)
@wraparound(False)
cdef int find_span_linear_c(int degree, vector[double] knot_vector, int num_ctrlpts, double knot):
    """ Finds the span of a single knot over the knot vector using linear search."""
    cdef int span = degree + 1  # Knot span index starts from zero
    while span < num_ctrlpts and knot_vector[span] <= knot:
        span += 1

    return span - 1


def find_span_linear(int degree, vector[double] knot_vector, int num_ctrlpts, double knot):
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
    return find_span_linear_c(degree, knot_vector, num_ctrlpts, knot)


@boundscheck(False)
@wraparound(False)
cpdef vector[int] find_spans(int degree, vector[double] knot_vector, int num_ctrlpts, vector[double] knots):
    """Finds spans of a list of knots over the knot vector.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param num_ctrlpts: number of control points, :math:`n + 1`
    :type num_ctrlpts: int
    :param knots: list of knots or parameters
    :type knots: list, tuple
    :param func: function for span finding, e.g. linear or binary search
    :type func: SpanFunc
    :return: list of spans
    :rtype: list
    """
    cdef int i
    cdef vector[int] spans

    for i in range(len(knots)):
        spans.push_back(find_span_linear_c(degree, knot_vector, num_ctrlpts, knots[i]))

    return spans


def find_multiplicity(double knot, vector[double] knot_vector, **kwargs):
    """Finds knot multiplicity over the knot vector.

    Keyword Arguments:
        * ``tol``: tolerance (delta) value for equality checking

    :param knot: knot or parameter, :math:`u`
    :type knot: float
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :return: knot multiplicity, :math:`s`
    :rtype: int
    """
    # Get tolerance value
    cdef double tol = kwargs.get("tol", 10e-8)

    cdef int mult = 0  # initial multiplicity
    cdef int i

    for i in range(len(knot_vector)):
        if abs(knot - knot_vector[i]) <= tol:
            mult += 1

    return mult


@boundscheck(False)
@wraparound(False)
@cdivision(True)
cdef double* basis_function_c(int degree, vector[double] knot_vector, int span, double knot):
    cdef double *left = <double *>PyMem_Malloc((degree + 1) * sizeof(double))
    cdef double *right = <double *>PyMem_Malloc((degree + 1) * sizeof(double))
    cdef double *N = <double *>PyMem_Malloc((degree + 1) * sizeof(double))

    # Initialize N[0] to 1.0 by definition
    N[0] = 1.0

    cdef int j, r
    cdef double temp, saved

    for j in range(1, degree + 1):
        left[j] = knot - knot_vector[span + 1 - j]
        right[j] = knot_vector[span + j] - knot
        saved = 0.0

        for r in range(j):
            temp = N[r] / (right[r + 1] + left[j - r])
            N[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp

        N[j] = saved

    PyMem_Free(left)
    PyMem_Free(right)
    return N


@boundscheck(False)
@wraparound(False)
cpdef vector[double] basis_function(int degree, vector[double] knot_vector, int span, double knot):
    """Computes the non-vanishing basis functions for a single parameter.

    Implementation of Algorithm A2.2 from The NURBS Book by Piegl & Tiller.
    Uses recurrence to compute the basis functions, also known as Cox - de
    Boor recursion formula.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param span: knot span, :math:`i`
    :type span: int
    :param knot: knot or parameter, :math:`u`
    :type knot: float
    :return: basis functions
    :rtype: list
    """
    cdef double *result = basis_function_c(degree, knot_vector, span, knot)
    # Convert the C array to a Python list before returning
    cdef vector[double] result_list = [result[i] for i in range(degree + 1)]
    PyMem_Free(result)  # Free the dynamically allocated memory
    return result_list


cdef double* basis_function_one_c(int degree, vector[double] knot_vector, int span, double knot):
    # Special case at boundaries
    if (
        (span == 0 and knot == knot_vector[0])
        or (span == len(knot_vector) - degree - 2)
        and knot == knot_vector[len(knot_vector) - 1]
    ):
        return [1.0]

    # Knot is outside of span range
    if knot < knot_vector[span] or knot >= knot_vector[span + degree + 1]:
        return [0.0]

    cdef int j, k
    cdef double Uleft, Uright
    cdef double saved, temp
    cdef double *N = <double *>PyMem_Malloc((degree + span + 1) * sizeof(double))

    for j in range(degree + span + 1):
        N[j] = 0.0
    # Initialize the zero th degree basis functions
    for j in range(degree + 1):
        if knot_vector[span + j] <= knot < knot_vector[span + j + 1]:
            N[j] = 1.0

    # Computing triangular table of basis functions
    for k in range(1, degree + 1):
        # Detecting zeros saves computations
        saved = 0.0
        if N[0] != 0.0:
            saved = ((knot - knot_vector[span]) * N[0]) / (knot_vector[span + k] - knot_vector[span])

        for j in range(degree - k + 1):
            Uleft = knot_vector[span + j + 1]
            Uright = knot_vector[span + j + k + 1]

            # Zero detection
            if N[j + 1] == 0.0:
                N[j] = saved
                saved = 0.0
            else:
                temp = N[j + 1] / (Uright - Uleft)
                N[j] = saved + (Uright - knot) * temp
                saved = (knot - Uleft) * temp

    return N


cpdef double basis_function_one(int degree, list knot_vector, int span, double knot):
    """Computes the value of a basis function for a single parameter.

    Implementation of Algorithm 2.4 from The NURBS Book by Piegl & Tiller.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector
    :type knot_vector: list, tuple
    :param span: knot span, :math:`i`
    :type span: int
    :param knot: knot or parameter, :math:`u`
    :type knot: float
    :return: basis function, :math:`N_{i,p}`
    :rtype: float
    """
    cdef double *result_list = basis_function_one_c(degree, np.array(knot_vector, dtype=np.float64), span, knot)
    cdef double result = result_list[0]
    PyMem_Free(result_list)  # Free the dynamically allocated memory
    return result


@boundscheck(False)
@wraparound(False)
cpdef vector[vector[double]] basis_functions(int degree, vector[double] knot_vector, vector[int] spans,
                                             vector[double] knots):
    """Computes the non-vanishing basis functions for a list of parameters.

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
    cdef vector[vector[double]] basis
    cdef size_t i, n = len(spans)
    for i in range(n):
        basis.push_back(basis_function(degree, knot_vector, spans[i], knots[i]))
    return basis


cpdef list basis_function_all(int degree, list knot_vector, int span, double knot):
    """Computes all non-zero basis functions of all degrees from 0 up to the input degree for a single parameter.

    A slightly modified version of Algorithm A2.2 from The NURBS Book by Piegl & Tiller.
    Wrapper for :func:`.helpers.basis_function` to compute multiple basis functions.
    Uses recurrence to compute the basis functions, also known as Cox - de Boor
    recursion formula.

    For instance; if ``degree = 2``, then this function will compute the basis function
    values of degrees **0, 1** and **2** for the ``knot`` value at the input knot ``span``
    of the ``knot_vector``.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector:  knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param span: knot span, :math:`i`
    :type span: int
    :param knot: knot or parameter, :math:`u`
    :type knot: float
    :return: basis functions
    :rtype: list
    """
    cdef list N = [[None for _ in range(degree + 1)] for _ in range(degree + 1)]
    cdef int i, j
    cdef double *b_func
    for i in range(0, degree + 1):
        b_func = basis_function_c(i, np.array(knot_vector, dtype=np.float64), span, knot)
        for j in range(0, i + 1):
            N[j][i] = b_func[j]
        PyMem_Free(b_func)  # free the dynamically allocated memory
    return N


@boundscheck(False)
@wraparound(False)
@cdivision(True)
cpdef vector[vector[double]] basis_function_ders(int degree, vector[double] knot_vector,
                                                 int span, double knot, int order):
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
    cdef int i, j, k, r
    cdef int s1, s2, j1, j2, pk, rk
    cdef double saved, temp, d
    cdef double *left = <double *>PyMem_Malloc((degree + 1) * sizeof(double))
    cdef double *right = <double *>PyMem_Malloc((degree + 1) * sizeof(double))
    for i in range(degree + 1):
        left[i] = 1.0
        right[i] = 1.0

    cdef vector[vector[double]] ndu = vector[vector[double]](degree + 1, vector[double](degree + 1, 1.0))
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
    PyMem_Free(left)
    PyMem_Free(right)
    # Load the basis functions
    cdef vector[vector[double]] ders = vector[vector[double]]((min(degree, order) + 1), vector[double](degree + 1,
                                                                                                       0.0))
    for j in range(0, degree + 1):
        ders[0][j] = ndu[j][degree]

    # Start calculating derivatives
    cdef vector[vector[double]] a = vector[vector[double]](2, vector[double](degree + 1, 1.0))
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

    # Multiply through by the correct factors
    cdef double f = float(degree)
    for k in range(1, order + 1):
        for j in range(0, degree + 1):
            ders[k][j] *= f
        f *= (degree - k)

    # Return the basis function derivatives list
    return ders


cpdef list basis_function_ders_one(int degree, list knot_vector, int span, double knot, int order):
    """Computes the derivative of one basis functions for a single parameter.

    Implementation of Algorithm A2.5 from The NURBS Book by Piegl & Tiller.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot_vector, :math:`U`
    :type knot_vector: list, tuple
    :param span: knot span, :math:`i`
    :type span: int
    :param knot: knot or parameter, :math:`u`
    :type knot: float
    :param order: order of the derivative
    :type order: int
    :return: basis function derivatives
    :rtype: list
    """
    cdef int j, k, jj
    cdef double Uleft, Uright
    cdef double saved, temp
    cdef list ND
    cdef list ders = [0.0 for _ in range(0, order + 1)]

    # Knot is outside of span range
    if (knot < knot_vector[span]) or (knot >= knot_vector[span + degree + 1]):
        for k in range(0, order + 1):
            ders[k] = 0.0

        return ders

    cdef list N = [[0.0 for _ in range(0, degree + 1)] for _ in range(0, degree + 1)]

    # Initializing the zeroth degree basis functions
    for j in range(0, degree + 1):
        if knot_vector[span + j] <= knot < knot_vector[span + j + 1]:
            N[j][0] = 1.0

    # Computing all basis functions values for all degrees inside the span
    for k in range(1, degree + 1):
        saved = 0.0
        # Detecting zeros saves computations
        if N[0][k - 1] != 0.0:
            saved = ((knot - knot_vector[span]) * N[0][k - 1]) / (knot_vector[span + k] - knot_vector[span])

        for j in range(0, degree - k + 1):
            Uleft = knot_vector[span + j + 1]
            Uright = knot_vector[span + j + k + 1]

            # Zero detection
            if N[j + 1][k - 1] == 0.0:
                N[j][k] = saved
                saved = 0.0
            else:
                temp = N[j + 1][k - 1] / (Uright - Uleft)
                N[j][k] = saved + (Uright - knot) * temp
                saved = (knot - Uleft) * temp

    # The basis function value is the zeroth derivative
    ders[0] = N[0][degree]

    # Computing the basis functions derivatives
    for k in range(1, order + 1):
        # Buffer for computing the kth derivative
        ND = [0.0 for _ in range(0, k + 1)]

        # Basis functions values used for the derivative
        for j in range(0, k + 1):
            ND[j] = N[j][degree - k]

        # Computing derivatives used for the kth basis function derivative

        # Derivative order for the k-th basis function derivative
        for jj in range(1, k + 1):
            if ND[0] == 0.0:
                saved = 0.0
            else:
                saved = ND[0] / (knot_vector[span + degree - k + jj] - knot_vector[span])

            # Index of the Basis function derivatives
            for j in range(0, k - jj + 1):
                Uleft = knot_vector[span + j + 1]
                # Wrong in The NURBS Book: -k is missing.
                # The right expression is the same as for saved with the added j offset
                Uright = knot_vector[span + j + degree - k + jj + 1]

                if ND[j + 1] == 0.0:
                    ND[j] = (degree - k + jj) * saved
                    saved = 0.0
                else:
                    temp = ND[j + 1] / (Uright - Uleft)

                    ND[j] = (degree - k + jj) * (saved - temp)
                    saved = temp

        ders[k] = ND[0]

    return ders


cpdef list basis_functions_ders(int degree, list knot_vector, list spans, list knots, int order):
    """Computes derivatives of the basis functions for a list of parameters.

    Wrapper for :func:`.helpers.basis_function_ders` to process multiple span
    and knot values.

    :param degree: degree, :math:`p`
    :type degree: int
    :param knot_vector: knot vector, :math:`U`
    :type knot_vector: list, tuple
    :param spans: list of knot spans
    :type spans:  list, tuple
    :param knots: list of knots or parameters
    :type knots: list, tuple
    :param order: order of the derivative
    :type order: int
    :return: derivatives of the basis functions
    :rtype: list
    """
    cdef int span
    cdef double knot
    cdef list basis_ders = []
    for span, knot in zip(spans, knots):
        basis_ders.append(basis_function_ders(degree, knot_vector, span, knot, order))
    return basis_ders


@boundscheck(False)
@wraparound(False)
def build_coeff_matrix(int degree, vector[double] knotvector, double[:] params, size_t num_points):
    """Builds the coefficient matrix for global interpolation.

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
    cdef double[:, :] matrix_a = np.zeros((num_points, num_points), dtype=np.double)
    cdef int span
    cdef double* basis_func
    cdef size_t i
    cdef size_t j
    for i in range(num_points):
        span = find_span_linear_c(degree, knotvector, num_points, params[i])
        basis_func = basis_function_c(degree, knotvector, span, params[i])
        for j in range(span - degree, span + 1):
            matrix_a[i][j] = basis_func[j - (span - degree)]
    PyMem_Free(basis_func)
    # Return coefficient matrix
    return matrix_a


cpdef list knot_insertion_kv(list[float] knotvector, float u, int span, int r):
    """ Computes the knot vector of the rational/non-rational spline after knot insertion.

    Part of Algorithm A5.1 of The NURBS Book by Piegl & Tiller, 2nd Edition.

    :param knotvector: knot vector
    :param u: knot
    :param span: knot span
    :param r: number of knot insertions
    :return: updated knot vector
    """
    # Initialize variables
    cdef int kv_size = len(knotvector)
    cdef list[float] kv_updated = [0.0 for _ in range(kv_size + r)]

    # Compute new knot vector
    cdef int i
    for i in range(0, span + 1):
        kv_updated[i] = knotvector[i]
    for i in range(1, r + 1):
        kv_updated[span + i] = u
    for i in range(span + 1, kv_size):
        kv_updated[i + r] = knotvector[i]

    # Return the new knot vector
    return kv_updated


def evaluate_curve(dict datadict, double start: float = 0.0, double stop: float = 1.0):
    cdef int degree = datadict["degree"]
    cdef list knotvector = datadict["knotvector"]
    cdef tuple ctrlpts = datadict["control_points"]
    cdef int size = datadict["size"]
    cdef int sample_size = datadict["sample_size"]
    cdef bint rational = datadict["rational"]
    cdef int dimension = datadict["dimension"] + 1 if rational else datadict["dimension"]
    cdef int precision = datadict["precision"]

    if rational:
        return evaluate_curve_rational(degree, knotvector, ctrlpts, size, sample_size,
                                       dimension, precision, start, stop)
    return evaluate_curve_c(degree, knotvector, ctrlpts, size, sample_size, dimension, precision, start, stop)


cdef list evaluate_curve_c(int degree, vector[double] knotvector, vector[vector[double]] ctrlpts, int size,
                           int sample_size, int dimension, int precision, double start, double stop):
    """Evaluates the curve.

    Keyword Arguments:
        * ``start``: starting parametric position for evaluation
        * ``stop``: ending parametric position for evaluation

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :return: evaluated points
    :rtype: list
    """
    # Algorithm A3.1
    cdef list knots = helpers.linspace(start, stop, sample_size, decimals=precision)
    cdef list spans = find_spans(degree, knotvector, size, knots)
    cdef list basis = basis_functions(degree, knotvector, spans, knots)

    cdef list eval_points = []
    cdef int i, idx
    cdef list crvpt
    cdef double crv_p, ctl_p
    for idx in range(len(knots)):
        crvpt = [0.0 for _ in range(dimension)]
        for i in range(0, degree + 1):
            crvpt[:] = [
                crv_p + (basis[idx][i] * ctl_p) for crv_p, ctl_p in zip(crvpt, ctrlpts[spans[idx] - degree + i])
            ]

        eval_points.append(crvpt)

    return eval_points


def derivatives_curve(dict datadict, double parpos, int deriv_order):
    cdef int degree = datadict["degree"]
    cdef list knotvector = datadict["knotvector"]
    cdef tuple ctrlpts = datadict["control_points"]
    cdef int size = datadict["size"]
    cdef bint rational = datadict["rational"]
    cdef int dimension = datadict["dimension"] + 1 if rational else datadict["dimension"]

    if rational:
        return derivatives_curve_rational(degree, knotvector, ctrlpts, size, dimension, parpos, deriv_order)
    return derivatives_curve_c(degree, knotvector, ctrlpts, size, dimension, parpos, deriv_order)


@boundscheck(False)
@wraparound(False)
cdef vector[vector[double]] derivatives_curve_c(int degree, vector[double] knotvector, vector[vector[double]] ctrlpts,
                                                int size, int dimension, double parpos, int deriv_order):
    """Evaluates the n-th order derivatives at the input parametric position.

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :param parpos: parametric position where the derivatives will be computed
    :type parpos: list, tuple
    :param deriv_order: derivative order; to get the i-th derivative
    :type deriv_order: int
    :return: evaluated derivatives
    :rtype: list
    """

    # Algorithm A3.2
    cdef int du = min(degree, deriv_order)

    cdef vector[vector[double]] CK = vector[vector[double]](deriv_order + 1, vector[double](dimension, 0.0))

    cdef int span = find_span_linear_c(degree, knotvector, size, parpos)
    cdef vector[vector[double]] bfunsders = basis_function_ders(degree, knotvector, span, parpos, du)

    cdef size_t k, j, i
    for k in range(0, du + 1):
        for j in range(0, degree + 1):
            for i in range(dimension):
                CK[k][i] = CK[k][i] + (bfunsders[k][j] * ctrlpts[span - degree + j][i])

    # Return the derivatives
    return CK


cdef evaluate_curve_rational(int degree, vector[double] knotvector, vector[vector[double]] ctrlpts, int size,
                             int sample_size, int dimension, int precision, double start, double stop):
    """Evaluates the rational curve.

    Keyword Arguments:
        * ``start``: starting parametric position for evaluation
        * ``stop``: ending parametric position for evaluation

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :return: evaluated points
    :rtype: list
    """
    # Algorithm A4.1
    cdef list crvptw = evaluate_curve_c(degree, knotvector, ctrlpts, size, sample_size,
                                        dimension, precision, start, stop)

    # Divide by weight
    cdef list eval_points = []
    cdef list pt, cpt
    for pt in crvptw:
        cpt = [float(c / pt[-1]) for c in pt[0: (dimension - 1)]]
        eval_points.append(cpt)

    return eval_points


@boundscheck(False)
@wraparound(False)
@cdivision(True)
cdef vector[vector[double]] derivatives_curve_rational(int degree, vector[double] knotvector,
                                                       vector[vector[double]] ctrlpts,
                                                       int size, int dimension, double parpos, int deriv_order):
    """Evaluates the n-th order derivatives at the input parametric position.

    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :param parpos: parametric position where the derivatives will be computed
    :type parpos: list, tuple
    :param deriv_order: derivative order; to get the i-th derivative
    :type deriv_order: int
    :return: evaluated derivatives
    :rtype: list
    """

    # Call the parent function to evaluate A(u) and w(u) derivatives
    cdef vector[vector[double]] CKw = derivatives_curve_c(degree, knotvector, ctrlpts,
                                                          size, dimension, parpos,
                                                          deriv_order)
    # Algorithm A4.2
    cdef vector[vector[double]] CK = vector[vector[double]](deriv_order + 1, vector[double](dimension - 1, 0.0))
    cdef int k, i, j
    cdef vector[double] v
    v.reserve((dimension - 1))
    cdef double binomial_coeff

    for k in range(0, deriv_order + 1):
        for j in range(dimension - 1):
            v[j] = CKw[k][j]
        for i in range(1, k + 1):
            binomial_coeff = binomial_coefficient_c(k, i)
            for j in range(dimension - 1):
                v[j] -= binomial_coeff * CKw[i][dimension - 1] * CK[k - i][j]
        binomial_coeff = 1.0 / CKw[0][dimension - 1]
        for j in range(dimension - 1):
            CK[k][j] = v[j] * binomial_coeff

    # Return C(u) derivatives
    return CK


def evaluate_surface(dict datadict, **kwargs):
    cdef int[2] degree = datadict["degree"]
    cdef list knotvector = datadict["knotvector"]
    cdef tuple ctrlpts = datadict["control_points"]
    cdef int[2] size = datadict["size"]
    cdef int[2] sample_size = datadict["sample_size"]
    cdef bint rational = datadict["rational"]
    cdef int dimension = datadict["dimension"] + 1 if rational else datadict["dimension"]
    cdef int precision = datadict["precision"]

    # Keyword arguments
    cdef double[2] start = kwargs.get("start", [0.0, 0.0])
    cdef double[2] stop = kwargs.get("stop", [1.0, 1.0])

    if rational:
        return evaluate_surface_rational(degree, knotvector, ctrlpts, size, sample_size,
                                       dimension, precision, start, stop)
    return evaluate_surface_c(degree, knotvector, ctrlpts, size, sample_size, dimension, precision, start, stop)


cdef vector[vector[double]] evaluate_surface_c(int[2] degree, vector[double] knotvector, tuple ctrlpts, int[2] size,
                                int[2] sample_size, int dimension, int precision, double[2] start, double[2] stop):
    """
    Evaluates surface.
    """
    # Algorithm A3.5
    cdef vector[double] knots_u = linalg.linspace(start[0], stop[0], sample_size[0], decimals=precision)
    cdef vector[int] spans_u = helpers.find_spans(degree[0], knotvector[0], size[0], knots_u, self._span_func)
    cdef vector[vector[double]] basis_u = helpers.basis_functions(degree[0], knotvector[0], spans_u, knots_u)

    cdef vector[double] knots_v = linalg.linspace(start[1], stop[1], sample_size[1], decimals=precision)
    cdef vector[int] spans_v = helpers.find_spans(degree[1], knotvector[1], size[1], knots_v, self._span_func)
    cdef vector[vector[double]] basis_v = helpers.basis_functions(degree[1], knotvector[1], spans_v, knots_v)

    cdef list eval_points = []
    cdef int i, j, k, m
    cdef int idx_u, idx_v
    cdef list spt, temp
    cdef double pt, tmp
    for i in range(len(spans_u)):
        idx_u = spans_u[i] - degree[0]
        for j in range(len(spans_v)):
            idx_v = spans_v[j] - degree[1]
            spt = [0.0 for _ in range(dimension)]
            for k in range(0, degree[0] + 1):
                temp = [0.0 for _ in range(dimension)]
                for m in range(0, degree[1] + 1):
                    temp[:] = [tmp + (basis_v[j][m] * cp) for tmp, cp in
                               zip(temp, ctrlpts[idx_v + m + (size[1] * (idx_u + k))])]
                spt[:] = [pt + (basis_u[i][k] * tmp) for pt, tmp in zip(spt, temp)]

            eval_points.append(spt)
    return eval_points


cdef vector[vector[double]] evaluate_surface_rational(int[2] degree, vector[double] knotvector, tuple ctrlpts,
                                                      int[2] size, int[2] sample_size, int dimension, int precision,
                                                      double[2] start, double[2] stop):


    """Evaluates the rational surface.
    
    Keyword Arguments:
        * ``start``: starting parametric position for evaluation
        * ``stop``: ending parametric position for evaluation
    
    :param datadict: data dictionary containing the necessary variables
    :type datadict: dict
    :return: evaluated points
    :rtype: list
    """
    # Algorithm A4.3
    cdef vector[vector[double]] cptw = evaluate_surface_c(degree, knotvector, ctrlpts, size, sample_size, dimension,
                                                          precision, start, stop)

    # Divide by weight
    cdef list eval_points = []
    cdef double c
    cdef list pt, cpt
    for pt in cptw:
        cpt = [float(c / pt[-1]) for c in pt[0: (dimension - 1)]]
        eval_points.append(cpt)

    return eval_points
