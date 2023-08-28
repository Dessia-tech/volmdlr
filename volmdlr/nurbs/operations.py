"""
Nurbs main operations algorithms.
"""

from functools import lru_cache
import numpy as np

import volmdlr
from volmdlr.nurbs import core


def knot_insertion(degree, knotvector, ctrlpts, u, **kwargs):
    """
    Computes the control points of the rational/non-rational spline after knot insertion.

    Part of Algorithm A5.1 of The NURBS Book by Piegl & Tiller, 2nd Edition.

    Keyword Arguments:
        * ``num``: number of knot insertions. *Default: 1*
        * ``s``: multiplicity of the knot. *Default: computed via :func:`.find_multiplicity`*
        * ``span``: knot span. *Default: computed via :func:`.find_span_linear`*

    :param degree: degree
    :type degree: int
    :param knotvector: knot vector
    :type knotvector: list, tuple
    :param ctrlpts: control points
    :type ctrlpts: list
    :param u: knot to be inserted
    :type u: float
    :return: updated control points
    :rtype: list

    """
    # Get keyword arguments
    num_insertions = kwargs.get('num', 1)  # number of knot insertions
    knot_multiplicity = kwargs.get('s', core.find_multiplicity(u, knotvector))  # multiplicity
    knot_span = kwargs.get('span', core.find_span_linear(degree, knotvector, len(ctrlpts), u))  # knot span

    # Initialize variables
    num_ctrlpts_orig = len(ctrlpts)
    num_ctrlpts_new = num_ctrlpts_orig + num_insertions

    # Initialize new control points array (control points may be weighted or not)
    ctrlpts_new = [np.empty_like(ctrlpts[0]) for _ in range(num_ctrlpts_new)]

    # Initialize a local array of length p + 1
    temp = [np.empty_like(ctrlpts[0]) for _ in range(degree + 1)]

    # Save unaltered control points
    ctrlpts_new[:knot_span - degree + 1] = ctrlpts[:knot_span - degree + 1]
    ctrlpts_new[knot_span - knot_multiplicity + num_insertions:] = ctrlpts[knot_span - knot_multiplicity:]

    # Start filling the temporary local array which will be used to update control points during knot insertion
    for i in range(degree - knot_multiplicity + 1):
        temp[i] = np.array(ctrlpts[knot_span - degree + i])

    # Insert knot "num_insertions" times
    for j in range(1, num_insertions + 1):
        new_knot_index = knot_span - degree + j
        for i in range(degree - j - knot_multiplicity + 1):
            alpha = knot_insertion_alpha(u, tuple(knotvector), knot_span, i, new_knot_index)
            temp[i] = alpha * temp[i + 1] + (1.0 - alpha) * temp[i]
        ctrlpts_new[new_knot_index] = temp[0]
        ctrlpts_new[knot_span + num_insertions - j - knot_multiplicity] = \
            temp[degree - j - knot_multiplicity]

    # Load remaining control points
    new_knot_index = knot_span - degree + num_insertions
    ctrlpts_new[new_knot_index + 1:knot_span - knot_multiplicity] = \
        temp[1:knot_span - knot_multiplicity - new_knot_index + 1]

    # Return control points after knot insertion
    return ctrlpts_new


@lru_cache(maxsize=128)
def knot_insertion_alpha(u, knotvector, span, idx, leg):
    """
    Computes :math:`\\alpha` coefficient for knot insertion algorithm.

    :param u: knot
    :type u: float
    :param knotvector: knot vector
    :type knotvector: tuple
    :param span: knot span
    :type span: int
    :param idx: index value (degree-dependent)
    :type idx: int
    :param leg: i-th leg of the control points polygon
    :type leg: int
    :return: coefficient value
    :rtype: float
    """
    return (u - knotvector[leg + idx]) / (knotvector[idx + span + 1] - knotvector[leg + idx])


def knot_insertion_kv(knotvector, u, span, num_insertions):
    """
    Computes the knot vector of the rational/non-rational spline after knot insertion.

    Part of Algorithm A5.1 of The NURBS Book by Piegl & Tiller, 2nd Edition.

    :param knotvector: knot vector
    :type knotvector: list, tuple
    :param u: knot
    :type u: float
    :param span: knot span
    :type span: int
    :param num_insertions: number of knot insertions
    :type num_insertions: int
    :return: updated knot vector
    :rtype: list
    """
    # Initialize variables
    kv_size = len(knotvector)
    kv_updated = np.zeros(kv_size + num_insertions)

    # Compute new knot vector
    kv_updated[:span + 1] = knotvector[:span + 1]
    kv_updated[span + 1:span + num_insertions + 1] = u
    kv_updated[span + num_insertions + 1:] = knotvector[span + 1:]

    # Return the new knot vector
    return kv_updated


def insert_knot_curve(obj, param, num, **kwargs):
    """
    Inserts knots n-times to a spline geometry.

    The following code snippet illustrates the usage of this function:

    .. code-block:: python

        # Insert knot u=0.5 to a curve 2 times
        operations.insert_knot(curve, [0.5], [2])

        # Insert knot v=0.25 to a surface 1 time
        operations.insert_knot(surface, [None, 0.25], [0, 1])

        # Insert knots u=0.75, v=0.25 to a surface 2 and 1 times, respectively
        operations.insert_knot(surface, [0.75, 0.25], [2, 1])

        # Insert knot w=0.5 to a volume 1 time
        operations.insert_knot(volume, [None, None, 0.5], [0, 0, 1])

    Please note that input spline geometry object will always be updated if the knot insertion operation is successful.

    Keyword Arguments:
        * ``check_num``: enables/disables operation validity checks. *Default: True*

    :param obj: spline geometry
    :type obj: abstract.SplineGeometry
    :param param: knot(s) to be inserted in [u, v, w] format
    :type param: list, tuple
    :param num: number of knot insertions in [num_u, num_v, num_w] format
    :type num: list, tuple
    :return: updated spline geometry

    """
    # Get keyword arguments
    check_num = kwargs.get('check_num', True)  # can be set to False when the caller checks number of insertions

    if check_num:
        # Check the validity of number of insertions
        if not isinstance(num, (list, tuple)):
            raise ValueError("The number of insertions must be a list or a tuple")

        if len(num) != 1:
            raise ValueError("The length of the num array must be equal to the number of parametric dimensions")

        for val in num:
            if val < 0:
                raise ValueError('Number of insertions must be a positive integer value')

    if param[0] is not None and num[0] > 0:
        # Find knot multiplicity
        knot_multiplicity = core.find_multiplicity(param[0], obj.knotvector)

        # Check if it is possible add that many number of knots
        if check_num and num[0] > obj.degree - knot_multiplicity:
            raise ValueError("Knot " + str(param[0]) + " cannot be inserted " + str(num[0]) + " times")

        # Find knot span
        span = core.find_span_linear(obj.degree, obj.knotvector, len(obj.ctrlpts), param[0])

        # Compute new knot vector
        kv_new = knot_insertion_kv(obj.knotvector, param[0], span, num[0])

        # Compute new control points
        cpts = obj.ctrlptsw if obj.rational else obj.ctrlpts
        cpts_tmp = knot_insertion(obj.degree, obj.knotvector, cpts, param[0],
                                          num=num[0], s=knot_multiplicity, span=span)
        weights = None
        if obj.rational:
            cpts_tmp, weights = separate_ctrlpts_weights(cpts_tmp)

        # Create new curve
        knots = list(sorted(set(kv_new)))
        knot_multiplicities = [core.find_multiplicity(knot, kv_new) for knot in knots]
        point_name = "Point" + obj.__class__.__name__[-2:]
        cpts_tmp = [getattr(volmdlr, point_name)(*point) for point in cpts_tmp]
        obj = obj.__class__(obj.degree, cpts_tmp, knot_multiplicities, knots, weights)
        # obj.control_points = cpts_tmp
        # obj.knotvector = kv_new
    # Return new spline geometry
    return obj

def split_curve(obj, param, **kwargs):
    """
    Splits the curve at the input parametric coordinate.

    This method splits the curve into two pieces at the given parametric coordinate, generates two different
    curve objects and returns them. It does not modify the input curve.

    Keyword Arguments:
        * ``find_span_func``: FindSpan implementation. *Default:* :func:`.helpers.find_span_linear`
        * ``insert_knot_func``: knot insertion algorithm implementation. *Default:* :func:`.operations.insert_knot`

    :param obj: Curve to be split
    :type obj: abstract.Curve
    :param param: parameter
    :type param: float
    :return: a list of curve segments
    :rtype: list

    """
    # # Validate input
    # if not isinstance(obj, abstract.Curve):
    #     raise GeomdlException("Input shape must be an instance of abstract.Curve class")

    if param in set(obj.domain):
        raise ValueError("Cannot split from the domain edge")

    # Find multiplicity of the knot and define how many times we need to add the knot
    knot_span = core.find_span_linear(obj.degree, obj.knotvector, len(obj.ctrlpts), param) - obj.degree + 1
    knot_multiplicity = core.find_multiplicity(param, obj.knotvector)
    insertion_count = obj.degree - knot_multiplicity

    # Insert knot
    temp_obj = insert_knot_curve(obj, [param], num=[insertion_count], check_num=False)

    # Knot vectors
    knot_span_new = core.find_span_linear(temp_obj.degree, temp_obj.knotvector, len(temp_obj.ctrlpts), param) + 1
    curve1_kv = list(temp_obj.knotvector[0:knot_span_new])
    curve1_kv.append(param)
    curve2_kv = list(temp_obj.knotvector[knot_span_new:])
    for _ in range(0, temp_obj.degree + 1):
        curve2_kv.insert(0, param)

    # Control points
    control_points = temp_obj.control_points
    curve1_ctrlpts = control_points[0:knot_span + insertion_count]
    curve2_ctrlpts = control_points[knot_span + insertion_count - 1:]
    if obj.rational:
        curve1_weights = temp_obj.weights[0:knot_span + insertion_count]
        curve2_weights = temp_obj.weights[knot_span + insertion_count - 1:]
    else:
        curve1_weights = None
        curve2_weights = None

    knots_1 = list(sorted(set(curve1_kv)))
    knot_multiplicities_1 = [core.find_multiplicity(knot, curve1_kv) for knot in knots_1]
    # Create a new curve for the first half
    curve1 = temp_obj.__class__(temp_obj.degree, curve1_ctrlpts, knot_multiplicities_1, knots_1, curve1_weights)

    knots_2 = list(sorted(set(curve2_kv)))
    knot_multiplicities_2 = [core.find_multiplicity(knot, curve2_kv) for knot in knots_2]
    # Create another curve for the second half
    curve2 = temp_obj.__class__(temp_obj.degree, curve2_ctrlpts, knot_multiplicities_2, knots_2, curve2_weights)

    # Return the split curves
    ret_val = [curve1, curve2]
    return ret_val


def separate_ctrlpts_weights(ctrlptsw):
    """
    Divides weighted control points by weights to generate unweighted control points and weights vector.

    This function is dimension agnostic, i.e. control points can be in any dimension but the last element of the array
    should indicate the weight.

    :param ctrlptsw: weighted control points
    :type ctrlptsw: list, tuple
    :return: unweighted control points and weights vector
    :rtype: list
    """
    ctrlpts = []
    weights = []
    for ptw in ctrlptsw:
        temp = [float(pw / ptw[-1]) for pw in ptw[:-1]]
        ctrlpts.append(temp)
        weights.append(ptw[-1])

    return [ctrlpts, weights]
