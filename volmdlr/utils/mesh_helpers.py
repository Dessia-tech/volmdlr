"""
Common operations for objects that manipulates meshes.
"""
from numpy.typing import NDArray
import pyfqmr


def perform_decimation(
        vertices: NDArray[float],
        triangles: NDArray[int],
        target_count: int,
        update_rate: int = 5,
        aggressiveness: float = 7.0,
        max_iterations: int = 100,
        verbose: bool = False,
        lossless: bool = False,
        threshold_lossless: float = 1e-3,
        alpha: float = 1e-9,
        k: int = 3,
        preserve_border: bool = True,
):
    """
    Decimate triangular mesh given its vertices and the list of vertice index that forms the triangle of the mesh.

    Vertices of the mesh should be merged (and maybe rounded) for efficient decimation.

    Note: threshold = alpha * pow(iteration + k, aggressiveness)

    :param vertices: An array of 3D vertices specifying the 3D mesh.
    :type vertices: NDArray[float]
    :param triangles: An array of triangles representing the connectivity of the 3D mesh.
    :type triangles: NDArray[int]
    :param target_count: Target number of triangles. Not used if `lossless` is True.
    :type target_count: int
    :param update_rate: Number of iterations between each update. If `lossless` flag is set to True, rate is 1.
    :type update_rate: int
    :param aggressiveness: Parameter controlling the growth rate of the threshold at each iteration when `lossless`
        is False.
    :type aggressiveness: float
    :param max_iterations: Maximal number of iterations.
    :type max_iterations: int
    :param verbose: Control verbosity.
    :type verbose: bool
    :param lossless: Use the lossless simplification method.
    :type lossless: bool
    :param threshold_lossless: Maximal error after which a vertex is not deleted. Only for `lossless` method.
    :type threshold_lossless: float
    :param alpha: Parameter for controlling the threshold growth.
    :type alpha: float
    :param k: Parameter for controlling the threshold growth.
    :type k: int
    :param preserve_border: Flag for preserving vertices on open border.
    :type preserve_border: bool
    """
    # pylint: disable=too-many-arguments
    simplifier = pyfqmr.Simplify()
    simplifier.setMesh(vertices, triangles)
    simplifier.simplify_mesh(
        target_count=target_count,
        update_rate=update_rate,
        aggressiveness=aggressiveness,
        max_iterations=max_iterations,
        verbose=verbose,
        lossless=lossless,
        threshold_lossless=threshold_lossless,
        alpha=alpha,
        K=k,
        preserve_border=preserve_border,
    )

    vertices, faces, _ = simplifier.getMesh()
    return vertices, faces
