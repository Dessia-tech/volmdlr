import os
import time
import re

import matplotlib.pyplot as _plt

scripts = [
    # Core.py
    "core/points.py",
    "core/frames.py",
    "core/grid2d_with_direction.py",
    "core/points_axial_symmetry.py",
    # geometry
    "geometry.py",
    # Edges
    "edges/arcs2D.py",
    "edges/arcs3D.py",
    "edges/bspline.py",
    "edges/bspline2.py",
    "edges/areas_moment_cog_check.py",
    "edges/lines2d.py",
    # Wires
    "wires/roundedlines.py",
    "wires/polygon2D.py",
    "wires/polygon3d.py",
    "wires/triangle2D.py",
    "wires/axial_symmetry.py",
    # Primitives
    "primitives/extrusion.py",
    "primitives/sweep.py",
    "primitives/bspline_sweep.py",
    "primitives/revolved_profile.py",
    "primitives/block3d.py",
    "primitives/sphere.py",
    "primitives/cone.py",
    "primitives/cylinder.py",
    "primitives/hollow_cylinder.py",
    # Faces
    "faces/triangle3d.py",
    "faces/bspline.py",
    "faces/bspline_bark.py",
    "faces/union.py",
    "faces/blocks_set_operations.py",
    "faces/surface2d.py",
    # 'faces/export_to_step.py', #TO BE USED WHEN THE EXPORT IS IMPROVED
    "faces/BSplineSurface/bspline_curves.py",
    "faces/BSplineSurface/bspline_curves_point_belongs.py",
    "faces/BSplineSurface/bspline_surface_interpolation.py",
    "faces/BSplineSurface/bspline_surface_merge.py",
    "faces/BSplineSurface/bspline_surface_split.py",
    "faces/BSplineSurface/from_cylindrical_to_bspline_surface.py",
    "faces/BSplineSurface/bspline_surface_definition.py",
    "faces/BSplineSurface/bspline_surfaces_grid3d.py",
    # 'faces/faces_with_inner_contours.py', #TO BE USED WHEN HOLES IS MERGED
    # Shells
    "shells/operations.py",
    # Cloud
    "cloud/sewing_two_polygons.py",
    "cloud/sewing_stl.py",
    # Contours
    "contours/contour_split.py",
    "contours/contour_merge_with_1.py",
    "contours/contour_merge_with_2.py",
    "contours/cut_by_line.py",
    "contours/contour_cut_by_wire.py",
    # Steps
    "step/read_steps.py",
    # Stl
    "stl/read_stls.py",
    # Distance
    "distance/arc3D_arc3D.py",
    "distance/arc3D_ls3D.py",
    # Showcases
    "showcases/simple_shapes.py",
    "showcases/casing.py",
    "showcases/vessel.py",
    # Mesh
    "mesh/read_msh_file.py",
    # 'mesh/geo_file_1.py',
    # 'mesh/geo_file_2.py',
    # 'mesh/geo_file_3.py',
    # display.py
    "display/mesh3d.py",
    "display/mesh_decimation.py",
    # cad_simplification
    "cad_simplification/voxelization_simplify.py",
    "cad_simplification/triple_extrusion_simplify.py",
    "cad_simplification/triangle_decimation_simplify.py",
    "cad_simplification/alpha_wrap_simplify.py",
    # Voxelization
    "discrete_representation/voxelization/interference/interference_surface_showcase.py",
    "discrete_representation/voxelization/interference/interference_volume_showcase.py",
    "discrete_representation/voxelization/step_file_voxelization.py",
    "discrete_representation/voxelization/stl_file_voxelization.py",
    "discrete_representation/voxelization/volume_model_matrix_based_voxelization.py",
    "discrete_representation/voxelization/volume_model_point_based_voxelization.py",
    "discrete_representation/voxelization/voxelization_filling.py",
    "discrete_representation/voxelization/voxelization_inverse.py",
    "discrete_representation/voxelization/voxelization_moving.py",
    "discrete_representation/voxelization/octree/inner_growing_showcase.py",
    "discrete_representation/voxelization/octree/non_homogeneous_showcase.py",
    "discrete_representation/voxelization/octree/intersecting_faces_combinations.py",
    # Pixelization
    "discrete_representation/pixelization/pixelization_boolean_operations.py",
    "discrete_representation/pixelization/pixelization_filling.py",
    "discrete_representation/pixelization/pixelization_from_closed_polygon.py",
    "discrete_representation/pixelization/pixelization_inverse.py",
    # Others
    "grid.py",
]


# Maximum time for a script
CONTROLED_TIMES = {"showcases/casing.py": 15, "primitives/sweep.py": 15}

# Executing scripts
print("Executing scripts for CI:")
total_time = time.time()
top_level_dir = os.getcwd()
# top_level_dir = os.sep.join(__file__.split(os.sep)[:-1])

times = {}
for script_name in scripts:
    if not os.path.isfile(os.path.join(top_level_dir, script_name)):
        raise FileNotFoundError(f"Script {script_name} does not exists in CI scripts")

    print(f"\t* {script_name}")
    # Reset dir
    os.chdir(top_level_dir)
    # Change cwd
    if "/" in script_name:
        script_folder = "/".join(script_name.split("/")[:-1])
        if script_folder:
            script_folder = os.path.join(top_level_dir, script_folder)
            os.chdir(script_folder)
    file_name = script_name.split("/")[-1]
    time_start_script = time.time()

    with open(file_name, "r", encoding="utf-8") as script:
        # Call babylon_data instead of babylonjs to avoid opening browser windows
        pattern = r'\.babylonjs\(.*?\)'
        replacement = r'.volmdlr_volume_model().babylon_data()'
        exec(re.sub(pattern, replacement, script.read()))

    time_start_script = time.time() - time_start_script
    times[script_name] = time_start_script
    _plt.close("all")

print("Computation times:")
for script_name, t in sorted(times.items(), key=lambda x: x[1]):
    print(f"* script {script_name}: {round(t, 3)} seconds ")
    if script_name in CONTROLED_TIMES and t > CONTROLED_TIMES[script_name]:
        raise RuntimeError(f"This script {script_name} should take less than {CONTROLED_TIMES[script_name]}")

total_time = time.time() - total_time
print(f"Total time for CI scripts: {total_time}")
