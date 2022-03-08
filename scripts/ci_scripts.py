import os
import time

scripts = [
            # Core.py
            'core/points.py',
            'core/grid2d_with_direction.py',
            # geometry
            'geometry.py',
            # Edges
            'edges/arcs2D.py',
            'edges/arcs3D.py',
            'edges/bspline.py',
            'edges/bspline2.py',
            'edges/areas_moment_cog_check.py',
            'edges/lines2d.py',
            # Wires
            'wires/roundedlines.py',
            'wires/polygon2D.py',
            'wires/triangle2D.py',
            # Primitives
            'primitives/extrusion.py',
            'primitives/sweep.py',
            'primitives/revolved_profile.py',
            'primitives/block3d.py',
            'primitives/sphere_to_point.py',
            'primitives/cone.py',
            'primitives/cylinders.py',
            # Faces
            'faces/triangle3d.py',
            'faces/bspline.py',
            'faces/bspline_bark.py',
            'faces/union.py',
            'faces/blocks_set_operations.py',
            'faces/BSplineSurface/bspline_curves.py',
            'faces/BSplineSurface/bspline_curves_point_belongs.py',
            'faces/BSplineSurface/bspline_surface_interpolation.py',
            'faces/BSplineSurface/bspline_surface_merge.py',
            'faces/BSplineSurface/bspline_surface_split.py',
            'faces/BSplineSurface/bspline_surface_to_2d.py',
            'faces/BSplineSurface/from_cylindrical_to_bspline_surface.py',
            'faces/BSplineSurface/bspline_surface_definition.py',
            'faces/BSplineSurface/bspline_surfaces_grid3d.py',
            # Cloud
            'cloud/sewing_two_polygons.py',
            'cloud/sewing_stl.py',
            # Contours
            'contours/contour_split.py',
            'contours/contour_merge_with_1.py',
            'contours/contour_merge_with_2.py',
            'contours/cut_by_line.py',
            'contours/contour_cut_by_wire.py',
            # Steps
            'step/read_steps.py',
            # Stl
            'stl/read_stls.py',
            # Distance
            'distance/arc3D_arc3D.py',
            'distance/arc3D_ls3D.py',
            # Showcases
            'showcases/simple_shapes.py',
            'showcases/casing.py'
            ]

# Testing if all scripts exists before launching them
for script_name in scripts:
    if not os.path.isfile(script_name):
        raise FileNotFoundError(f'Script {script_name} does not exists in CI scripts')

# Executing scripts
print('Executing scripts for CI:')
total_time = time.time()
top_level_dir = os.getcwd()
times = {}
for script_name in scripts:
    print(f'\t* {script_name}')
    # Reset dir
    os.chdir(top_level_dir)
    # Change cwd
    if '/' in script_name:
        script_folder = '/'.join(script_name.split('/')[:-1])
        if script_folder:
            script_folder = os.path.join(top_level_dir, script_folder)
            os.chdir(script_folder)
    file_name = script_name.split('/')[-1]
    t = time.time()
    with open(file_name, 'r', encoding='utf-8') as script:
        exec(script.read())
    t = time.time() - t
    times[script_name] = t

print('Computation times:')
for script_name, t in sorted(times.items(), key=lambda x:x[1]):
    print(f'* script {script_name}: {round(t, 3)} seconds ')
    
total_time = time.time() - total_time
print(f'Total time for CI scripts: {total_time}')
