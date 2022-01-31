
import os

scripts = [
            # Edges
            'edges/arcs2D.py', 'edges/arcs3D.py',
            # Wires
            'wires/roundedlines.py','wires/polygon2D.py',
            'wires/triangle2D.py',
            'demo2D.py',
            # Primitives
            'primitives/extrusion.py', 
            'primitives/sweep.py',
            'primitives/revolved_profile.py', 'edges/areas_moment_cog_check.py',
            'primitives/block3d.py',
            'primitives/sphere_to_point.py',
            # Faces
            'faces/triangle3d.py',
            'faces/bspline.py',
            'faces/bspline_bark.py',
            'faces/union.py',
            'faces/blocks_set_operations.py',
            # Cloud
            'cloud/sewing_two_polygons.py',
            'cloud/sewing_stl.py',
            # Steps
            'read_steps.py',
            # Stl
            'stl/stl_reading.py',
            # Distance
            'distance/arc3D_arc3D.py','distance/arc3D_ls3D.py',
            # Showcases
            'showcases/simple_shapes.py',
            'showcases/casing.py',
            ]
#  'cyl_cyl.py', 'cyl_pf.py',
# 'ls3D_ls3D.py', 'sweep_sweep.py', 'tore_cyl.py','tore_pf.py'
# 'tore_tore.py'

for script_name in scripts:
    print('\n## Executing script {}'.format(script_name))

    exec(open(script_name).read())
