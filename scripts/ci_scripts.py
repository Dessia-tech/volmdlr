
import os

scripts = [
            'edges/arcs2D.py', 'edges/arcs3D.py', 'showcases/simple_shapes.py',
            'wires/roundedlines.py','wires/polygon2D.py',
            'wires/triangle2D.py',
            'primitives/extrusion.py', 'demo2D.py',
            'showcases/casing.py',
            'primitives/sweep.py',
            'primitives/revolved_profile.py', 'edges/areas_moment_cog_check.py',
            'distance/arc3D_arc3D.py','distance/arc3D_ls3D.py',
            'primitives/block3d.py',
            'faces/triangle3d.py',
            'faces/union.py',
            'faces/blocks_set_operations.py',
            'cloud/sewing_two_polygons.py',
            'read_steps.py',
            'cloud/sewing_stl.py',
            'primitives/sphere_to_point.py'
            ]
#  'cyl_cyl.py', 'cyl_pf.py',
# 'ls3D_ls3D.py', 'sweep_sweep.py', 'tore_cyl.py','tore_pf.py'
# 'tore_tore.py'

for script_name in scripts:
    print('\n## Executing script {}'.format(script_name))

    exec(open(script_name).read())
