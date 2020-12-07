
import os

scripts = ['arcs2D.py', 'arcs3D.py', 'block3d.py', 'simple_shapes.py',
           'roundedlines.py','polygon2D.py', 'polygon2D_2.py',
           'extrusion.py', 'demo2D.py', 'casing.py', 'sweep.py',
           'revolved_profile.py']

for script_name in scripts:
    print('Executing script {}'.format(script_name))

    exec(open(script_name).read())

distance_scripts = ['arc3D_arc3D.py','arc3D_ls3D.py', 'cyl_cyl.py', 'cyl_pf.py',
                   'ls3D_ls3D.py', 'sweep_sweep.py', 'tore_cyl.py','tore_pf.py'
                   'tore_tore.py']

for script_name in distance_scripts:
    print('Executing script {}'.format(script_name))
    exec(open(os.path.join('distance', script_name)))
