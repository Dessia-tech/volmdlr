

scripts = [
            # Core.py
            'core/points.py',
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
            # Steps
            'read_steps.py',
            # Stl
            'stl_reading.py',
            # Distance
            'distance/arc3D_arc3D.py',
            'distance/arc3D_ls3D.py',
            # Showcases
            'showcases/simple_shapes.py',
            'showcases/casing.py'

            ]
#  'cyl_cyl.py', 'cyl_pf.py',
# 'ls3D_ls3D.py', 'sweep_sweep.py', 'tore_cyl.py','tore_pf.py'
# 'tore_tore.py'

for script_name in scripts:
    print('\n## Executing script {}'.format(script_name))

    exec(open(script_name).read())
