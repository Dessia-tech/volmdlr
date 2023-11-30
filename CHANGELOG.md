# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.16.0 [future]

### New Features
-

### Fixed

#### curves.py
- Ellipse2D/3D: mutualize length method.

#### edges.py
- BSplineCurve: handles exceptions in simplify method.
- BSplineCurve: Consider overlaping curves also as periodic.
- BSplineCurve.simplify: handles exceptions.
- 
#### surface.py
- PeriodicalSurface: handles exceptions in connect_contours method.
- ExtrusionSurface3D: fullarcellipse3d_to_2d
- ExtrusionSurface3D: generalization of the _repair_points_order method to repair the order of parametric points of edges after transformation.
- ToroidalSurface3D: increases precision of point3d_to_2d.

### Refactor
- Big refactor to improve and simplify complex and long methods in various modules.

#### surfaces.py
- contour3d_to_2d: It also returns a dictionary with the correspondence between the parametric and 3D primitives.

### Changed
- Edge.split_between_two_points -> trim

### Unittests
- 

## v0.15.0

### New Features

#### core_compiled.py
- Point2D/Point3D: allow users to use a point or a list of points direct inside a numpy array. ex.: np.array(volmdlr.O3D)
- Point2D/Point3D: in_list. ** ATTENTION:** -> use in_list instead of volmdlr.core.point_in_list.
- cad_simplification: VoxelizationSimplify, TripleExtrusionSimplify, TriangleDecimationSimplify.

#### surfaces.py
- ToroidalSurface3D: line_intersections, linesegment_intersections, plane_intersections
- ToroidalSurface3D: cylindricalSurface_intersections, circle_intersections, fullarc_intersections, dict_to_object, conicalsurface_intersections, sphericalsurface_intersections
- ToroidalSurface3D: Handles degenerated surfaces (major_radius < minor_radius).
- CylindricalSurface3D: circle_intersections, sphericalsurface_intersections, cylindricalsurface_intersections
- ToroidalFace3D: PlaneFace3D intersectios.
- SphericalSurface3D: circle_intersections, arc_intersections, ellipse_intersections, arcellipse_intersections, sphericalsurface_intersections
- ConicalSurface3D: sphericalsurface_intersections
- General improvements on sufaces\` parametric operations.

#### edges.py
- BsplineCurve3D: circle_intersections.
- ArcEllipse3D/FullArcEllipse3D: line_intersections.
#### curves.py
- Circle3D: point_distance.
#### shell.py
- OpenTriangleShell3D: triangle decimation
- ClosedTriangleShell3D: turn_normals_outwards, are_normals_pointing_outwards, turn_normals_inwards, are_normals_pointing_inwards
- DiplayTriangleShell3D: concatenate

#### core.py
- BoundingBox: is_close, scale
- BoundingBox: triangle_intersects_voxel, is_intersecting_triangle
#### discrete_representation.py
- Voxelization: from_mesh_data
- OctreeBasedVoxelization


#### step.py
- Support to Datakit CrossCadWare STEP file format.

### Fixed
- Drone : run generate sdist and generate bdist_wheel only on master

#### core.py
- VolumeModel: get_mesh_lines (change tolerance to 1e-5)

#### edges.py 
- Arc2D: direction conservation in rotation / translation / frame_mapping.

#### surfaces.py
- ToroidalSurface3D: line_intersections, linesegment_intersections, plane_intersections 
- ConicalSurface3D: circle_generatrixes direction.

#### faces.py
- ToroidalFace3D: PlaneFace3D intersections.
- PlaneFace3D: circle_intersections. planeface_intersections
- BsplineFace3D: adjacent_direction_uu
- PlaneFace3D: project_faces (check first if surfaces are coincident)

#### wires.py
- delete remaining inplace methods in wires.py

#### shells.py
- Fixes to boolean operations. Added some tolerance parameters to some methods. 
- Shell3D: get_geo_lines (consider edge.inverse in get_edge_index_in_list check), is_shell_open
- DisplayTriangleShell3D: eq, data_eq, hash, data_hash, concatenate

#### surfaces.py 
- SphericalSurface3D: use circle 3d instead of polygon3D for plotting. 

#### utils
- common_operations separate_points_by_closeness: consider more than two cluster groups.

#### curves.py
- Circle3D: circle_intersectios when the circle are coplanar.
- Circle2D: Now, it needs a Frame2D and a radius instead of a Center and a Radius. This allows to easily control the circle's direction (clockwise/counterclockwise)

#### surfaces.py
- ExtrusionSurface3D: enhance parametric operations.

#### edges.py
- bsplineCurve: line_intersections. 

#### discrete_representation.py
- MatrixBasedVoxelization: _logical_operation
- Remove inheritance from ABC for platform usage

#### cad_simplification.py
- Remove inheritance from ABC for platform usage

### Refactor
- Face3D: create a generic method for calculating intersections between two faces: _generic_face_intersections.
- Voxelization: refactor class methods

#### core.py
- babylon_data: avoid using bounding_box for performance
- BoundingBox: uses numpy to improve performance.

#### core_compiled
- Frame2D: fix rotation, now it has an optional parameter rotate_basis, set to False by default option, so the user can specify if he wants to rotate also the basis of the frame.

#### edges.py
- Circle2D: Now, it needs a Frame2D and a radius instead of a Center and a Radius. This allows to easily control the circle's direction (clockwise/counterclockwise)
- Arc2D: Arc 2D now must follow the same rotation direction of its circle.
- LineSegment2D/3D: The line attribute from which the line segment was defined was converted to a property, for performance and memory efficiency reasons.
- BSplineCurve: improve line_intersections performance.

#### faces.py
- Face3D: create a generic method for calculating intersections between two faces: _generic_face_intersections.

#### primitives3d.py
- Sweep: accepts an optional parameter starting_frame that can control the orientation of the profile.
- Block: get_bounding_box

#### shells.py
- boolean operations - now works also for triangle meshed objects, containing coincident faces.
#### surfaces.py
- ExtrusionSurface3D: Uses edge abscissa as u parameter.
- ExtrusionSurface3D: general improvements in parametric operations.


### Changed
- ToroidalSurface3D: init param tore_radius and small_radius changed to major_radius and minor_radius respectevely.
- ToroidalSurface3D: plots now use Circles 3D instead of ClosedPolygon3D. Performance improved.
- CylindricalSurface3D: More comprehesive plot
- BoundingBox: from_bounding_boxes
- BSplineCurve: improve line_intersections performance.
- core_compiled.pyx: update typing because Point2D, Point3D, Vector2D and Vector3D are now extension types (C structures.)
- BSplineCurve: improve line_intersections performance.
- SphericalSurface3D: enhance bsplinecurve3d_to_2d.


### Unittests
#### curves 
- Circle3D: new case to test_circle_intersections, new test: test_point_distance.
#### surfaces
- ToroidalSurface3D: test_line_intersections, test_plane_intersections, test_cylindrical_surface_intersections, test_circle_intersections
- CylindricalSurface3D:  test_circle_intersections.
#### faces
- ToroidalFace3D: PlaneFace3D intersectios.
- SphericalSurface3D: circle_intersections, arc_intersections, arcellipse_intersections
- PlaneFace3D: point_belongs
#### core
- BoundingBox: is_close, scale
#### primitives3d
- Block: from_bounding_box, get_bounding_box

## v0.14.0

### New Features
- DisplayTriangleShell3D: a TriangleShell3D optimized for performance of display / saving / loading.
- BSplineSurface3D: from_points_interpolation, from_points_approximation.
- nurbs module.
- New curves classes: Hyperbola2D and Hyperbola3D.
- Line: closest_point_on_line, from_point_and_vector
- Line2D: get_slope, get_y_intersection.
- New curves classes: Parabola2D/3D.
- ConicalSurface3D: line/line_segment intersections, perpendicular_plane_intersection
- ConicalSurface3D: line/line_segment intersections, perpendicular_plane_intersection, parallel_plane_intersections, concurent_plane_intersections, plane_intersections.
- Hyperbola2D/3D and Parabola2D/3D: split
- PlaneFace3D: conicalface_intersections
- CylindricalSurface3D: conicalsurface_intersections
- CylindricalFace3D: conicalface_intersections
- Curve: general_method curve_intersections
- Parabola2d/3D / Hyperbola2D/3D: point_belongs, tangent
- BSplineCurve: point_to_parameter, abscissa_to_parameter.
- Basis3D: is_normilized, is_orthogonal, is_orthonormal.
- BSplineSurface3D: fullarcellipse3d_to_2d
- ClosedPolygon2D: points_in_polygon

### Fixed
- add missing name attributes to classmethods.
- fixed circular imports
- BSplineSurface3D: from_points_interpolation, from_points_approximation.
- ConicalFace3D: point_belongs
- nurbs.core: find_multiplicity, evaluate_curve.
- LineSegment3d: line_intersections.
- Circle2D: line_intersections
- Step.read_lines: handles name with # character in name.
- ExtrusionSurface3D: enhance 3D to parametric operations.
- BSplineCurve: direction_vector, point_at_abscissa, abscissa, trim
- ConicalSurface3D and RevolutionSurface3D: bsplinecurve3d_to_2d when start or and points are at surface singularity
- ClosedCurves: discretization_points
- ArcEllipse3D: is_close
- LineSegment3D: revolution
- FullArcEllipse3D, FullArcEllipse2D: discretization_points
- ConicalSurface3D: linesegment2d_to_3d
- BSplineSurface3D: bsplinecurve3d_to_2d, prevents code execution from stopping when point3d_to_2d does not converge
- BSplineSurface3D: derivatives
- BSplineCurve: split
- Matrix based discrete representation: boolean operations
- read the docs settings
- fix: move code complexity at end
- ClosedPolygon2D: points_in_polygon, fix include_edge_points
- ClosedShell3D: is_face_between_shells

### Refactor
- TriangleShell3D: various improvement such as get_bounding_box, to_mesh_data, from_mesh_data, to_dict, dict_to_object

### Changed
- Cache BSplineCurve points into a numpy array to reduce memory usage.
- Vector2D, Vector3D: __repr__
- core_compiled: cdef functions' names.
- Vector2D, Vector3D, Point2D, Point3D: transformed into extension types for memory performance
- limit warning on step reading
- BSplineSurface3D: point3d_to_2d

### Unittests
- Hyperbola2D/3D: line_intersections
- Parabola2D/3D: line_intersections
- ConicalSurface3D: test_line_intersections, test_plane_intersections.

## v0.13.0

### New Features
- Line: reverse.
- BSplineCurve: Remove dependencies from the geomdl library.
- perf: to_dict/dict_to_obj of OpenTriangleShell3D
- Cylinder / Cone / HollowCylinder: from_center_point_and_axis
- Cone: remove inheritance from RevolvedProfile
- Ellipse2D: point_distance, bounding rectangle, ellipse_intersections
- Curve: local_discretization
- Ellipse3D: line_intersections, linesegment_intersections, ellipse_intersections
- ArcEllipse3D : Linesegment_intersections, arcellipse_intersections
- Circle3D: circle_intersections, ellipse_intersections
- Circle2D: ellipse_intersections.
- Arc3D: arc_intersections, arcellipse_intersections
- Wire3D/Contour3D: edge_intersections, wire_intersections
- BSpline3D: arc_intersections
- New module: discrete_representation for voxelization of 3D geometries and pixelization of 2D geometries
- BSplineSurface3D: partial removal of dependencies on geomdl objects

### Fixed
- Sweep with non smoth path
- plot of vector3D.
- BSplineSurface3D: point3d_to_2d, improve inital condition.
- EdgeCollection3D: babylon_meshes.
- BSplineCurve3D: trim
- FullArc3D: hash
- SphericalSurface3D: enhance repair_periodicity_method
- CylindricalSurface3D: concurrent_plane_intersection
- BSplineFace3D: fix neutral_fiber
- Step: assembly import
- BSplineFace3D: fix bounding_box.
- Ellipse3D: from_step
- edges.py: general improvements.
- ExtrusionSurface3D: point3d_to_2d.
- ExtrusionSurface3D: enhance parametric operations when the surface is periodic.
- BSplineFace3D: fix neutral_fiber
- BSplineSurface3D: improve bsplinecurve3d_to_2d.
- BSplineSurface3D: improve bsplinecurve3d_to_3d.
- Circle2D: plot
- Line3D: fix Line3D plot()
- Vector2D: plot()
- fix RevolutionFace3D init parameter wire to edge.
- Update documentation
- fix Sweep: bug when first primitive is an arc.
- fix closedshell3d volume
- Step.py: enhance step import/export
- VolumeModel: get_shells
- step.py uses deque in stack based algorithms
- VolumeModel: get_shells
- add error protection stl
- Sweep - add raise ValueError if section too big in comparision to arc radiuses
- Update cython version requirement in setup.py
- Step import: handles when there is an empty assembly in the file.
- Ellipse2D: point_at_abscissa
- ultis.common_operations: get_edge_distance_to_point and get_get_abscissa_discretization from edges so it can be used in curves too.
- edges.Edge._generic_minimum_distance
- LineSegment3D: distance_linesegment
- BSpline3D: linesegment_intersections

### Refactor
- refator some classes' init in primitives3D. 
- Shells: refactor.
- Composite_primitives
- Surface3D: enhance repair_primitives_periodicity method.
- volmdlr.utils.intersections:
- BSplineCurve: replace periodic bool parameter with verification inside from_points_intepolation method.
- Wire3D: removes heritage from volmdlr.core.CompositePrimitive3D
- BSplineCurve3D: bounding_box
- edges: minimum_distance.
- BSplineSurface3D: bsplinecurve3d_to_2d
- BSplineCurve: transform some attributs into lazy evaluation and Caching
- BSplineSurface3D: transform some attributs into lazy evaluation and Caching
- BSplineSurface3D: store control_points as numpy array for memory efficiency
- PlaneFace3D: distance_to_point -> point_distance
- remove normalize() methods for Vectors. Replaced by unit_vector(), it returns a new normalized vector.
- Cylinder / Cone / HollowCylinder: docstrings, typings, style, coherence
- BSplineSurface3D: point3d_to_2d performance improvements.


### Changed
- Moves functions from step.py to volmdlr.utils.step_reader
- Cylinder / HollowCylinder: `from_extremal_points` is now depracted. Use `from_end_points` instead (for lexical reason)

### Unittests
- Cylinder / Cone / HollowCylinder
- Ellipse2D: point_distance
- Ellipse3D: test_ellipse_intersections, test_linesegment_intersections
- ArcEllipse3D : Linesegment_intersections, arcellipse_intersections
- Circle3D: circle_intersections.
- Arc3D: arc_intersections, arcellipse_intersections, test_minimum_distance_bspline
- BSplineCurve3D: test_bspline_linesegment_minimum_distance, test_bspline_linesegment_intersections
- Contour3D: test_edge_intersections

## v0.12.0


### New Features
- New module: cad_simplification - OctreeBlockSimplify, TrippleExtrusionSimplify
- shells.py : function to performe union operations for a given list of shells.
- ClosedShell3D: is_face_intersecting, is_intersecting_with
- BoundingBox: get_points_inside_bbox, size
- Vector3D: unit_vector
- Face3D: split_inner_contour_intersecting_cutting_contours
- Shell3D: get_ray_casting_line_segment
- WireMixin: get_connected_wire, is_sharing_primitives_with
- OpenShell3D: faces_graph
- Plane3D: arc_intersections, bsplinecurve_intersections
- common_operations: split_wire_by_plane
- SphericalSurface3D: line_intersections, linesegment_intersections.
- Sweep with muitiform profile contour.
- minimum_distance: face-to-face, shell-to-shell
- OpenShell3D: from_faces (using faces graph)
- SphericalFace3D: from_contours3d_and_rectangular_cut
- RevolutionSurface3D: Translation
- wires.WireMixin: from_circle
- curves.CircleMixin: trim
- Face3D: point_distance
- BSplineCurve3D: revolution method.

### Fixed
- ClosedShell3D: is_face_inside, get_subtraction_valid_faces, valid_intersection_faces, point_belongs
- ContourMixin: delete_shared_contour_section, reorder_contour_at_point, are_extremity_points_touching
- RevolutionSurface3D: fix some special cases whiling transforming from 3D space to parametric domain.
- fix drone python version
- BSplineFace3D: neutral_fiber
- BSplineSurface3D: arc3d_to_2d, removes repeated parametric points if any.
- surfaces.Plane3D: linesegment_intersections
- Step export
- Face3D: is_linesegment_crossing.
- Edge: fix orientation of edges commig from step.
- BSplineCurve3D: from_step.
- Export to step file
- Step import
- Edge: fix orientation of edges commig from step.
- Sphere: point_belongs, inherits from ClosedShell3D instead of RevolvedProfile
- Step import.
- PeriodicalSurface: linesegment3d_to_2d, takes into account small 3D line segments that should be actually 3D arcs
- babylondata: removes empty objects.
- ClosedPolygon2D: point_belongs.
- Fullarc: get_reverse.
- Arc2D: point_belongs
- ArcEllipse2D: point_at_abscissa
- Frame3D: import/export step.
- BSplineFace3D: neutral_fiber.
- Step: read_lines, take into account the space character in step entity names
- Circle3D: fix trim.
- Edge: from_step trim of periodic curves with different orientation of original edge
- Arc3D: fix abscissa, fix get_arc_point_angle
- add missing toleraces to some methods.
- Arc3D: line_intersections
- Line3D: minimum_distance_points
- remove arcellipse handleling for bspline2d_3d.
- plot of vector3D
- Ellipse3D: discretization_points.

### Refactor
- ClosedShell3D: point_belongs, get_non_intersecting_faces
- BoundingBox: bbox_intersection
- Face3D: get_face_cutting_contours
- parametric.py: fix numerical instability in some functions used in Arc3D to parametric surface domain transformation.
- intersections: get_bsplinecurve_intersections generalization, so it can also be used
to calculate intersections between a plane 3d and bsplinecurve3d.
- Big refactor: New module curves.py containing classes as Line, Circle and Ellipse.
Most edges will now be formed by a curve and a start and end points. Unittests for all these classes have been created.
All adequations have been done for all tests and existing scripts.

- bspline_compiled: refactor binomial_coefficient for performance.
- Improve step translator.
- Delete inplace methods: rotation, translation and frame_mapping. replace by juste the rotation, translation and frame_mapping. objects are no longer changed inplace, a new transformed object is returned each time.
- OpenShell3D: faces_graph.
- RevolutionSurface3D: Improve init and methods

### Changed
- OpenShell3D: faces_graph is now vertices_graph. faces_graph method now represents the faces' topology of the shell.

### Unittests
- FullArc2D: split_between_two_points
- Face3D: set_operations_new_faces
- ClosedShell3D: point_belongs
- Plane3D: arc_intersections, bsplinecurve_intersections
- common_operations: split_wire_by_plane
- SphericalSurface3D: line_intersections, linesegment_intersections.

## v0.11.0


### New Features
- BSplineCurve, Edge: simplify
- Plane3D: angle_between_planes, plane_betweeen_two_planes
- Edge: intersections, crossings, validate_crossings
- Arc2D: bsplinecurve_intersections, arc_intersections, arcellipse_intersections.
- ArcEllipse2D: bsplinecurve_intersections
- get_circle_intersections added to volmdlr.utils.intersections, so it can be used to calculate intersections between two arcs 2d.
- get_bsplinecurve_intersections added to volmdlr.utils.intersections. Used to calculate intersection between a bspline and another edge.
- Wire2D: edge_intersections, wire_intersections, edge_crossings, edge_intersections, validate_edge_crossings, validate_wire_crossings
- Contour2D: split_contour_with_sorted_points, intersection_contour_with
- CylindricalSurface3D: point_projection, point_distance
- ToroidalSurface3D: point_projection
- BsplineCurve: point_distance, point_belongs
- ContourMixin: is_adjacent
- Wire2D: area
- Circle2D: bsplinecurve_intersections.
- add tolerance param to many methods from edges and wires.
- Surface3D: add contour healing into face_from_contours3d method.
- ExtrusionSurface3D: implement missing cases for linesegment2d_to_3d method.
- BSplineSurface3D: to_plane3d
- BSplineFace3D: to_planeface3d
- BSplineCurve, Arc, LineSegment: is_close
- Core: get_edge_index_in_list, edge_in_list
- mesh: TetrahedralElementQuadratic 
- GmshParser: define_quadratic_tetrahedron_element_mesh
- GmshParser: to_vtk (consider quadratic tetrahedron element)
- VolumeModel: to_msh (consider both order 1 and 2)
- Assembly: define a volmdlr Assembly object.
- Edge: direction_independent_is_close
- Arcellipse2D, 3D: complementary, translation
- Arcellipse2D, 3D: complementary
- Face3D: is_linesegment_crossing, linesegment_intersections_approximation.
- Assembly: define a volmdlr Assembly object.
- Contour2D: copy
- LineSegment2D: copy
- FullArcEllipse3D: split
- ArcEllipse3D: split, point_at_abscissa
- Vector: is_perpendicular_to
- babylonjs: add nested meshes
- CylindricalFace3D, ConicalFace3D, ToroidalFace3D, BSplineFace3D: neutral_fiber
- VolumeModel: get_shells
- WireMixin: wires_from_edges
- DisplayMesh3D: triangulation_faces
- Woodpecker CI setup
- ContourMixin: primitive_section_over_contour.
- Face3D: split_by_plane

### Fixed
- 2D conversion: create 2D function name in core_compiled
- LineSegment, Arc, BSplineCurve: get_shared_section()
- bSpline2D: linesegment_intersections
- BsplineCurve: from_points_interpolation
- Coverage: use coverage rc to enable cython coverage
- ClosedShel3D: cut_by_plane
- ClosedShell3D: union
- BSplineSurface3D: take into account oppened contour while using face_from_contours3d
- BsplineCurve: simplify
- Dessiaobject inheritance up-to-date
- Edge: unit_direction_vector, unit_normal_vector, split_between_two_points
- VolumeModel: get_mesh_lines (change tolerance 1e-20 to 1e-6)
- RevolutionSurface: fix some parametric operations.
- ClosedShel3D: intersection method
- Fix: plots
- add some fixes to pydocstyle errors
- ToroidalSurface3D: fix some parametric operations.
- Node2D, Node3D: is_close
- SphericalSurface3D: enhance arc3d_to_2d and bsplinecurve3d_to_2d.
- BSplineface3D: linesegment2d_to_3d, bsplinecurve2d_to_3d.
- OpenShell3D: get_geo_lines (use primitive.is_close)
- Basis3D: normalize
- Contour3D: from_step removes repeated edges from primitives list
- Face3D: add fixes to divide_face.
- ExtrusionSurface3D: linesegment2d_to_3d.
- Surface3D: repair_primitive_periodicity
- BSplineSurface3D: ban useless attr in serialization 
- utils.parametric: fix contour2d_healing
- BSplineSurface3D: ban useless attr in serialization
- BSplineCurve: simplify
- SphericalSurface3D: contour3d_to_2d
- WireMixin: to_wire_with_linesegments (use new methods, for 2D and 3D)
- ArcEllipse2d: point_belongs, abscissa, init.
- Face3D: face_inside - now considers inners_contours
- BoundingBox: point_belongs now considers bounds.
- ContourMixin: delete_shared_contour_section
- PlaneFace3D: merge_faces
- Contour2D: divide
- Step: raise NotimplementedError when it's not possible to instatiate assembly object.
- STL: handle mutiple space as separator
- fix: protect gmsh import

### Refactor
- Contour2D: cut_by_wire
- Contour2D: extract_with_points displaced to WireMixin
- Contour2D: extract_contour displaced to WireMixin and renamed to extract
- Contour2D: split_contour_with_sorted_points displaced to WireMixin and renamed to split_with_sorted_points
- Contour2D: get_divided_contours
- FullArc2D, FullArc3D: create FullArc Abstract class.
- Contour2D: ordering_contour
- WireMixin: order_wire
- Contour2D: delete cut_by_linesegments
- split faces.py into surfaces.py, faces.py and shells.py 
- ContourMixin: from_points
- ClosedShell3D: improve performance for boolean operations
- Face3D: reduce the triangulation discretization resolution of Toroidal and Cylindrical to improve redering performance.
- Cylinder: inheritance directly from ClosedShell3D
- Edges: cache middle_points and unit_direction_vector 
- Arc: add optional parameter center
- unittests: find dynamicly the folder for the json
- Arc: point_distance
- BSplineCurve: is_close
- CompositePrimitive3D: babylon_points
- WireMixin: split_with_sorted_points -> if a wire, and given points are start and end, return self directly.
- ContourMixin: contours_from_edges
- ExtrusionSurface3D: simplify bsplinecurve3d_to_2d method

### Changed
- better surface3d plots
- sphere methods renamed in_points & to_point_skin to inner points & skin_points
- Improve CylincricalFace3D and ToroidalFace3D rendering mesh.
- remove useless attribute in Bspline serialization
- Change python suport version from >=3.7 to >= 3.9
- LICENSE changed from GPL to Lesser GPL 
- Readme logo updated
- CI: do not check quality on tag

### Unittests
- Arc2D: test_arc_intersections
- TestEdge2DIntersections: test intersections for all edges.
- Circle2D: test_circle_intersections
- Contour2D: test_crossings, test_intersection_contour_with
- BSplineCurve: get_intersection_sections
- BSplineCurve2D: edge_intersections, arc_intersections, bsplinecurve_intersections
- CylindricalFace3D: test_triangulation_quality
- CylindricalSurface3D: test_point_projection
- BSplineCurve: point_projection
- ClosedShel3D: cut_by_plane
- Arc3D.minimum_distance_points_line
- New unittests for plane3d.
- ClosedShel3D: intersection
- Arcellipse2D: complementary
- Contour2D: contours_from_edges.
- PlaneFace3D: merge_faces
- Contour2D: divide.
- BSplineFace3D: test_linesegment_intersections_approximation.
- CylindricalFace3D: split_by_plane.

v0.10.0 [Released 20/04/2023]

### New Features
* Write .msh file (with stream)
* Arc: reverse
* BSplineCurve2D: offset
* Circle2D: bsplinecurve_intersections, point_distance
* ConicalSurface3D, CylindricalSurface3D: plot method
* BSplineCurve3D: minimum distance
* volmdlr.edge: FullArcEllipse
* BSplineCurve: evaluate_single
* Wire2: hash
* Contour3D: hash
* LineSegment3D, LineSegment2D, Arc3D, Arc2D, BSpline3D, BSpline2D: get_shared_section(), delete_shared_section()
* Contour2D: closest_point_to_point2, get_furthest_point_to_point2
* Block: octree, quadtree, subdivide_block

### Fixed
* Bspline in sweep
* Plane3D: plane_intersections
* fixes to step assemblies
* LineSegment3D: matrix_distance
* fixes to wire
* Arc: split. Case when spliting point is the start or end point.
* BplineCurve2D: tangent, vector_direction, normal_vector
* BSplineCurve: abscissa, line_intersections
* Add some important fixes to unittests: missing two __init__py files.
* Contour2D, Contour3D: merge_with()
* Edge: change unit_direction_vector and unit_normal_vector to concrete methods
* stl: add _standalone_in_db to Stl class
* BSplineSurface3D: merge_with
* Documentation: Add introduction to volmdlr technology
* BSplineSurface3D: refactor bsplinecurve3d_to_2d to take into account periodic behavior
* OpenedRoundedLineSegments2D/ClosedRoundedLineSegments2D: fix radius type
* Surface3D: debug some special cases while using face_from_contours3d.
* Step: debug some special cases while reading step file.
* BSplineSurface3D: fix simplify_surface method.
* Improve pylint code quality.
* PeriodicalSurface: enhance some parametric transformations.

### Removed
- stl: remove default value in from_stream method

### Changed

- argument convexe in volmdlr.cloud has been renamed to convex
- Add some missing docstrings in volmdlr.faces
- Using full arcs for Circles primitives

### Performance improvements
- BSplineCurve: compilation of some functions used by from_points_interpolation classmethod.
- BSplineSurface3D: compilation of some functions used in the evaluation of a parametric point.
- eq & hash: Some eq and hash methods have been fixed. starting from clases Point and Vector.
- BSplinecurve2D: point_belongs
- lighten some dicts with optional name
- Step reader: refactor to_volume_model. Remove the dependency of the method of creating a graph.

### Refactorings
- ContourMixin: to_polygon (for both 2D and 3D)
- BSplineCurve2D.point_distance 
- new dataclass EdgeStyle: to be used in several plot methods. simplifying its structure.


### Unittests
* BSplineCurve2D: offset, point_distance, point_belongs
* Circle2D: bspline_intersections, point_distance
* Unittests for Vector2D
* Unittests for Point2D
* Unittests for Vector3D
* Unittests for Point3D
* LineSegment3D: test_matrix_distance
* LineSegment3D, LineSegment2D, Arc3D, Arc2D, BSpline3D, BSpline2D: get_shared_section(), delete_shared_section()
* Contour3D: merge_with()
* Contour2D: closest_point_to_point2, get_furthest_point_to_point2

## v0.9.3

- build: bump dessia common to 0.10.0
- build: remove useless jsonschema dep
- build: update package.xml for freecad

## v0.9.1

### Fixed
- build: manifest was not shipping bspline_compiled
- fixed many pylint errors: 13/03/2023
- fix contour2d: divide

### Documentation
 - typo in CONTRIBUTING.md
 - typo in README.md

## v0.9.0 [released 03/26/2023]

### New Features
* Unit coversion factor parameter added to the end of the from_step arguments parameter (So we can convert the units correctly)
* SphericalSurface3D: rotation, translation, frame_mapping
* read steps: Identify assemblies in a step file.
* ClosedTriangleShell3D: to_trimesh method
* PointCloud3D: add method shell_distances to compute distances from triangular mesh in PointCloud3D
* BSplineSurface3D: Now the plot method uses u and v curves
* Create .geo and .msh files (Mesh geometries with GMSH)
* RevolutionSurface3D: point3d_to_2d, point2d_to_3d, plot, rectangular_cut, from_step
* RevolutionFace3D
* WiriMixin: from points: general method for Wire3D and 2D and for Contour2D and 3D.
* Added package.xml metadata in order to be listed in the FreeCAD Addon Manager
* Edge: local_discretization
* ArcEllipse2d: point_at_abscissa, translation, split, point_distance.

### Fixed

* WireMixin: abscissa (add tolerance as parameter)
* OpenRoundedLineSegment2D: deleted discretization_points() so it uses the one from WireMixin.
* Contour2D: moved bounding_rectangle and get_bounding_rectangle to Wire2D.
* BSplineCurve: from_points_interpolation, uses centripedal method for better fitting.
* Conical, Cylindrical and Toroidal Surfaces 3D: fix face_from_contours - bug when step file doesnot follow a standard.
* BSplineSurface3D: debug linesegment2d_to_3d method.
* Parametric operations with BSpline curves.
* OpenTriangleShell3D: fix from_mesh_data method.
* PeriodicalSurface: fix face from contours.
* LineSegment2D.line_intersections: verify if colinear first.
* Cylinder: to_dict, min_distance_to_other_cylinder.
* Step_assemblies: consider when no transformation is needed.
* fix some pydocstyle errors
* Script/step/workflow: Update Workflow, use last version of dessia_common
* LineSegment3D: Rotation method update due to points attribute deletion
* ConicalSurface3D: fix from_step class method by adding the angle convertion factor
* fix f string usage
* Add some typings
* Step: Step translator now handles some EDGE_LOOP inconsistencies coming from step files
* Arc2d: point_belongs, abscissa.


### Removed

- edges: remove attributes points from lines & linesegments for performance purpose


### Performance improvements

- wires.py's 2D objects: chache bounding_rectangle results
- faces.py's Triangle3D objects: subdescription points and triangles
- EdgeCollection3D: new object for displaying series of edges
- BSplineSurface3D: compile BSplineSurface3D.derivatives
- Contour2D.area(): save area in a cache variable.
- Contour2D.__eq__(): verify contour length first, when verify if two contours are the same.
- Contour2D.is_inside(): verify first if the area of the contour2 is not smaller that contour 1.
- Disabling pointer in to_dict for most primitives
- Better hash for shells, contours & wires 


### Refactorings
- Remove usage of deprecated method old_coordinates and new_coordinates
- Indicate 'inplace' methods as deprecated
* Wire: extract_with_points

### Documentation
- BoundingBox docstrings

### Unittests
- ConicalSurface3D: face_from_contours, bsplinecurve3d_to_2d.
- CompositePrimitive2D: rotation, translation, frame_mapping
- core.py: delete_double_point, step_ids_to_str
- CompositePrimitive3D: plot
- BoundingRectangle: bounds, plot, area, center, b_rectangle_intersection, is_inside_b_rectangle, point_belongs,
intersection_area, distance_to_b_rectangle, distance_to_point
- BoundingBox: center, add, to_dict, points, from_bounding_boxes, from_points, to_frame, volume, bbox_intersection,
is_inside_bbox, intersection_volume, distance_to_bbox, point_belongs, distance_to_point, plot
* VolumeModel: eq, volume, rotation, translation, frame_mapping, bounding_box, plot
* Wire: extract_with_points, split_with_two_points
* Arc2d: point_belongs, abscissa.
* ArcEllipse2d: point_belongs, abscissa, init, translation, split, point_at_abscissa, point_distance.

### CI
- add spell check to pylint with pyenchant
- make code_pydocstyle more explicit
- upload html coverage to cdn.dessia.tech
- limit time effect on master & testing

## v0.8.0 [Released 26/01/2023]

### New Features

- PlaneFace3D: project_faces
- OpenShell3D: project_coincident_faces_of
- GmshParser: to_vtk
- BSplineCurve: derivatives
- ClosedPolygon2D: point_belongs, now the user can choose whether points on the edge of the polygon
            should be considered inside or not.
- ArcEllipse2D: line_intersections, frame_mapping, linesegment_intersections
- Line2D: point_belongs, frame_mapping()
- New Class wires.Ellipse2D
- Ellipse2D: point_over_ellipse(), line_intersections(), linesegment_intersections(), discretization_points(),
abscissa(), point_angle_with_major_dir(), area(), rotation(), tranlation(), frame_mapping()
- Plane3D: is_parallel, fullarc_intersections
- Arc2D: cut_betweeen_two_points
- Contour3D: linesegment_intersections, line_intersections
- Circle3D: primitives: [Arc3D, Arc3D], get_primitives, abscissa, linesegment_intersections
- Arc3D: line_intersections, linesegment_intersections
- new module utils: intersections -> circle_3d_linesegment_intersections
- hash for Frame2D
- Ellipse3D: point_belongs, abscissa, length, to_2d
- CylindricalSurface3D: point_on_surface, is_coincident, arcellipse3d_to_2d
- BSplineSurface3D: derivatives

### Fixed

- PlaneFace3D: cut_by_coincident_face (consider self.inner_contours inside face)
- Contour2D: bounding_rectangle (specify number_points for discretization_points), point_belongs
- Line2D: line_intersections
- BSplineCurve2D: line_intersections
- PlaneFace3D: cut_by_coincident_face (consider self.inner_contours inside face)
- BSplineCurve2D: bounding_rectangle (specify number_points for discretization_points)
- Mesh: delete_duplicated_nodes
- BSplineSurface3D: fix arc3d_to_2d method
- Frame3D : fix from_point_and_vector method ( error for the case vector=main_axis)
- BSplineCurve2D: linesegment_intersections
- Contour2D: merge_primitives_with
- BSplineCurve: fix to take into account weighted B-spline curves.
- Step: fix reading of rational BSpline curves and surfaces from step file.
- BSplineCurve2D: tangent (use position/length)
- Babylon: some scene settings for better rendering
- Arc2D: fix get_center: name referenced before assignement
- SphericalSurface3D : enhancement of primitives parametrization on surface parametric domain.
- BSplineSurface3D: debug linesegment2d_to_3d method.
- Parametric operations with BSpline curves.
- OpenTriangleShell3D: fix from_mesh_data method
- pydocstyle fixes
- bounding box: fix for cylindrical and BSplineCurve3D
- contour2d: ordering_primitives, order_primitives
- Plane3D: plane_intersections, is_coindident
- contour2d: ordering_primitives, order_primitives
- Linesegment2D: infinite_primitive
- Arc2D: point_belongs
- Arc2D: infinite_primitive
- Wire2D: infinite_intersections
- infinite primitive offset of linesegment
- Ellispe3D: discretization_points
- BSplineSurface: Improved surface periodicity calculation

### Removed

- babylon script remaining functions

### Performance improvements
- ClosedPolygon2D: triangulation
- Cylinder: min_distance_to_other_cylinder
- BSplineCurve: discretization_points
- Face3D: triangulation
- triangulation performance by use of Node2D instead of points (x15 on casing)
- cache variable self._polygon_point_belongs_100, to avoid recalculating each
time we have to verify if a point is inside
- Improvements in BSplineSurface3D.point3d_to_2d performance
- Triangle3D serialization speed-up
- Serialization without memo for faces
- Custom serialization for BsplineCurves

### Refactorings

- Basis2D, Basis3D, Frame2D, Frame3D: old_coordinates and new_coordinates method are now deprecated.
local_to_global_coordinates and global_to_local_coordinates are the new more explicit ones.
- Line3D: intersections

### Unittests

- Contour2D: point_belongs
- Basis2D, Basis3D, Frame2D, Frame3D: local_to_global_coordinates and global_to_local_coordinates
- ArcEllipse2D: linesegment_intersections
- LineSegment2D: to_wire
- Line2D: point_belongs
- BSplineCurve2D: line_intersections
- Ellipse2D.point_over_ellipse()
- Ellipse2D.line_intersections()
- Ellipse2D.linesegment_intersections()
- Ellipse2D.discretization_points()
- Ellipse2D.abscissa()
- Ellipse2D.point_angle_with_major_dir()
- Ellipse2D.area()
- Ellipse2D.rotation()
- Ellipse2D.tranlation()
- Ellipse2D.frame_mapping()
- Line2D.frame_mapping()
- Plane3D: plane_intersections, fullarc_intersections, is_parallel, is_coincident
- Contour2D: offset
- ArcEllipse3D.to_2d()
- Circle3D: point_belongs
- Circle3D: discretization_points
- Arc3D: line_intersections, linesegment_intersections
- Contour2D: ordering_contour, is_ordered, order_contour
- Ellipse3D: point_belongs, abscissa, length, to_2d, discretization_points
- CylindricalSurface3D: point_on_surface, is_coincident

### CI

- Mandatory CHANGELOG.md update for PR
- pre-commit checks with cython-lint

## v0.7.0

### New Features

- Open/Closed TriangleShells: ability to implement specific algorithm to triangles
- Block: faces_center (calculate directly point in the middle of the faces)
- Circle2D: split_by_line
- BoundingRectangle: bounds, plot, area, center, b_rectangle_intersection, is_inside_b_rectangle, point_belongs, intersection_area, distance_to_b_rectangle, distance_to_point
- Cylinder: random_point_inside, interference_volume_with_other_cylinder, lhs_points_inside
- CylindricalSurface3D: line_intersections, linesegment_intersections, plane_intersection
- Line2D: point_distance
- Line3D: to_2d
- Line3D: skew_to (verifies if two Line3D are skew)
- LineSegment3D: line_interserctions
- ArcEllipse3D: discretization_points
- FullArc3D: linesegment_intersections
- Line: sort_points_along_line
- Line2D: point_belongs
- ArcEllipse2D: length, point_belongs, abscissa, bounding_rectangle, straight_line_area, discretization_points, reverse

### Fixed

- Contour2D: point_belongs
- BsplineCurve: abscissa (use different start point between 0 and length)
- Arc3D: plot
- Cylinder: point_belongs
- FullArc3D: plot (use discretization_points instead of discretise)
- Face3D: line_intersections: consider borders
- STL: from stream (use BinaryFile and StringFile instead of io.BinaryIO and FileIO)
- Step: from stream (use BinaryFile instead of io.BinaryIO)
- Contour: is_overlapping (consider intersecting_points is empty)
- LineSegment2D: to_wire (use discretization_points instead of discretise)
- ArcEllipse2D: to_3d
- Fix boolean operations when faces are 100% coincident
- Fix some to_step methods from edges.py and faces.py


### Performance improvements

- Avoid unneeded bbox computation


### Refactorings

- cleanup of ClosedShell (double methods with Openshells)
- LineSegment3D: intersections
- Line2D: sort_points_along_line



### Unittests

- PlaneFace3D: line_intersections
- BsplineCurve: abscissa
- Circle2D: split_by_line
- BoundingRectangle: area, center, intersection, is_inside, point_belongs, intersection_area, distance_to_point, distance_to_b_rectangle
- Cylinder: point_belongs, random_point_inside, interference_volume_with_other_cylinder, min_distance_to_other_cylinder, is_intersecting_other_cylinder, lhs_points_inside
- CylindricalFace3D: linesegment_intersections
- CylindricalSurface3D: line_intersections
- Line3D: line_distance
- Line3D: skew_to
- Line3D: intersections
- LineSegment3D: line_intersections
- LineSegment3D: linesegment_intersections
- Contour: is_overlapping
- LineSegment2D: line_intersections
- ArcEllipse3D: discretization_points
- FullArc3D: linesegment_intersections
- Line2D: sort_points_along_line
- Line3D: sort_points_along_line
- ArcEllipse2D: length, point_belongs, abscissa, bounding_rectangle, straight_line_area, discretization_points, reverse


## v0.6.1 [12/13/2022]

### Changes

- Import from dessia_common are now performed from dessia_common.core

### Fixed
- infinite primitive offset of linesegment

## v0.6.0 [11/7/2022]

### New Features

- Stl:load_from_file, to_volume_model
- Surface2D: copy (specific method)
- GmshParser: read_file (.msh) and related methods, define_triangular_element_mesh, define_tetrahedron_element_mesh
- Circle2D: primitives (defined with 2 Arc2D)
- Node2D/3D, TriangularElement, QuadrilateralElement2D, TriangularElement3D
- ElementsGroup: nodes, elements_per_node
- Mesh: bounding_rectangle, delete_duplicated_nodes
- PlaneFace3D: cut_by_coincident_face
- Vector2D: to_step
- BSplineCurve2D: to_step
- LineSegment3D: to_bspline_curve
- BSplineCurve3D: from_geomdl_curve
- Surface2D: line_crossings
- Surface2D: from_contour
- BSplineSurface3D: simpifly_surface - verifies if BSplineSurface3D could be a Plane3D
- OpenShell3D: to_step_face_ids
- Contour2D: repair_cut_contour
- Circle2D: cut_by_line

### Fixed

- Contour3D: average_center_point (use edge_polygon.points instead of points)
- Contour: edges_order_with_adjacent_contour
- Arc2D: translate_inplace
- Arc2D: point_belongs
- Arc2D: abscissa (consider point2d == arc2d.start/end)
- Arc2D: split (how to choose the interior point)
- Wire: extract_primitives (consider point1 and point2 belong to the same primitive, REMOVE Contour.extract_primitives)
- LineSegment: abcissa (consider point2d == arc2d.start/end)
- Contour2D: cut_by_wire
- Contour2D: point_belongs (bug when contour has only one primitive, like FullArc2D)
- Contour: contours_from_edges
- PlaneFace3D: face_intersections
- Edge: insert_knots_and_mutiplicity
- BSplineCurve3D: from_step
- Surface2D: cut_by_line
- Circle3D: to_step
- ArcEllipse3D.to_2d()
- infinite primitive offset of linesegment
- Contour3D: order_contour.

### Performance improvements

- Improve reading STEP files (Faster BSplineCurve3D.look_up_table, Better info when _edges not following eachother_ )
- Improve multiple substractions
- Speedup Contour2D.point_belongs using bounding_rectangle
- Custom to dicts for Shells and primitives inheriting


### Refactorings

- Normalize STL methods regarding STEP
- Refacor and update old code in mesh.py
- Define a Parent class 'Triangle' for Triangle2D/3D


### Unittests

- Wire: extract_primitives, extract_without_primitives


## v0.5.0

### New Features

- Contour: is_overlapping, is_supperposing
- Point, Edges and Wires: axial_symmetry
- Surface2D: rotation, rotation_inplace
- Wire2D: bsplinecurve_crossings,  bsplinecurve_intersections
- Cylinder: min_distance_to_other_cylinder, is_intersecting_other_cylinder
- New point_distance method for Wire3D

### Fixed

- Wire3D.babylonjs
- BSplineSurface3D.merge_with (consider overlapping, intersecting surfaces)
- Wire.extract_primitives (consider point1 & point2 belong to the same primitive)
- Wire.extract_without_primitives (consider the primitives’ order to choose the primitives)
- Contour.shared_primitives_with (consider contours sharing a lot of primitives groups)
- Contour2D.contour_intersections (check if the point is not already in the lis)
- Line.is_between_points (consider point1==point2)
- BSplineCurve2D.split (consider point==start/end)
- Contour3D.bounding_box (use _utd_bounding_box to be defined as a property)
- BSplineSurface3D.grid2d_deformed (add more constraints to compute surface deformation)
- BSplineSurface3D.from_cylindrical_faces (consider **kwargs parameters)
- Duplicated methods cleaned
- triangulation of planar faces
- Wire3D: fix Bounding box
- Wire3D: Bounding box
- Arc2D: primitives bad calculation (arc2d)
- Update plotdata in setup.py
- add some fixes pydocstyle

### Performance improvements

- Remove Copy param from movement of primitives and add inplace methods
- Improve union operations
- Return the same result type (a boolean) in Contour.is_sharing_primitives_with
- Add hidden attribute _bounding_rectangle for Contour2D
- Add hidden attribute _length for BSplineCurve2D/3D
- Consider different types of primitives in Wire.wire_intersections/wire_crossings
- Add hidden attribute _length for Edge

### Refactorings

- Define _eq_ in Contour (to be used for both 2D and 3D)
- Use Grid2D object in different BSplineSurface3D methods (especially: to_2d_with_dimension)
- Define length in LineSegment (to be used for both 2D and 3D)
- Delete diplicated methods (length and point_at_abscissa) from Contour3D (inherit from Wire)
- Define a Parent class 'Bsplinecurve' to mutulize Bsplinecurve2D/3D methods
- Clean duplicated methods
- Define length in LineSegment (to be used for both 2D and 3D)
- Delete diplicated methods (length and point_at_abscissa) from Contour3D (inherit from Wire)
- Define a Parent class 'Bsplinecurve' to mutulize Bsplinecurve2D/3D methods


## v0.4.0
### Fixed
- various fixes in cuts of wires and contours
- Fix of missing face in Union
- following dessia_common v0.7.0


## v0.3.0

### New Features
- Bspline with dimensions
- cut_by_line for Surface2D
- Bspline merge

### Fixed
- Various Steps improvement
- Bspline periodicity in step reading
- sewing improvements
- Substraction of shells

## v0.2.10

### New Features

- union of shells (only with planeface for the moment
- Sewing of polygon3D
- Concav hull of PointCloud2D

## v0.2.9

### New Features

- support STL import & export
- point cloud2D & cloud3D

## v0.2.8

### New Features

- support stringIO in step save

### Fixes

- depack of point2D
- to_vector2D

### Performance improvements

- better bounding box for cylindrical face


## [v0.2.7]
### Changed
- direction vector of linesegments are now normalized

### New Features

- straight line area for BsplineCurve2D
- split of circleby start end
- closedpolygon2d is_trigo
- Auto-adaptative camera/edge width babylonjs
- splitting of bsplinecurve2d
- BezierSurface3D implemented
- added rotation and translation for faces
- new classes BezierCurve2D and BezierCurve3D
- spherical surface
- (core): update plot_data method
- update plot_data methods in wires and edges
- step almost working for cylindrical, conical toroidal
- difference between intersections and crossings
- plot_data version set to 0.3.8 or above

### Fixes

- support of mixed vector point in to step
- remove debug mode babylonjs
- remove sci notation in step export
- use stable cdn for babylonjs
- sweep extrusion length
- line circle intersection with tolerance, normal and dir vector for arc
- offset of wire
- remove useless non serializable attr
- secondmoment area from straight lines
- reversed faces in extrusion correction
- enhancement of rotation/translation of shells
- bug fix BezierCurve2D and 3D
- eq and hash for basis and frames
- shell and frame mapped shell correctly read
- small try except added for step reading
- all SHAPE_REPRESENTATION are now read
- Arc3D from step full debug
- arc3d to 2d in bspline3d surface
- missing faces at end of sweep
- splitting faces and arcs
- perf in display nodes and toroidal aspect
- setup.py requires plot_data>=0.3.9
- (primitives2d): serialization
- debug of shell method
- porting shells methods
- Debug of conical faces
- Porting cylinders and hollow
- porting from missing from_contour3d for planeface
- reading steps, but artefact on faces
- Correcting arc from_step

### Performance improvements

- LineSegment2D.points is non serializable attribute
- ClosedPolygon2D.line_segment is non_serializable_attributes
- Optimization of mesh generation

#### Refactorings
- (edges): put data argument back into Arc2D.plot_data()
- (edges): redefined Arc2D.plot_data()

## v0.2.6

### Changed
- debugs on frame 2D

### Optimized
- babylon data generation speed up

## v0.2.5

### Added
- translation and rotation for various primitives

### Changed
- Frame3D rotation takes also into account origin
- following plot_data v0.5.3

## v0.2.4
### Added
- handle spherical surfaces
- positionning of parts in STEP reading

## v0.2.1
### Added
- step export

## v0.2

### Changed
- modules -2D or *3D renamed in *2d, *3d
- point and vector declared with their x, y, z vm.Point2D((0, 0)) -> vm.Point2D(0, 0)
- separating in new modules: display, wires, edges...
- PEP8: method names
- PointAtCurvilinearAbscissa changed to point_at_abscissa
- MPLPlot changed to plot()
- plot now returns only ax instead of fig, ax

## v0.1.11

### Added
- Calculate the distance between LineSegment3D/LS3D, Arc3D/LS3D, Arc3D/Arc3D and between CylindricalFace3D too.
- Use PlaneFace3D with contours2D in a classic way and use it with contours3D with a 'from_contours3d' as CylindricalFace3D does.
- Calculate the distance between CylindricalFace3D and PlaneFace3D.
- Calculate the distance between CylindricalFace3D, PlaneFace3D and ToroidalFace3D.
- contours2d.tessel_points which gives all points of a contour2d, and .points the end points of primitives.
- Implementation of ConicalFace3D in Core and RevolvedProfile.
- Implementation of SphericalFace3D in Core.
- BSplineFace3D works.

### Changed
- cut_contours in Face3D which take all points from a Contour2D, not one side like before. Furthermore, it is light and quick.

## [v0.1.10]
- typings
- workflow to instanciate point

## [v0.1.9]

### Added
- mesh module

## [v0.1.8]

### Added
- color and alpha options for various primitives
- line segments intersection

### Debug
- arcs: is_trigo and angle were sometimes false

## [v0.1.7]

### Added
- random vector and points
- dashed line option in babylon of LineSegment3D
- Measure2D
- babylon_data: a dict language to describe models to be unpacked by a babylonjs unpacker

### Removed
- constants o2D, x2D, y2D...: use O2D, X2D...

### Changed
- Mesure -> Measure3D
