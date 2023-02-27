# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.10.0 [Unreleased yet]

### New Features

* Write .msh file (with stream)
* Arc: reverse
* Circle2D: bsplinecurve_intersections, point_distance

### Fixed

### Removed

### Performance improvements
* BSplinecurve2D: point_belongs

### Refactorings

* ContourMixin: to_polygon (for both 2D and 3D)
* BSplineCurve2D.point_distance

### Unittests
*BSplineCurve2D: point_distance, point_belongs
Circle2D: bspline_intersections, point_distance

## v0.9.0 [Testing]

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
* ConicalSurface3D, CylindricalSurface3D: plot method


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

### Removed

* edges: remove attributes points from lines & linesegments for performance purpose


### Performance improvements

* wires.py's 2D objects: chache bounding_rectangle results
* faces.py's Triangle3D objects: subdescription points and triangles
* EdgeCollection3D: new object for displaying series of edges
* BSplineSurface3D: compile BSplineSurface3D.derivatives
* Contour2D.area(): save area in a cache variable.
* Contour2D.__eq__(): verify contour length first, when verify if two contours are the same.
* Contour2D.is_inside(): verify first if the area of the contour2 is not smaller that contour 1.
* Disabling pointer in to_dict for most primitives
* Better hash for shells, contours & wires 


### Refactorings
- Remove usage of deprecated method old_coordinates and new_coordinates
- Indicate 'inplace' methods as deprecated


### Documentation
- BoundingBox docstrings

### Unittests
* ConicalSurface3D: face_from_contours, bsplinecurve3d_to_2d.
* CompositePrimitive2D: rotation, translation, frame_mapping
* core.py: delete_double_point, step_ids_to_str
* CompositePrimitive3D: plot
* BoundingRectangle: bounds, plot, area, center, b_rectangle_intersection, is_inside_b_rectangle, point_belongs,
intersection_area, distance_to_b_rectangle, distance_to_point
* BoundingBox: center, add, to_dict, points, from_bounding_boxes, from_points, to_frame, volume, bbox_intersection,
is_inside_bbox, intersection_volume, distance_to_bbox, point_belongs, distance_to_point, plot
* VolumeModel: eq, volume, rotation, translation, frame_mapping, bounding_box, plot

### CI
- add spell check to pylint with pyenchant
- make code_pydocstyle more explicit
- upload html coverage to cdn.dessia.tech


## v0.8.0 [Released 26/01/2023]

### New Features

* PlaneFace3D: project_faces
* OpenShell3D: project_coincident_faces_of
* GmshParser: to_vtk
* BSplineCurve: derivatives
* ClosedPolygon2D: point_belongs, now the user can choose whether points on the edge of the polygon
            should be considered inside or not.
* ArcEllipse2D: line_intersections, frame_mapping, linesegment_intersections
* Line2D: point_belongs, frame_mapping()
* New Class wires.Ellipse2D
* Ellipse2D: point_over_ellipse(), line_intersections(), linesegment_intersections(), discretization_points(),
abscissa(), point_angle_with_major_dir(), area(), rotation(), tranlation(), frame_mapping()
* Plane3D: is_parallel, fullarc_intersections
* Arc2D: cut_betweeen_two_points
* Contour3D: linesegment_intersections, line_intersections
* Circle3D: primitives: [Arc3D, Arc3D], get_primitives, abscissa, linesegment_intersections
* Arc3D: line_intersections, linesegment_intersections
* new module utils: intersections -> circle_3d_linesegment_intersections
* hash for Frame2D
* Ellipse3D: point_belongs, abscissa, length, to_2d
* CylindricalSurface3D: point_on_surface, is_coincident, arcellipse3d_to_2d
* BSplineSurface3D: derivatives

### Fixed

* PlaneFace3D: cut_by_coincident_face (consider self.inner_contours inside face)
* Contour2D: bounding_rectangle (specify number_points for discretization_points), point_belongs
* Line2D: line_intersections
* BSplineCurve2D: line_intersections
* PlaneFace3D: cut_by_coincident_face (consider self.inner_contours inside face)
* BSplineCurve2D: bounding_rectangle (specify number_points for discretization_points)
* Mesh: delete_duplicated_nodes
* BSplineSurface3D: fix arc3d_to_2d method
* Frame3D : fix from_point_and_vector method ( error for the case vector=main_axis)
* BSplineCurve2D: linesegment_intersections
* Contour2D: merge_primitives_with
* BSplineCurve: fix to take into account weighted B-spline curves.
* Step: fix reading of rational BSpline curves and surfaces from step file.
* BSplineCurve2D: tangent (use position/length)
* Babylon: some scene settings for better rendering
* Arc2D: fix get_center: name referenced before assignement
* SphericalSurface3D : enhancement of primitives parametrization on surface parametric domain.
* BSplineSurface3D: debug linesegment2d_to_3d method.
* Parametric operations with BSpline curves.
* OpenTriangleShell3D: fix from_mesh_data method
* pydocstyle fixes
* bounding box: fix for cylindrical and BSplineCurve3D
* contour2d: ordering_primitives, order_primitives
* Plane3D: plane_intersections, is_coindident
* contour2d: ordering_primitives, order_primitives
* Linesegment2D: infinite_primitive
* Arc2D: point_belongs
* Arc2D: infinite_primitive
* Wire2D: infinite_intersections
* infinite primitive offset of linesegment
* Ellispe3D: discretization_points
* BSplineSurface: Improved surface periodicity calculation

### Removed

* babylon script remaining functions

### Performance improvements
* ClosedPolygon2D: triangulation
* Cylinder: min_distance_to_other_cylinder
* BSplineCurve: discretization_points
* Face3D: triangulation
* triangulation performance by use of Node2D instead of points (x15 on casing)
* cache variable self._polygon_point_belongs_100, to avoid recalculating each
time we have to verify if a point is inside
* Improvements in BSplineSurface3D.point3d_to_2d performance
* Triangle3D serialization speed-up
* Serialization without memo for faces
* Custom serialization for BsplineCurves

### Refactorings

* Basis2D, Basis3D, Frame2D, Frame3D: old_coordinates and new_coordinates method are now deprecated.
local_to_global_coordinates and global_to_local_coordinates are the new more explicit ones.
* Line3D: intersections

### Unittests

* Contour2D: point_belongs
* Basis2D, Basis3D, Frame2D, Frame3D: local_to_global_coordinates and global_to_local_coordinates
* ArcEllipse2D: linesegment_intersections
* LineSegment2D: to_wire
* Line2D: point_belongs
* BSplineCurve2D: line_intersections
* Ellipse2D.point_over_ellipse()
* Ellipse2D.line_intersections()
* Ellipse2D.linesegment_intersections()
* Ellipse2D.discretization_points()
* Ellipse2D.abscissa()
* Ellipse2D.point_angle_with_major_dir()
* Ellipse2D.area()
* Ellipse2D.rotation()
* Ellipse2D.tranlation()
* Ellipse2D.frame_mapping()
* Line2D.frame_mapping()
* Plane3D: plane_intersections, fullarc_intersections, is_parallel, is_coincident
* Contour2D: offset
* ArcEllipse3D.to_2d()
* Circle3D: point_belongs
* Circle3D: discretization_points
* Arc3D: line_intersections, linesegment_intersections
* Contour2D: ordering_contour, is_ordered, order_contour
* Ellipse3D: point_belongs, abscissa, length, to_2d, discretization_points
* CylindricalSurface3D: point_on_surface, is_coincident

### CI

* Mandatory CHANGELOG.md update for PR
* pre-commit checks with cython-lint

## v0.7.0 

### New Features

* Open/Closed TriangleShells: ability to implement specific algorithm to triangles
* Block: faces_center (calculate directly point in the middle of the faces)
* Circle2D: split_by_line
* BoundingRectangle: bounds, plot, area, center, b_rectangle_intersection, is_inside_b_rectangle, point_belongs, intersection_area, distance_to_b_rectangle, distance_to_point
* Cylinder: random_point_inside, interference_volume_with_other_cylinder, lhs_points_inside
* CylindricalSurface3D: line_intersections, linesegment_intersections, plane_intersection
* Line2D: point_distance
* Line3D: to_2d
* Line3D: skew_to (verifies if two Line3D are skew)
* LineSegment3D: line_interserctions
* ArcEllipse3D: discretization_points
* FullArc3D: linesegment_intersections
* Line: sort_points_along_line
* Line2D: point_belongs
* ArcEllipse2D: length, point_belongs, abscissa, bounding_rectangle, straight_line_area, discretization_points, reverse

### Fixed

* Contour2D: point_belongs
* BsplineCurve: abscissa (use different start point between 0 and length)
* Arc3D: plot
* Cylinder: point_belongs
* FullArc3D: plot (use discretization_points instead of discretise)
* Face3D: line_intersections: consider borders
* STL: from stream (use BinaryFile and StringFile instead of io.BinaryIO and FileIO)
* Step: from stream (use BinaryFile instead of io.BinaryIO)
* Contour: is_overlapping (consider intersecting_points is empty)
* LineSegment2D: to_wire (use discretization_points instead of discretise)
* ArcEllipse2D: to_3d
* Fix boolean operations when faces are 100% coincident
* Fix some to_step methods from edges.py and faces.py


### Performance improvements

* Avoid unneeded bbox computation


### Refactorings

* cleanup of ClosedShell (double methods with Openshells)
* LineSegment3D: intersections
* Line2D: sort_points_along_line



### Unittests

* PlaneFace3D: line_intersections
* BsplineCurve: abscissa
* Circle2D: split_by_line
* BoundingRectangle: area, center, intersection, is_inside, point_belongs, intersection_area, distance_to_point, distance_to_b_rectangle
* Cylinder: point_belongs, random_point_inside, interference_volume_with_other_cylinder, min_distance_to_other_cylinder, is_intersecting_other_cylinder, lhs_points_inside
* CylindricalFace3D: linesegment_intersections
* CylindricalSurface3D: line_intersections
* Line3D: line_distance
* Line3D: skew_to
* Line3D: intersections
* LineSegment3D: line_intersections
* LineSegment3D: linesegment_intersections
* Contour: is_overlapping
* LineSegment2D: line_intersections
* ArcEllipse3D: discretization_points
* FullArc3D: linesegment_intersections
* Line2D: sort_points_along_line
* Line3D: sort_points_along_line
* ArcEllipse2D: length, point_belongs, abscissa, bounding_rectangle, straight_line_area, discretization_points, reverse


## v0.6.1 [12/13/2022]

### Changes

* Import from dessia_common are now performed from dessia_common.core

### Fixed
* infinite primitive offset of linesegment

## v0.6.0 [11/7/2022]

### New Features

* Stl:load_from_file, to_volume_model
* Surface2D: copy (specific method)
* GmshParser: read_file (.msh) and related methods, define_triangular_element_mesh, define_tetrahedron_element_mesh
* Circle2D: primitives (defined with 2 Arc2D)
* Node2D/3D, TriangularElement, QuadrilateralElement2D, TriangularElement3D
* ElementsGroup: nodes, elements_per_node
* Mesh: bounding_rectangle, delete_duplicated_nodes
* PlaneFace3D: cut_by_coincident_face
* Vector2D: to_step
* BSplineCurve2D: to_step
* LineSegment3D: to_bspline_curve
* BSplineCurve3D: from_geomdl_curve
* Surface2D: line_crossings
* Surface2D: from_contour
* BSplineSurface3D: simpifly_surface - verifies if BSplineSurface3D could be a Plane3D
* OpenShell3D: to_step_face_ids
* Contour2D: repair_cut_contour
* Circle2D: cut_by_line

### Fixed

* Contour3D: average_center_point (use edge_polygon.points instead of points)
* Contour: edges_order_with_adjacent_contour
* Arc2D: translate_inplace
* Arc2D: point_belongs
* Arc2D: abscissa (consider point2d == arc2d.start/end)
* Arc2D: split (how to choose the interior point)
* Wire: extract_primitives (consider point1 and point2 belong to the same primitive, REMOVE Contour.extract_primitives)
* LineSegment: abcissa (consider point2d == arc2d.start/end)
* Contour2D: cut_by_wire
* Contour2D: point_belongs (bug when contour has only one primitive, like FullArc2D)
* Contour: contours_from_edges
* PlaneFace3D: face_intersections
* Edge: insert_knots_and_mutiplicity
* BSplineCurve3D: from_step
* Surface2D: cut_by_line
* Circle3D: to_step
* ArcEllipse3D.to_2d()
* infinite primitive offset of linesegment

### Performance improvements

* Improve reading STEP files (Faster BSplineCurve3D.look_up_table, Better info when _edges not following eachother_ )
* Improve multiple substractions
* Speedup Contour2D.point_belongs using bounding_rectangle
* Custom to dicts for Shells and primitives inheriting


### Refactorings

* Normalize STL methods regarding STEP
* Refacor and update old code in mesh.py
* Define a Parent class 'Triangle' for Triangle2D/3D


### Unittests

* Wire: extract_primitives, extract_without_primitives


## v0.5.0

### New Features

* Contour: is_overlapping, is_supperposing
* Point, Edges and Wires: axial_symmetry
* Surface2D: rotation, rotation_inplace
* Wire2D: bsplinecurve_crossings,  bsplinecurve_intersections
* Cylinder: min_distance_to_other_cylinder, is_intersecting_other_cylinder
* New point_distance method for Wire3D

### Fixed

* Wire3D.babylonjs
* BSplineSurface3D.merge_with (consider overlapping, intersecting surfaces)
* Wire.extract_primitives (consider point1 & point2 belong to the same primitive)
* Wire.extract_without_primitives (consider the primitives’ order to choose the primitives)
* Contour.shared_primitives_with (consider contours sharing a lot of primitives groups)
* Contour2D.contour_intersections (check if the point is not already in the lis)
* Line.is_between_points (consider point1==point2)
* BSplineCurve2D.split (consider point==start/end)
* Contour3D.bounding_box (use _utd_bounding_box to be defined as a property)
* BSplineSurface3D.grid2d_deformed (add more constraints to compute surface deformation)
* BSplineSurface3D.from_cylindrical_faces (consider **kwargs parameters)
* Duplicated methods cleaned
* triangulation of planar faces
* Wire3D: fix Bounding box
* Wire3D: Bounding box
* Arc2D: primitives bad calculation (arc2d)
* Update plotdata in setup.py
* add some fixes pydocstyle

### Performance improvements

* Remove Copy param from movement of primitives and add inplace methods
* Improve union operations
* Return the same result type (a boolean) in Contour.is_sharing_primitives_with
* Add hidden attribute _bounding_rectangle for Contour2D
* Add hidden attribute _length for BSplineCurve2D/3D
* Consider different types of primitives in Wire.wire_intersections/wire_crossings
* Add hidden attribute _length for Edge

### Refactorings

* Define _eq_ in Contour (to be used for both 2D and 3D)
* Use Grid2D object in different BSplineSurface3D methods (especially: to_2d_with_dimension)
* Define length in LineSegment (to be used for both 2D and 3D)
* Delete diplicated methods (length and point_at_abscissa) from Contour3D (inherit from Wire)
* Define a Parent class 'Bsplinecurve' to mutulize Bsplinecurve2D/3D methods
* Clean duplicated methods
* Define length in LineSegment (to be used for both 2D and 3D)
* Delete diplicated methods (length and point_at_abscissa) from Contour3D (inherit from Wire)
* Define a Parent class 'Bsplinecurve' to mutulize Bsplinecurve2D/3D methods


## v0.4.0
### Fixed
* various fixes in cuts of wires and contours
* Fix of missing face in Union
* following dessia_common v0.7.0


## v0.3.0

### New Features
* Bspline with dimensions
* cut_by_line for Surface2D
* Bspline merge

### Fixed
* Various Steps improvement
* Bspline periodicity in step reading
* sewing improvements
* Substraction of shells

## v0.2.10

### New Features

* union of shells (only with planeface for the moment
* Sewing of polygon3D
* Concav hull of PointCloud2D

## v0.2.9

### New Features

* support STL import & export
* point cloud2D & cloud3D

## v0.2.8

### New Features

* support stringIO in step save

### Fixes

* depack of point2D
* to_vector2D

### Performance improvements

* better bounding box for cylindrical face


## [v0.2.7]
### Changed
* direction vector of linesegments are now normalized

### New Features

* straight line area for BsplineCurve2D
* split of circleby start end
* closedpolygon2d is_trigo
* Auto-adaptative camera/edge width babylonjs
* splitting of bsplinecurve2d
* BezierSurface3D implemented
* added rotation and translation for faces
* new classes BezierCurve2D and BezierCurve3D
* spherical surface
* (core): update plot_data method
* update plot_data methods in wires and edges
* step almost working for cylindrical, conical toroidal
* difference between intersections and crossings
* plot_data version set to 0.3.8 or above

### Fixes

* support of mixed vector point in to step
* remove debug mode babylonjs
* remove sci notation in step export
* use stable cdn for babylonjs
* sweep extrusion length
* line circle intersection with tolerance, normal and dir vector for arc
* offset of wire
* remove useless non serializable attr
* secondmoment area from straight lines
* reversed faces in extrusion correction
* enhancement of rotation/translation of shells
* bug fix BezierCurve2D and 3D
* eq and hash for basis and frames
* shell and frame mapped shell correctly read
* small try except added for step reading
* all SHAPE_REPRESENTATION are now read
* Arc3D from step full debug
* arc3d to 2d in bspline3d surface
* missing faces at end of sweep
* splitting faces and arcs
* perf in display nodes and toroidal aspect
* setup.py requires plot_data>=0.3.9
* (primitives2d): serialization
* debug of shell method
* porting shells methods
* Debug of conical faces
* Porting cylinders and hollow
* porting from missing from_contour3d for planeface
* reading steps, but artefact on faces
* Correcting arc from_step

### Performance improvements

* LineSegment2D.points is non serializable attribute
* ClosedPolygon2D.line_segment is non_serializable_attributes
* Optimization of mesh generation

#### Refactorings
* (edges): put data argument back into Arc2D.plot_data()
* (edges): redefined Arc2D.plot_data()

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
- modules *2D or *3D renamed in *2d, *3d
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
