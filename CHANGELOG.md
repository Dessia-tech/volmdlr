# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## Unrealeased

## v0.6.0 [11/7/2022]

### New Features

* Block: faces_center (calculate directly point in the middle of the faces)
* Circle2D: split_by_line
* BoundingRectangle: bounds, plot, area, center, b_rectangle_intersection, is_inside_b_rectangle, point_belongs, intersection_area, distance_to_b_rectangle, distance_to_point


### Fixed

* BsplineCurve: abscissa (use different start point between 0 and length)
* Arc3D: plot


### Performance improvements



### Refactorings



### Unittests

* PlaneFace3D: line_intersections
* BsplineCurve: abscissa
* Circle2D: split_by_line
* BoundingRectangle: area, center, intersection, is_inside, point_belongs, intersection_area, distance_to_point, distance_to_b_rectangle


## v0.6.0 [Unrealeased]

### New Features

* Stl:load_from_file, to_volume_model
* Surface2D: copy (specific method)
* GmshParser: read_file (.msh) and related methods, define_triangular_element_mesh, define_tetrahedron_element_mesh
* Circle2D: primitives (defined with 2 Arc2D)
* Node2D/3D, TriangularElement, QuadrilateralElement2D, TriangularElement3D
* ElementsGroup: nodes, elements_per_node
* Mesh: bounding_rectangle, delete_duplicated_nodes
* PlaneFace3D: cut_by_coincident_face


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
